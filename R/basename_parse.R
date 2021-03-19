###### Known issue: 
## 1. Some sites (US-Ha1) provide top-layer as gap-filled, those variables are stripped out in the initial parsing
## 2. Fix issue parsing FETCH_XX_H_V_R & VAR_#_N & VAR_#_SD
## 3. Fix issue parsing _H_V_A_SD & _H_V_A_N cases
## 4. create a basename_parse_test_list for a list of possible variable names for testing

basename_parse <- function(var.name,                             ## list of names to be parsed
                            FP.ls,                                ## reference FP variable names
                            warning.only = T,                     ## return only warning
                            echo = F,
                            gapfill_postfix = "_PI_F"             ## can be specified as _F for QC Combined files
  ) {
                            
  ## used to parse numbers within string
  Numextract <- function(string) {
    unlist(regmatches(string, gregexpr(
      "[[:digit:]]+\\.*[[:digit:]]*", string
    )))
  }
  
  basename_decode <- data.frame(variable_names = var.name,  # original variable name 
                                working_names = NA,         # working variable name, dropped in output
                                basename = NA,              # associated basename, w/o qualifier
                                proposed_name = NA,         # proposed var name after network-aggregation 
                                qualifier_gf = NA,          # qualifier associated with gap-filling   
                                qualifier_pi = NA,          # qualifier associated with PI version, excluding gap-filling   
                                qualifier_pos = NA,         # qualifier associated with position    
                                qualifier_ag = NA,          # qualifier associated with layer-aggregation, e.g., _N, _SD
                                need_aggregate = NA,        # does var need further aggregation 
                                aggregate_method = NA,      # aggregation method
                                                            #   1. equal: unique variable, no aggregation needed, 
                                                            #             -> use base name 
                                                            #             e.g., ignore any _1_1_1 or _1 if provided
                                                            #   2. naive: already layer-aggregated & not unique, 
                                                            #             -> use provided layer index
                                                            #   3. 1step: already replicate-averaged & not unique, 
                                                            #             -> need to aggregate to layer
                                                            #   4. 2step: quadruplet provided & not unique, 
                                                            #             -> need 1 or 2 steps aggregation 
                                                            #   5. manual: mix of aggregation levels
                                                            #             -> need inspection before aggregation
                                layer = NA,                 # layer index provided, if any
                                H_index = NA,               # H index provided, if any 
                                V_index = NA,               # V index provided, if any
                                R_index = NA,               # R index provided, if any
                                is_correct_basename = NA,   # is the parsed basename recognized in FP-Standard 
                                is_unique = NA,             # is this a unique (no other position-siblings)
                                is_PI_provide = NA,         # is this a PI provided variable e.g., _PI 
                                is_gapfill = NA,            # is this a gap-filled variable, _PF_F or _F
                                is_fetch = NA,              # is this a fetch quantile variable, e.g., FETCH_70...FETCH_90
                                                            #    i.e., exceptional case where _<number> exist 
                                                            #          other than position qualifier  
                                is_layer_integrate = NA,    # is this a layer-integrated var, i.e., _#
                                is_layer_SD = NA,           # is this a standard deviation of layer-integrated var, i.e., spatial variability
                                is_layer_number = NA,       # is this a number of samples of layer-integrated var, i.e., spatial variability
                                is_replicate_average = NA,  # is this a replicate-averaged var, i.e., _<H>_<V>_A
                                is_replicate_SD = NA,       # is this a standard deviation of replicate-averaged var, i.e., _<H>_<V>_A_SD
                                is_replicate_number = NA,   # is this a number of samples of replicate-averaged var, i.e., _<H>_<V>_A_N
                                is_quadruplet = NA          # is this a quadruplet, i.e., _<H>_<V>_<R>
                                )         
                              
  
  ## find gap-filled ones
  basename_decode$is_gapfill <- grepl(gapfill_postfix,
                                      basename_decode$variable_names,
                                      perl = TRUE)
  basename_decode$working_names <- sub(gapfill_postfix,
                                       "",
                                       basename_decode$variable_names)
  basename_decode$qualifier_gf <- ifelse(basename_decode$is_gapfill,
                                         gapfill_postfix,
                                         NA)
  
  ## take out _PI 
  basename_decode$is_PI_provide <- grepl("_PI",
                                         basename_decode$working_names,
                                         perl = TRUE)
  basename_decode$working_names <- sub("_PI",
                                       "",
                                       basename_decode$working_names)
  basename_decode$qualifier_pi <- ifelse(basename_decode$is_PI_provide,
                                         "_PI",
                                         NA)
  
  ## pre-screen for FETCH_70, _80, _90
  basename_decode$is_fetch <- grepl("FETCH_[[:digit:]]+",
                                    basename_decode$working_names,
                                    perl = TRUE)
  ## [code note]: currently, no case with FETCH_XX_H_V_R, ignore this exception for now
  
  ## find quadruplet
  basename_decode$is_quadruplet <- (!basename_decode$is_fetch & 
                                      grepl("_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+",
                                            basename_decode$working_names,
                                            perl = TRUE)) | grepl(
                                              "FETCH_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+",
                                              basename_decode$working_names,
                                              perl = TRUE)
  
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_quadruplet,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE) + ifelse(basename_decode$is_fetch, 3, 0),
                                                 stop = nchar(basename_decode$working_names)),
                                          basename_decode$qualifier_pos)
                                          

  for (i1 in 1:nrow(basename_decode)) {
    if (basename_decode$is_quadruplet[i1]) {
      basename_decode[i1, c("H_index", "V_index", "R_index")] <-
        Numextract(basename_decode$qualifier_pos[i1])
    }
  }
  
  ## find replicate aggregated SD
  basename_decode$is_replicate_SD <- grepl("_[[:digit:]]+_[[:digit:]]+_A_SD",
                                       basename_decode$working_names,
                                       perl = TRUE)
  
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_replicate_SD,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+_[[:digit:]]+_A_SD",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE),
                                                 stop = regexpr("_SD",
                                                                basename_decode$working_names,
                                                                perl = TRUE) - 1),
                                          basename_decode$qualifier_pos)
  
  for (i2 in 1:nrow(basename_decode)) {
    if (basename_decode$is_replicate_SD[i2]) {
      basename_decode[i2, c("H_index", "V_index", "R_index")] <-
        c(Numextract(basename_decode$qualifier_pos[i2]), "A")
      basename_decode[i2, c("qualifier_ag")] <- c("_SD")
    }
  }
  
  ## find replicate aggregated N
  basename_decode$is_replicate_number <- grepl("_[[:digit:]]+_[[:digit:]]+_A_N",
                                           basename_decode$working_names,
                                           perl = TRUE)
  
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_replicate_number,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+_[[:digit:]]+_A_N",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE),
                                                 stop = regexpr("_N",
                                                                basename_decode$working_names,
                                                                perl = TRUE) - 1),
                                          basename_decode$qualifier_pos)
  
  for (i3 in 1:nrow(basename_decode)) {
    if (basename_decode$is_replicate_number[i3]) {
      basename_decode[i3, c("H_index", "V_index", "R_index")] <-
        c(Numextract(basename_decode$qualifier_pos[i3]), "A")
      basename_decode[i3, c("qualifier_ag")] <- c("_N")
    }
  }
  
  ## find replicate averaged
  basename_decode$is_replicate_average <- (grepl("_[[:digit:]]+_[[:digit:]]+_A",
                                                 basename_decode$working_names,
                                                 perl = TRUE) & 
                                             !basename_decode$is_replicate_number &
                                             !basename_decode$is_replicate_SD)
    
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_replicate_average,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+_[[:digit:]]+_A",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE),
                                                 stop = nchar(basename_decode$working_names)),
                                          basename_decode$qualifier_pos)
  for (i4 in 1:nrow(basename_decode)) {
    if (basename_decode$is_replicate_average[i4]) {
      basename_decode[i4, c("H_index", "V_index", "R_index")] <-
        c(Numextract(basename_decode$qualifier_pos[i4]), "A")
    }
  }
  
  ## find layer aggregated SD
  basename_decode$is_layer_SD <- grepl("_[[:digit:]]+_SD",
                                       basename_decode$working_names,
                                       perl = TRUE)
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_layer_SD,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+_SD",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE),
                                                 stop = regexpr("_SD",
                                                                basename_decode$working_names,
                                                                perl = TRUE) - 1),
                                          basename_decode$qualifier_pos)
  
  for (i5 in 1:nrow(basename_decode)) {
    if (basename_decode$is_layer_SD[i5]) {
      basename_decode[i5, c("layer")] <-
        c(Numextract(basename_decode$qualifier_pos[i5]))
      basename_decode[i5, c("qualifier_ag")] <- c("_SD")
    }
  }
  
  ## find layer aggregated Number
  basename_decode$is_layer_number <- grepl("_[[:digit:]]+_N",
                                           basename_decode$working_names,
                                           perl = TRUE)
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_layer_number,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+_N",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE),
                                                 stop = regexpr("_N",
                                                                basename_decode$working_names,
                                                                perl = TRUE) - 1),
                                          basename_decode$qualifier_pos)
  
  for (i6 in 1:nrow(basename_decode)) {
    if (basename_decode$is_layer_number[i6]) {
      basename_decode[i6, c("layer")] <-
        c(Numextract(basename_decode$qualifier_pos[i6]))
      basename_decode[i6, c("qualifier_ag")] <- c("_N")
    }
  }
  
  ## find layer aggregated variables
  basename_decode$is_layer_integrate <- (grepl("_[[:digit:]]+",
                                               basename_decode$working_names,
                                               perl = TRUE) &
                                           !basename_decode$is_fetch &
                                           !basename_decode$is_quadruplet &
                                           !basename_decode$is_replicate_average &
                                           !basename_decode$is_replicate_number &
                                           !basename_decode$is_replicate_SD &
                                           !basename_decode$is_layer_SD &
                                           !basename_decode$is_layer_number
                                         ) | (grepl("FETCH_[[:digit:]]+_[[:digit:]]+",
                                                    basename_decode$working_names,
                                                    perl = TRUE) & 
                                                basename_decode$is_fetch &
                                                !basename_decode$is_quadruplet &
                                                !basename_decode$is_replicate_average &
                                                !basename_decode$is_replicate_number &
                                                !basename_decode$is_replicate_SD &
                                                !basename_decode$is_layer_SD &
                                                !basename_decode$is_layer_number)
  
  basename_decode$qualifier_pos <- ifelse(basename_decode$is_layer_integrate,
                                          substr(basename_decode$working_names,
                                                 start = regexpr("_[[:digit:]]+",
                                                                 basename_decode$working_names,
                                                                 perl = TRUE) + ifelse(basename_decode$is_fetch, 3, 0),
                                                 stop = nchar(basename_decode$working_names)),
                                          basename_decode$qualifier_pos)
  
  for (i7 in 1:nrow(basename_decode)) {
    if (basename_decode$is_layer_integrate[i7]) {
      basename_decode[i7, c("layer")] <-
        c(Numextract(basename_decode$qualifier_pos[i7]))
    }
  }
  
  ## parse basename, w/o all qualifiers
  basename_decode$basename <- basename_decode$working_names
  for (i8 in 1:nrow(basename_decode)) {
    if (!is.na(basename_decode$qualifier_pos[i8])) {
      basename_decode$basename[i8] <- sub(basename_decode$qualifier_pos[i8],
                                          "",
                                          basename_decode$working_names[i8])
      }
    if (!is.na(basename_decode$qualifier_ag[i8])) {
      basename_decode$basename[i8] <- sub(paste0(basename_decode$qualifier_pos[i8],
                                                 basename_decode$qualifier_ag[i8]),
                                          "",
                                          basename_decode$working_names[i8])
      }
  }
  
  ## check if parsed basename valid in FP
  #   return waning if any unrecognized basename
  for (i9 in 1:nrow(basename_decode)) {
    basename_decode$is_correct_basename[i9] <- ifelse(length(which(
      FP.ls == paste(basename_decode$basename[i9])
      )) == 1,
      TRUE, FALSE)
  }
  
  if (echo & sum(!basename_decode$is_correct_basename) > 0) {
    print("[Warning] Invalid basename parsing")
    print(basename_decode$basename[!basename_decode$is_correct_basename])
  }
  
  ########################################################################################################################################
  ##### the following part works separately on var in gap-filled and non-filled groups
  ##     check if each variable is unique (after ripping out position qualifier)
  ##     for non-filled one, work to decode the aggregation level/info, and propose potential aggregation method & var name
  
  ### Work on gap-filled variables,
  basename_decode_gf <- basename_decode[basename_decode$is_gapfill, ]
  
  #   Specify if unique variable
  if (nrow(basename_decode_gf) > 0) {
    ### Specify unique variable
    basename_decode_gf$is_unique <- TRUE
    dupl.ls <-
      basename_decode_gf$basename[duplicated(basename_decode_gf$basename)]
    
    for (i6 in 1:length(dupl.ls)) {
      basename_decode_gf$is_unique[which(basename_decode_gf$basename == paste(dupl.ls[i6]))] <-
        FALSE
    }
  }
  
  ### Work on non-filled variables,
  basename_decode_nogf <-
    basename_decode[!basename_decode$is_gapfill, ]
  
  #   Specify unique variable
  if (nrow(basename_decode_nogf) > 0) {
    basename_decode_nogf$is_unique <- TRUE
    dupl.ls2 <-
      basename_decode_nogf$basename[duplicated(basename_decode_nogf$basename)]
    
    for (i10 in 1:length(dupl.ls2)) {
      basename_decode_nogf$is_unique[which(basename_decode_nogf$basename == paste(dupl.ls2[i10]))] <-
        FALSE
    }
    
    # suggest potential aggregation, based on if variable unique
    #    e.g., ignore any _1_1_1 or _1 if provided
    basename_decode_nogf$need_aggregate <-
      ifelse(basename_decode_nogf$is_unique, FALSE, TRUE)
    basename_decode_nogf$aggregate_method <-
      ifelse(basename_decode_nogf$is_unique, "equal", NA)
    basename_decode_nogf$proposed_name <-
      ifelse(basename_decode_nogf$is_unique,
             basename_decode_nogf$basename,
             NA)
    
    ## for non-unique var, check aggregation level/info (i.e., layer + duplicate + quadruplet)
    chk.var.ls <-
      unique(basename_decode_nogf$basename[!basename_decode_nogf$is_unique])
    
    for (i11 in 1:length(chk.var.ls)) {
      chk.var.loc <-
        which(basename_decode_nogf$basename == paste(chk.var.ls[i11]))
      
      ## identify layer-aggregated var
      if (sum(basename_decode_nogf$is_layer_integrate[basename_decode_nogf$basename ==
                                                      paste(chk.var.ls[i11])]) ==
          sum(basename_decode_nogf$basename == paste(chk.var.ls[i11]))) {
        if (!warning.only &
            echo)
          print(paste("[Info]", chk.var.ls[i11], "are multi-layered"))
        
        basename_decode_nogf$need_aggregate[chk.var.loc] <- FALSE
        basename_decode_nogf$aggregate_method[chk.var.loc] <- "naive"
        basename_decode_nogf$proposed_name[chk.var.loc] <-
          paste(basename_decode_nogf$basename[chk.var.loc],
                basename_decode_nogf$qualifier_pos[chk.var.loc],
                sep = "")
        
        ## identify replicate-avergaged var
      } else if (sum(basename_decode_nogf$is_replicate_average[basename_decode_nogf$basename ==
                                                               paste(chk.var.ls[i11])]) ==
                 sum(basename_decode_nogf$basename == paste(chk.var.ls[i11]))) {
        if (!warning.only &
            echo)
          print(paste("[Info]", chk.var.ls[i11], "are replicate-agregated"))
        
        basename_decode_nogf$need_aggregate[chk.var.loc] <- TRUE
        basename_decode_nogf$aggregate_method[chk.var.loc] <- "1step"
        
        ## identify quadruplet var
      } else if (sum(basename_decode_nogf$is_quadruplet[basename_decode_nogf$basename ==
                                                        paste(chk.var.ls[i11])]) ==
                 sum(basename_decode_nogf$basename == paste(chk.var.ls[i11]))) {
        if (!warning.only &
            echo)
          print(paste("[Info]", chk.var.ls[i11], "are quadruplet-specific"))
        
        basename_decode_nogf$need_aggregate[chk.var.loc] <- TRUE
        basename_decode_nogf$aggregate_method[chk.var.loc] <- "2step"
        
        ## all others as mixed, complicated cases
      } else{
        if (echo)
          print(paste("[Warning]", chk.var.ls[i11], "are mixed aggregated"))
        #print(basename_decode_nogf$variable_names[basename_decode_nogf$basename==paste(chk.var.ls[i11])])
        
        basename_decode_nogf$need_aggregate[chk.var.loc] <- TRUE
        basename_decode_nogf$aggregate_method[chk.var.loc] <- "manual"
      }
    }
  }
  
  ## prepare output
  basename_decode <- rbind.data.frame(basename_decode_nogf,
                                      basename_decode_gf)
  basename_decode <-
    basename_decode[, -which(colnames(basename_decode) == "working_names")]
  
  ## re-order to follow input order
  basename_decode <- merge.data.frame(x = data.frame(variable_names = var.name,
                                                       stringsAsFactors = FALSE),
                                      y = basename_decode,
                                      by = "variable_names",
                                      sort = FALSE)
  
  return(basename_decode)  
}



