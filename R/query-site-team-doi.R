rm(list = ls())

library(amerifluxr)

path <-
  "D:\\Housen\\Flux\\Data-exploring\\00_Synthesis_AmeriFlux_Time_Series_Cluster\\"
badm.path <- "D:\\AmeriFlux-Data\\00_olaf_data_ameriflux\\BADM\\"
path.out <- paste0(path, "miscellaneous\\")

## complete site list used
site.ls <-
  read.csv(
    paste0(path.out, "HarmonicParameter_20220701.csv"),
    header = T,
    stringsAsFactors = F
  )
na.sum <- function(x) sum(!is.na(x))
site.ls <- site.ls[which(apply(site.ls[, c(2:9)], 1, na.sum) > 0), ]

## BADM
data.in <- amerifluxr::amf_read_bif(paste0(badm.path,
                                           "AMF_AA-Flx_BIF_LEGACY_20230331.xlsx"))

###################################################################################################
## Work on Site team contact list
sgi_member <- amerifluxr::amf_extract_badm(data.in,
                                           select_group = "GRP_TEAM_MEMBER")
  
sgi_member <- sgi_member[sgi_member$SITE_ID %in% site.ls$Site & 
                           sgi_member$TEAM_MEMBER_ROLE == "PI", ]

## work on duplicate people on the list
alias.ls <-
  list(
    c(
      "Hank A. Margolis",
      "J. William Munger",
      "Russell Scott",
      "Andrew T. Black",
      "Walter Oechel",
      "A. Christopher Oishi",
      "Carl Bernacchi"
    ),
    c(
      "Hank Margolis",
      "J. William  Munger",
      "Russ Scott",
      "T. Andrew Black",
      "Walt Oechel",
      "Chris Oishi",
      "Carl J Bernacchi"
    )
  )

for (i in 1:length(alias.ls[[1]])) {
  sgi_member$TEAM_MEMBER_NAME[which(sgi_member$TEAM_MEMBER_NAME == alias.ls[[2]][i])] <-
    alias.ls[[1]][i]
}

fpt.member.ls <-
  data.frame(
    TEAM_MEMBER_NAME = names(table(sgi_member$TEAM_MEMBER_NAME)),
    TEAM_MEMBER_INSTITUTION = NA,
    TEAM_MEMBER_EMAIL = NA,
    SITE = NA,
    stringsAsFactors = F
  )
for (j in 1:nrow(fpt.member.ls)) {
  fpt.member.ls$TEAM_MEMBER_INSTITUTION[j] <-
    sgi_member$TEAM_MEMBER_INSTITUTION[which(sgi_member$TEAM_MEMBER_NAME ==
                                               fpt.member.ls$TEAM_MEMBER_NAME[j])][1]
  fpt.member.ls$TEAM_MEMBER_EMAIL[j] <-
    sgi_member$TEAM_MEMBER_EMAIL[which(sgi_member$TEAM_MEMBER_NAME == fpt.member.ls$TEAM_MEMBER_NAME[j])][1]
  fpt.member.ls$SITE[j] <-
    paste(sgi_member$SITE_ID[which(sgi_member$TEAM_MEMBER_NAME == fpt.member.ls$TEAM_MEMBER_NAME[j])], collapse =
            ";")
}

write.csv(
  fpt.member.ls,
  paste(
    path.out,
    "AMF_BADM_GRP_TEAM_MEMBER_coherence-20220701.csv",
    sep = ""
  ),
  row.names = F
)

#################################################################################
## Work on Data DOI
sgi_doi <- amerifluxr::amf_extract_badm(data.in,
                                        select_group = "GRP_DOI")

sgi_doi <- sgi_doi[sgi_doi$SITE_ID %in% site.ls$Site, ]

write.csv(
  sgi_doi,
  paste(
    path.out,
    "AMF_BADM_GRP_DOI_coherence-20220701.csv",
    sep = ""
  ),
  row.names = F
)

##################
write.csv(
  site.ls,
  paste(
    path.out,
    "TableS1_site_list_harmonic.csv",
    sep = ""
  ),
  row.names = F
)