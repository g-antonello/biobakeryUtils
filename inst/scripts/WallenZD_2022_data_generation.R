# WallenZD_2022 data generation

#' this dataset comes from the Wallen Nat. Comms paper in 2022 
#' (https://www.nature.com/articles/s41467-022-34667-x)
#' 
#' The processed source data for the analysis is published on Zenodo
#' (https://zenodo.org/records/7246185/), in a .xlsx file called 
#' `Source_Data_24Oct2022.xlsx`
#' 
#' The following code shows the steps needed to generate the source data
local_fileName <- file.path(tempdir(), "wallenZD_zenodo.xlsx")
download.file("https://zenodo.org/records/7246185/files/Source_Data_24Oct2022.xlsx?download=1", destfile = local_fileName)
# you can read sheet names like follows:
openxlsx::getSheetNames(local_fileName)

#' 01 - generation of `WallenZD_2022_metaphlan3_profiles.tsv`
#' md5sum: b76a7852a67fbb74d3b54ade931aa0a2
# UNKNOWN row gets lost with this method
mpa_profles <- openxlsx::read.xlsx(local_fileName, sheet = "metaphlan_rel_ab") 
readr::write_tsv(mpa_profles, "inst/extdata/WallenZD_2022_metaphlan3_profiles.tsv")
#' in the folder, run `bzip2 -z WallenZD_2022_metaphlan3_profiles.tsv`

#' 02 - generation of `WallenZD_2022_subjMetadata.tsv.bz2`
#' md5sum: 4f702d810a93785c1a61f8b6ccb433f3
subj_metadata <- openxlsx::read.xlsx(local_fileName, sheet = "subject_metadata") 
readr::write_tsv(subj_metadata, "inst/extdata/WallenZD_2022_subjMetadata.tsv")
#' in the folder, run `bzip2 -z WallenZD_2022_subjMetadata.tsv`