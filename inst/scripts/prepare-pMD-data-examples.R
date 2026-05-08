library(parkinsonsMetagenomicData)
data("sampleMetadata", package = "parkinsonsMetagenomicData")

data_types <- c("relative_abundance",
"pathabundance_unstratified", "genefamilies_unstratified")

example_metadata.df <- sampleMetadata[grepl("Bedarf",
sampleMetadata$study_name),]

# MetaPhlAn
Bedarf_pMD_raw_MetaPhlAn.tse <- returnSamples(example_metadata.df, data_type = data_types[1])

# HUMAnN pathways
Bedarf_pMD_raw_HUMAnN_pwy.tse <- returnSamples(example_metadata.df, data_type = data_types[2])

# save objects
usethis::use_data(Bedarf_pMD_raw_MetaPhlAn.tse, internal = FALSE)

usethis::use_data(Bedarf_pMD_raw_HUMAnN_pwy.tse, internal = FALSE)
