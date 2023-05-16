## code to prepare `exampleData` dataset goes here

urlfile = "https://raw.githubusercontent.com/CnrLwlss/Ahmed_2022/master/rawdat_a_filtered.csv"

rawData = readr::read_delim(url(urlfile), delim="\t")

channels = c("raw_CIV", "raw_porin", "raw_CI")

reducedData = rawData[,c("caseno", "Fibre", channels, "controls")]

exampleData = reducedData %>% tidyr::pivot_longer(!c("caseno", "Fibre", "controls"), names_to="channel")

exampleData = as.data.frame(exampleData)

colnames(exampleData) = c("sampleID", "fibreID", "control", "channel", "value")

exampleData = plyr::match_df(exampleData, exampleData[exampleData$value>0, ])

usethis::use_data(exampleData, overwrite = TRUE, name=exampleData)

