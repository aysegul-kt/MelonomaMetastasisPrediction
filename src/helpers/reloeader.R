dfMerged <- ReadFromLogFile("Main","MergedData.csv","D")
NormDfMerged <- read_csv("output/intermediateData/Main_FinalData.csv")
ttest.results <- ReadFromLogFile("Main","ttest.results.csv","L")
ttest.results_detailed <- ReadFromLogFile("Main","ttest.results_detailed.csv","L")
miRTarBase_MTI <- read_csv("data/raw/miRTarBase_MTI.csv")
hgnc_item <-   ReadFromLogFile("Main","hgnc_item.csv","L")
joind_map <-   ReadFromLogFile("Main","joind_map.csv","L")

dfMerged <- dfMerge[-c(1)]
NormDfMerged <- NormDfMerged[-c(1)]
ttest.results <- ttest.results[-c(1)]
ttest.results_detailed <- ttest.results_detailed[-c(1)]
miRTarBase_MTI <- miRTarBase_MTI[-c(1)]
hgnc_item <- hgnc_item[-c(1)]
joind_map <- joind_map[-c(1)]

# Methylation.data <-  ReadFromLogFile("MethylationService","BetaValuesData.csv","D")

rm(df, df.mirna, df.mrna, Final.mirna_list, Final.mrna_list,merged_mirna_significancy_results, merged_mrna_significancy_results, mirna_significancyTestResultst, mrna_significancyTestResultst,
   mrna_significancyTestResultst_clone, part1, part2, sample, sampleRow,search.results)