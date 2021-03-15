# animalcules.rds is an R list object that contains three dataframe: 
# count_table, gene_expression_table and metadata_table. These three
# tables come from a simulated dataset we built that contains 50 samples.
# All values are simulated using rnorm().


# TB_example_dataset.rds is an R MultiAssayExperiment object that contains
# two dataframe: count_table and metadata_table. These two tables come 
# from a published dataset from the paper:

# Botero LE, Delgado-Serrano L, Cepeda ML, Bustos JR, Anzola JM, 
# Del Portillo P, Robledo J, Zambrano MM. Respiratory tract clinical 
# sample selection for microbiota analysis in patients with pulmonary 
# tuberculosis. Microbiome. 2014 Aug 25;2:29. doi: 10.1186/2049-2618-2-29. 
# eCollection 2014. PubMed PMID: 25225609; PubMed Central PMCID: PMC4164332.

# Then we used a tool called PathoScope to process and generate the dataset.


# asthma_example_dataset.rds is an R MultiAssayExperiment object that 
# contains two dataframe: count_table and metadata_table. These two tables 
# come from a published dataset from the paper:

# Castro-Nallar E, Shen Y, Freishtat RJ, Pérez-Losada M, Manimaran S, Liu G, 
# Johnson WE, Crandall KA. Integrating microbial and host transcriptomics to 
# characterize asthma-associated microbial communities. BMC Med Genomics. 
# 2015 Aug 16;8:50. doi: 10.1186/s12920-015-0121-1. PubMed PMID: 26277095; 
# PubMed Central PMCID: PMC4537781.

# Then we used a tool called PathoScope to process and generate the dataset.


data_raw <- base::system.file("extdata/animalcules.rds", package = "animalcules") %>%
            base::readRDS()

se_mgx <- magrittr::use_series(data_raw, count_table) %>%
          base::data.matrix() %>%
          S4Vectors::SimpleList() %>%
          magrittr::set_names("MGX")

se_ge <- magrittr::use_series(data_raw, gene_expression_table) %>%
         base::data.matrix() %>%
         S4Vectors::SimpleList() %>%
         magrittr::set_names("GeneExpression")

se_colData <- magrittr::use_series(data_raw, metadata_table) %>%
              S4Vectors::DataFrame()

se_rowData <- magrittr::use_series(data_raw, tax_table) %>%
              base::data.frame() %>%
              dplyr::mutate_all(as.character) %>%
              dplyr::select(superkingdom, phylum, class, order, family, genus) %>%
              S4Vectors::DataFrame()

microbe_se <- SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                                         colData = se_colData,
                                                         rowData = se_rowData)

host_se <- SummarizedExperiment::SummarizedExperiment(assays = se_ge,
                                                      colData = se_colData)

mae_experiments <- S4Vectors::SimpleList(MicrobeGenetics = microbe_se, 
                                         HostGenetics = host_se)

MAE <- MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments, 
                                                  colData = se_colData)

saveRDS(MAE, "extdata/MAE.rds")

microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
tax_table <- as.data.frame(rowData(microbe)) # organism x taxlev
sam_table <- as.data.frame(colData(microbe)) # sample x condition
counts_table <- as.data.frame(assays(microbe))[,rownames(sam_table)] # organism x sample

toy_data <- list("tax_table"=tax_table,
                 "sam_table"=sam_table,
                 "counts_table"=counts_table)     

saveRDS(toy_data, "extdata/toy_data.rds")