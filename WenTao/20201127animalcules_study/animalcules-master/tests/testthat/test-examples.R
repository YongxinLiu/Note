data_dir = system.file("extdata/MAE.rds", package = "animalcules")
MAE = readRDS(data_dir)

# Summary Plot Top
test_that("filter_summary_bar_density() is working", {
    p <- filter_summary_bar_density(MAE,
                                    samples_discard = c("subject_2", "subject_4"),
                                    filter_type = "By Metadata",
                                    sample_condition = "AGE")
    
    expect_equal(length(p), 7)
})

# Summary Plot Bottom
test_that("filter_summary_bar_density() is working", {
    p <- filter_summary_bar_density(MAE,
                                    samples_discard = c("subject_2", "subject_4"),
                                    filter_type = "By Metadata",
                                    sample_condition = "SEX")
    expect_equal(length(p), 7)
})

# Categorize
test_that("filter_categorize() is working", {
    microbe <- MAE[['MicrobeGenetics']]
    samples <- as.data.frame(colData(microbe))
    result <- filter_categorize(samples,
                                sample_condition = "AGE",
                                new_label="AGE_GROUP",
                                bin_breaks=c(0,55,75,100),
                                bin_labels=c('Young','Adult',"Elderly"))
    expect_equal(dim(result$sam_table), c(50,5))
    expect_equal(length(result$plot.unbinned), 7)
    expect_equal(length(result$plot.binned), 7)
})

## Relative Abundance Stacked Bar Plot
test_that("relabu_barplot() is working", {
    p <- relabu_barplot(MAE,
                        tax_level="family",
                        order_organisms=c('Retroviridae'),
                        sort_by="organisms",
                        sample_conditions=c('SEX', 'AGE'),
                        show_legend=TRUE)
    expect_equal(length(p), 7)
})

## Relative Abundance Heatmap
test_that("relabu_heatmap() is working", {
    p <- relabu_heatmap(MAE,
                       tax_level="genus",
                       sort_by="conditions",
                       sample_conditions=c("SEX", "AGE"))
    expect_equal(length(p), 7)
})

## Relative Abundance Boxplot
test_that("relabu_boxplot() is working", {
    p <- relabu_boxplot(MAE,
                        tax_level="genus",
                        organisms=c("Escherichia", "Actinomyces"),
                        condition="SEX",
                        datatype="logcpm")
    expect_equal(length(p), 8)
})

## Alpha Diversity Boxplot
test_that("alpha_div_boxplot() is working", {
    p <- alpha_div_boxplot(MAE = MAE,
                      tax_level = "genus",
                      condition = "DISEASE",
                      alpha_metric = "shannon")
    expect_equal(length(p), 7)
})

## Alpha Diversity Statistical Test
test_that("do_alpha_div_test() is working", {
    p <- do_alpha_div_test(MAE = MAE,
                      tax_level = "genus",
                      condition = "DISEASE",
                      alpha_metric = "shannon",
                      alpha_stat = "T-test")
    pval_wil <- round(p[1,1,drop=TRUE], 2)
    pval_t <- round(p[1,2,drop=TRUE], 2)
    
    expect_equal(pval_wil, 0.33)
    expect_equal(pval_t, 0.42)
})

## Beta Diversity Heatmap
test_that("diversity_beta_heatmap() is working", {
    p <- diversity_beta_heatmap(MAE = MAE, 
                           tax_level = 'genus', 
                           input_beta_method = "bray",
                           input_bdhm_select_conditions = 'DISEASE',
                           input_bdhm_sort_by = 'condition')
    expect_equal(length(p), 7)
})

## Beta Diversity Boxplot
test_that("diversity_beta_boxplot() is working", {
    p <- diversity_beta_boxplot(MAE = MAE, 
                           tax_level = 'genus', 
                           input_beta_method = "bray",
                           input_select_beta_condition = 'DISEASE')
    expect_equal(length(p), 7)
})

## Beta Diversity Test
test_that("diversity_beta_test() is working", {
    p <- diversity_beta_test(MAE = MAE, 
                        tax_level = 'genus',
                        input_beta_method = "bray",
                        input_select_beta_condition =  'DISEASE',
                        input_select_beta_stat_method = 'PERMANOVA',
                        input_num_permutation_permanova = 999)
    expect_equal(p$Df[1], 1)
    expect_equal(p$Df[2], 48)
    expect_equal(p$Df[3], 49)
})

# PCA
test_that("dimred_pca() is working", {
    result <- dimred_pca(MAE,
                         tax_level="genus",
                         color="AGE",
                         shape="DISEASE",
                         pcx=1,
                         pcy=2,
                         datatype="logcpm")
    expect_equal(dim(result$table), c(50,4))
})

## PCoA
test_that("dimred_pcoa() is working", {
    result <- dimred_pcoa(MAE,
                          tax_level="genus",
                          color="AGE",
                          shape="DISEASE",
                          axx=1,
                          axy=2,
                          method="bray")
    expect_equal(dim(result$table), c(50,4))
})

## tSNE
test_that("dimred_tsne() is working", {
    result <- dimred_tsne(MAE,
                     tax_level="phylum",
                     color="AGE",
                     shape="GROUP",
                     k="3D",
                     initial_dims=30,
                     perplexity=10,
                     datatype="logcpm")
    expect_equal(length(result$plot), 8)
})

# Differential Analysis
test_that("differential_abundance() is working", {
    p <- differential_abundance(MAE,
                                tax_level="phylum",
                                input_da_condition=c("DISEASE"),
                                min_num_filter = 2,
                                input_da_padj_cutoff = 0.8)
    expect_equal(dim(p), c(8,8))
})

# Biomarker
test_that("find_biomarker() is working", {
    p <- find_biomarker(MAE,
                        tax_level="genus",
                        input_select_target_biomarker=c("SEX"),
                        nfolds = 3,
                        nrepeats = 3,
                        seed = 99,
                        percent_top_biomarker = 0.2,
                        model_name = "logistic regression")
    expect_equal(length(p), 3)
})

# inverse_simpson
test_that("inverse_simpson() is working", {
    p <- inverse_simpson(seq_len(10))
    p <- round(p, 2)
    expect_equal(p, 7.86)
})

# counts_to_relabu
test_that("counts_to_relabu() is working", {
    p <- counts_to_relabu(matrix(seq_len(12), 4))
    expect_equal(nrow(p), 4)
})

# counts_to_logcpm
test_that("counts_to_logcpm() is working", {
    p <- counts_to_logcpm(matrix(seq_len(12),4))
    expect_equal(nrow(p), 4)
})

# mae_pick_samples
test_that("mae_pick_samples() is working", {
    p <- mae_pick_samples(MAE, isolate_samples=c("subject_9", "subject_14"))
    expect_equal(length(p), 2)
})

# mae_pick_organisms
test_that("mae_pick_organisms() is working", {
    p <- mae_pick_organisms(MAE, isolate_organisms=c("ti|54005", "ti|73001"))
    expect_equal(length(p), 1)
})

# df_char_to_factor
test_that("df_char_to_factor() is working", {
    p <- df_char_to_factor(matrix(seq_len(12)))
    expect_equal(nrow(p), 12)
})

# percent
test_that("percent() is working", {
    p <- percent(c(0.42, 0.15, 0.4, 0.563, 0.2))
    expect_equal(p[1], "42.00%")
})

# is_categorical
test_that("is_categorical() is working", {
    p <- is_categorical(2)
    expect_true(p)
})

# is_integer0
test_that("is_integer0() is working", {
    p <- is_integer0(2)
    expect_false(p)
})

# is_integer1
test_that("is_integer1() is working", {
    p <- is_integer1(2)
    expect_false(p)
})

# pct2str
test_that("pct2str() is working", {
    p <- pct2str(0.23)
    expect_equal(p, "23.00")
})

# shannon
test_that("shannon() is working", {
    p <- shannon(seq_len(10))
    p <- round(p, 2)
    expect_equal(p, 2.15)
})

# gini_simpson
test_that("gini_simpson() is working", {
    p <- gini_simpson(seq_len(10))
    p <- round(p, 2)
    expect_equal(p, 0.87)
})

# grep_tid
test_that("grep_tid() is working", {
    p <- grep_tid("ti|700015|org|Coriobacterium_glomerans_PW2")
    expect_equal(p, "700015")
})

# find_taxonomy
test_that("find_taxonomy() is working", {
    p <- find_taxonomy(1200)
    expect_equal(p$Taxon$Taxon$TaxId, "131567")
})

# diversities
test_that("diversities() is working", {
    p <- diversities(matrix(seq_len(12), nrow = 3),index="shannon")
    expect_equal(round(p[1], 2), 1.01)
})

# diversities_help
test_that("diversities_help() is working", {
    p <- diversities_help(matrix(seq_len(12), nrow = 3),index="shannon")
    expect_equal(round(p[1], 2), 1.01)
})

# alpha_div_test
test_that("alpha_div_test() is working", {
  df_test <- data.frame(richness = seq_len(10),
                        condition = c(rep(1,5), 
                                      rep(0,5)))
  p <- alpha_div_test(df_test, alpha_stat="Wilcoxon rank sum test")
  expect_equal(round(as.numeric(as.character(p$output[2])), 4), 0.0011)
})
