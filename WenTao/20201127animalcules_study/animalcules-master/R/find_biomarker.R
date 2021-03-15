#' Identify biomarkers
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The taxon level used for organisms
#' @param input_select_target_biomarker Which condition is the target condition
#' @param nfolds number of splits in CV
#' @param nrepeats number of CVs with different random splits
#' @param seed for repeatable research
#' @param percent_top_biomarker Top importance percentage to pick biomarker
#' @param model_name one of 'logistic regression', 'random forest'
#'
#' @return A list
#'
#' @import MultiAssayExperiment
#' @import glmnet
#' @import DMwR
#' @import caret
#' @import plotROC
#' @import forcats
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 geom_col aes coord_flip theme_bw coord_equal annotate
#'
#' @examples
#' data_dir = system.file('extdata/MAE.rds', package = 'animalcules')
#' toy_data <- readRDS(data_dir)
#' p <- find_biomarker(toy_data,
#'                     tax_level='genus',
#'                     input_select_target_biomarker=c('DISEASE'),
#'                     nfolds = 3,
#'                     nrepeats = 3,
#'                     seed = 99,
#'                     percent_top_biomarker = 0.2,
#'                     model_name = 'logistic regression')
#' p
#'
#'
#' @export
find_biomarker <- function(MAE, 
                            tax_level, 
                            input_select_target_biomarker, 
                            nfolds = 3, 
                            nrepeats = 3, 
                            seed = 99, 
                            percent_top_biomarker = 0.2, 
                            model_name = c("logistic regression", 
                                            "random forest")) {
    ## SEED bioC not suggesst add set seed function in R code set.seed(seed)
    ## tables from MAE
    microbe <- MAE[["MicrobeGenetics"]]  #double bracket subsetting is easier
    tax_table <- as.data.frame(rowData(microbe))  # organism x taxlev
    sam_table <- as.data.frame(colData(microbe))  # sample x condition
    counts_table <- 
    as.data.frame(assays(microbe))[, rownames(sam_table)]  # organism x sample
    ## shiny UI input object
    # Sum counts by taxon level and return log10 cpm
    logcpm_table <- counts_table %>% 
            upsample_counts(tax_table, tax_level) %>% 
            counts_to_logcpm() %>% 
            base::t() %>% base::as.data.frame()
    # add target variable
    logcpm_table[, "y"] <- 
    sam_table %>% dplyr::pull(input_select_target_biomarker)
    # set up classification model prameters
    fitControl <- caret::trainControl(method = "repeatedcv", 
                                    number = nfolds, 
                                    repeats = nrepeats, 
                                    classProbs = TRUE, 
                                    summaryFunction = twoClassSummary, 
                                    sampling = "smote", 
                                    savePredictions = TRUE)
    # choose different model
    if (model_name == "logistic regression") {
        model_fit <- 
        caret::train(y ~ ., data = logcpm_table, method = "glmnet", 
            tuneLength = 5, trControl = fitControl, metric = "ROC")
    } else if (model_name == "svm") {
        model_fit <-
        caret::train(y ~ ., data = logcpm_table, method = "svmLinear", 
            tuneLength = 5, trControl = fitControl, metric = "ROC")
    } else if (model_name == "gbm") {
        model_fit <- 
        caret::train(y ~ ., data = logcpm_table, 
            method = "gbm", trControl = fitControl, 
            tuneLength = 5, metric = "ROC", verbose = FALSE)
    } else if (model_name == "random forest") {
        model_fit <- 
        caret::train(y ~ ., data = logcpm_table, method = "ranger", 
        trControl = fitControl, tuneLength = 5, 
        metric = "ROC", importance = "impurity")
    }
    # process the importance score
    if (model_name == "svm") {
        svm_importance <- caret::varImp(model_fit)$importance
        svm_importance[, 2] <- NULL
        colnames(svm_importance) <- "importance"
        biomarker <- svm_importance %>% rownames_to_column() %>% 
            dplyr::rename(biomarker = rowname) %>% 
            dplyr::arrange(importance) %>% 
            dplyr::filter(importance > quantile(importance, 
            1 - percent_top_biomarker)) %>% 
            dplyr::select(biomarker) %>% .$biomarker
        importance_plot <- svm_importance %>% 
            rownames_to_column() %>% 
            dplyr::rename(biomarker = rowname) %>% 
            dplyr::arrange(importance) %>% 
            dplyr::filter(importance > quantile(importance, 
            1 - percent_top_biomarker)) %>% 
            dplyr::mutate(biomarker = forcats::fct_inorder(biomarker)) %>% 
            ggplot2::ggplot() + geom_col(aes(x = biomarker, y = importance)) +
            coord_flip() + theme_bw()
    } else {
        biomarker <- caret::varImp(model_fit)$importance %>% 
        base::as.data.frame() %>% 
            rownames_to_column() %>% 
            dplyr::rename(importance = Overall) %>% 
            dplyr::rename(biomarker = rowname) %>% 
            dplyr::arrange(importance) %>% 
            dplyr::filter(importance > 
                quantile(importance, 1 - percent_top_biomarker)) %>% 
            dplyr::select(biomarker) %>% .$biomarker
        importance_plot <- caret::varImp(model_fit)$importance %>% 
            base::as.data.frame() %>% 
            rownames_to_column() %>% 
            dplyr::rename(importance = Overall) %>% 
            dplyr::rename(biomarker = rowname) %>% 
            dplyr::arrange(importance) %>% 
            dplyr::filter(importance > 
                    quantile(importance, 1 - percent_top_biomarker)) %>% 
            dplyr::mutate(biomarker = forcats::fct_inorder(biomarker)) %>% 
            ggplot2::ggplot() + 
            geom_col(aes(x = biomarker, y = importance)) + 
            coord_flip() + theme_bw()
    }
    # retrain the model using the biomarker
    logcpm_table <- logcpm_table %>% dplyr::select(biomarker, y)
    # choose different model
    if (model_name == "logistic regression") {
        model_fit <- caret::train(y ~ ., 
                data = logcpm_table, method = "glmnet", 
                tuneLength = 5, trControl = fitControl, metric = "ROC")
    } else if (model_name == "svm") {
        model_fit <- caret::train(y ~ ., 
                data = logcpm_table, method = "svmLinear", 
            tuneLength = 5, trControl = fitControl, metric = "ROC")
    } else if (model_name == "gbm") {
        model_fit <- caret::train(y ~ ., 
                data = logcpm_table, method = "gbm", trControl = fitControl, 
            tuneLength = 5, metric = "ROC", verbose = FALSE)
    } else if (model_name == "random forest") {
        model_fit <- caret::train(y ~ ., 
                data = logcpm_table, method = "ranger", 
            trControl = fitControl, tuneLength = 5, 
            metric = "ROC", importance = "impurity")
    }
    prob_pred <- as.numeric(model_fit$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    df_roc <- data.frame(m = model_fit$pred[, 
        which(colnames(model_fit$pred) == levels(model_fit$pred$obs)[2])], 
        d = prob_pred, stringsAsFactors = FALSE)
    g <- ggplot(df_roc, aes(m = m, d = d)) + 
        geom_roc(n.cuts = 0) + coord_equal() + 
        style_roc()
    roc_plot <- g + annotate("text", x = 0.75, y = 0.25, 
            label = paste("AUC =", round((calc_auc(g))$AUC, 
        4)))
    biomarker <- data.frame(biomarker_list = biomarker)
    # output a list
    list_output <- list(biomarker = biomarker, 
                        importance_plot = importance_plot, 
                        roc_plot = roc_plot)
    return(list_output)
}