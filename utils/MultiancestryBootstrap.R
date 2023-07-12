library(optparse)
library(PRROC)
source("../watershed.R")

# Helper function to remove NAs from vector
remove_na <- function(x) {
  return(x[!is.na(x)])
}

### Watershed parameters

nBootstrap = 100
number_of_dimensions = 4
n2_pair_pvalue_fraction = 0.02
binary_pvalue_threshold = 0.05



for(ancestry in c("African", "EastAsian", "NativeAmerican", "European")){
  
  ## We take input file for Watershed training, and learned parameters from the final Watershed model, in each population
  input_file <- paste0("Multiancestry/TOPMed/EURTraining-", ancestry, "Eval.txt.gz")
  watershed_model <- readRDS(paste0("Multiancestry/Results/TOPMED_", ancestry, "_Watershed_approximate_N2pair_0.02_pvalthresh_0.05_evaluation_model.rds"))
  
  
  data_input <- load_watershed_data(input_file, number_of_dimensions, n2_pair_pvalue_fraction, binary_pvalue_threshold)
  # Parse data_input for evaluation-related relevent fields
  feat_all <- data_input$feat
  discrete_outliers_all <- data_input$outliers_discrete
  binary_outliers_all <- data_input$outliers_binary
  fraction_binary_outliers_all <- data_input$fraction_outliers_binary
  N2_pairs <- data_input$N2_pairs
  
  ## Extract test data (ie N2 pairs)
  #(has to be done seperately for RIVER vs Watershed)
  if (number_of_dimensions == 1) {
    # Extraction of test data for RIVER
    feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
    discrete_outliers_test1 <- as.matrix(c(discrete_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)], discrete_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)]))
    discrete_outliers_test2 <- as.matrix(c(discrete_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)], discrete_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)]))
    binary_outliers_test1 <- as.matrix(c(binary_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)], binary_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)]))
    binary_outliers_test2 <- as.matrix(c(fraction_binary_outliers_all[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)], fraction_binary_outliers_all[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)]))
    # Absolute pvalues from test prediction data set (to be used for RNA-only analysis)
    real_valued_outliers_test1 <- -log10(abs(as.matrix(c(data_input$outlier_pvalues[!is.na(N2_pairs)][seq(from=1,to=sum(!is.na(N2_pairs)),by=2)], data_input$outlier_pvalues[!is.na(N2_pairs)][seq(from=2,to=sum(!is.na(N2_pairs)),by=2)]))) + 1e-7)
    
  } else {
    # Extraction of test data for Watershed
    feat_test <- rbind(feat_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], feat_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
    discrete_outliers_test1 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
    discrete_outliers_test2 <- rbind(discrete_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], discrete_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
    binary_outliers_test1 <- rbind(binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])
    binary_outliers_test2 <- rbind(fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),], fraction_binary_outliers_all[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),])
    # Absolute pvalues from test prediction data set (to be used for RNA-only analysis)
    real_valued_outliers_test1 <- -log10(abs(rbind(data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=1,to=sum(!is.na(N2_pairs)),by=2),], data_input$outlier_pvalues[!is.na(N2_pairs),][seq(from=2,to=sum(!is.na(N2_pairs)),by=2),])) + 1e-7)
  }
  
  #######################################
  ## Standardize Genomic Annotations (features)
  #######################################
  mean_feat <- apply(feat_all, 2, mean)
  sd_feat <- apply(feat_all, 2, sd)
  feat_test <- scale(feat_test, center=mean_feat, scale=sd_feat)
  nchoose <- nrow(feat_test) / 2
  bootstrapResult <- matrix(, nrow = nBootstrap, ncol = number_of_dimensions)
  
  ### We then sample half of the test set by nBootstrap times to estimate distribution of AUPRC in each population
  
  for(i in 1:nBootstrap){
    row.indices <- sample(1:nrow(feat_test), nchoose)
    
    posterior_info_test <- update_marginal_posterior_probabilities(feat_test[row.indices,], discrete_outliers_test1[row.indices,], watershed_model)
    posterior_prob_test <- posterior_info_test$probability  # Marginal posteriors
    posterior_pairwise_prob_test <- posterior_info_test$probability_pairwise  # Pairwise posteriors
    
    auprc <- c()
    
    for (dimension in 1:number_of_dimensions) {
      
      skip_to_next <- FALSE
      
      # Pseudo gold standard
      test_outlier_status <- binary_outliers_test2[row.indices, dimension]
      
      ## If there are too few evaluation instances (NA values), this analysis will fail, so we skip to the next bootstrapped samples
      tryCatch(watershed_pr_obj <- pr.curve(scores.class0 = remove_na(posterior_prob_test[,dimension][test_outlier_status==1 & !is.na(real_valued_outliers_test1[row.indices,dimension])]), 
                                            scores.class1 = remove_na(posterior_prob_test[,dimension][test_outlier_status==0 & !is.na(real_valued_outliers_test1[row.indices,dimension])]), curve = T), 
               error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next){
        auprc <- c(auprc, NA) 
      } else{
        auprc <- c(auprc, watershed_pr_obj$auc.integral)
      }
    }
    
    bootstrapResult[i,] <- auprc
    
  }
  
  saveRDS(bootstrapResult, file = paste0("Multiancestry/Results/TOPMED_", ancestry, "_watershedapproxiimate_bootstrap_", nBootstrap, "_result.RDS"))
  
}


