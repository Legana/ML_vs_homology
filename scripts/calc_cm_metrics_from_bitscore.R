# function to calculate confusion matrix and associated metrics

calc_cm_metrics_from_bitscore <- function(p_threshold, df) {
  
  
  TP <- df %>% filter((Label=="Pos")) %>% filter(bitscore > p_threshold) %>% n_distinct()
  FP <- df %>% filter((Label=="Neg")) %>% filter(bitscore > p_threshold) %>% n_distinct()
  TN <- df %>% filter((Label=="Neg")) %>% filter(bitscore < p_threshold) %>% n_distinct()
  FN <- df %>% filter((Label=="Pos")) %>% filter(bitscore < p_threshold) %>% n_distinct()
  
  Specificity <- TN / (TN + FP) #aka TNR
  Recall <- TP / (TP + FN) # aka sensitivity, TPR
  Precision <- TP / (TP + FP)  # positive predictive value
  FPR <- FP / (TN + FP)

  cm <- c(TP, FP, TN, FN, Specificity, Recall, Precision, FPR, p_threshold)
  names(cm) <-c("TP", "FP", "TN", "FN", "Specificity", "Recall", "Precision", "FPR", "p_threshold") 
  cm
}
