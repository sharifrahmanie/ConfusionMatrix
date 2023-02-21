ConfusionMatrix <- function(y_true, y_pred){
  # Designed by @biomedical_informatics. Edris Sharif Rahmani Feb 21, 2023
  
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' required but not installed.")
  }
  
  if(!is.null(y_true) && !is.null(y_pred)){
    y_true <- factor(y_true)
    y_pred <- factor(y_pred)
    cm <- table(y_true, y_pred)
    missing_1 <- which(!colnames(cm) %in% rownames(cm))
    missing_2 <- which(!rownames(cm) %in% colnames(cm))
    if(length(missing_1) > 0) {
      zero_list <- list()
      for(i in seq_along(missing_1)){
        zero_list[[colnames(cm)[missing_1][i]]] <- rep(0, NCOL(cm))
      }
      zero_df <- do.call(rbind, zero_list)
      cm <- rbind(zero_df, cm)
      cm <- cm[colnames(cm),]
    }
    else if(length(missing_2) > 0){
      zero_list <- list()
      for(i in seq_along(missing_2)){
        zero_list[[rownames(cm)[missing_2][i]]] <- rep(0, NROW(cm))
      }
      zero_df <- do.call(cbind, zero_list)
      cm <- cbind(zero_df, cm)
      cm <- cm[,rownames(cm)]
    } else {
      cm <- cm[,rownames(cm)]
    }

    if(NROW(cm) == 2 & NCOL(cm) == 2){
      Accuracy <- (cm[1,1] + cm[2,2])/(cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
      Precision <- (cm[1,1])/(cm[1,1] + cm[1,2])
      Sensitivity <- (cm[1,1])/(cm[1,1] + cm[2,1])
      Specificity <- (cm[2,2])/ (cm[2,2] + cm[1,2])
      AUC <- pROC::roc(as.numeric(y_true) ~ as.numeric(y_pred), quiet = T)$auc[1]
      f1_score <- 2*(Precision*Sensitivity)/(Precision + Sensitivity)
      result <- round(data.frame(Accuracy = Accuracy, 
                                 Precision = Precision, 
                                 Sensitivity = Sensitivity , 
                                 F1_Score = f1_score, 
                                 Specificity = Specificity, 
                                 AUC = AUC),3)
      return(result)
    }
    else if (NROW(cm) > 2 | NCOL(cm) > 2){
      TP <- list()
      for(i in 1:NROW(cm)){
        TP[[i]] <- cm[i,i]
      }
      TP <- data.frame(do.call(rbind,TP))
      FN <- rowSums(cm) - TP[,1]
      FP <- colSums(cm) - TP[,1]
      TN <- list()
      for(i in 1:NROW(cm)){
        TN[[i]] <- sum(cm) - sum(cm[i,]) - sum(cm[,i]) + cm[i,i]
      }
      TN <- data.frame(do.call(rbind, TN))
      con <- cbind(TP, FN, FP, TN)
      colnames(con) <- c("TP", "FN", "FP", "TN")
      rownames(con) <- rownames(cm)
      a <- list()
      for (i in 1:NROW(con)) {
        a[[i]] <- (con$TP[i] + con$TN[i])/(con$TP[i] + con$TN[i] + con$FN[i] + con$FP[i])
      }
      p <- list()
      for (i in 1:NROW(con)) {
        p[[i]] <- con$TP[i]/(con$TP[i] + con$FP[i])
      }
      se <- list()
      for (i in 1:NROW(con)) {
        se[[i]] <- con$TP[i]/(con$TP[i] + con$FN[i])
      }
      sp <- list()
      for (i in 1:NROW(con)) {
        sp[[i]] <- con$TN[i]/(con$TN[i] + con$FP[i])
      }
      a <- do.call(rbind, a)
      p <- do.call(rbind, p)
      se <- do.call(rbind, se)
      f1 <- 2*(p*se)/(p + se)
      sp <- do.call(rbind, sp)
      au <- pROC::multiclass.roc(as.numeric(y_true) ~ as.numeric(y_pred), quiet = T)$auc[1]
      au <- rep(au, length.out= NROW(con))
      result <- round(cbind(a, p, se, f1, sp, au),3)
      colnames(result) <- c("Accuracy", 
                            "Precision", 
                            "Sensitivity", 
                            "F1_Score", 
                            "Specificity", 
                            "AUC_average")
      rownames(result) <- rownames(cm)
      return(result)
    }
  }
}
