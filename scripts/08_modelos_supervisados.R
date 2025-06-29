#datos para los modelos
expr_matrix_all <- cbind(oligo_mat, micro_mat, astro_mat)
expr_matrix_all <- log2(expr_matrix_all + 1)

grupos <- c(
  rep("Oligo", ncol(oligo_mat)),
  rep("Microglia", ncol(micro_mat)),
  rep("Astro", ncol(astro_mat))
)

condiciones <- c(
  ifelse(grepl("AD", colnames(oligo_mat)), "AD", "Control"),
  ifelse(grepl("AD", colnames(micro_mat)), "AD", "Control"),
  ifelse(grepl("AD", colnames(astro_mat)), "AD", "Control")
)

gene_var <- apply(expr_matrix_all, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:200]

expr_subset <- expr_matrix_all[top_genes, ]
expr_df <- as.data.frame(t(expr_subset))

expr_df$Condicion <- condiciones

colnames(expr_df) <- make.names(colnames(expr_df))

dim(expr_df)
table(expr_df$Condicion)

#'Condicion' es factor
expr_df$Condicion <- as.factor(expr_df$Condicion)

#validación cruzada
set.seed(123)

ctrl <- trainControl(
method = "repeatedcv", 
number = 5,
repeats = 20,
classProbs = TRUE,
savePredictions = "final",
summaryFunction = twoClassSummary
)

#rf
model_rf <- train(Condicion ~ ., data = expr_df, method = "rf", trControl = ctrl, metric = "ROC")

#(SVM) lineal
model_svm <- train(Condicion ~ ., data = expr_df, method = "svmLinear", trControl = ctrl, metric = "ROC")

# Regresión logística penalizada (Lasso/GLM)
model_glm <- train(Condicion ~ ., data = expr_df, method = "glmnet", trControl = ctrl, metric = "ROC")

pred_rf <- model_rf$pred
pred_svm <- model_svm$pred
pred_glm <- model_glm$pred

pred_rf$obs <- factor(pred_rf$obs, levels = c("Control", "AD"))
pred_svm$obs <- factor(pred_svm$obs, levels = c("Control", "AD"))
pred_glm$obs <- factor(pred_glm$obs, levels = c("Control", "AD"))

roc_rf <- roc(pred_rf$obs, pred_rf$AD)
roc_svm <- roc(pred_svm$obs, pred_svm$AD)
roc_glm <- roc(pred_glm$obs, pred_glm$AD)

print(auc(roc_rf))
print(auc(roc_svm))
print(auc(roc_glm))

plot(roc_rf, col = "#AB82FF", main = "Curvas ROC - Modelos Comparados", legacy.axes = TRUE)
lines(roc_svm, col = "#B3EE3A")
lines(roc_glm, col = "palevioletred2")
legend("bottomright",
       legend = c(
         paste("Random Forest (AUC =", round(auc(roc_rf), 2), ")"),
         paste("SVM (AUC =", round(auc(roc_svm), 2), ")"),
         paste("GLMnet (AUC =", round(auc(roc_glm), 2), ")")
       ),
       col = c("#AB82FF", "#B3EE3A", "palevioletred2"), lwd = 2)

varImp(model_rf)
varImp(model_svm)
varImp(model_glm)

confusionMatrix(pred_rf$pred, pred_rf$obs)
confusionMatrix(pred_svm$pred, pred_svm$obs)
confusionMatrix(pred_glm$pred, pred_glm$obs)

#genes clave
rf_importance <- varImp(model_rf)
plot(rf_importance, top = 20, main = "Importancia de genes - Random Forest")

svm_importance <- varImp(model_svm)
plot(svm_importance, top = 20, main = "Importancia de genes - SVM")

coef(model_glm$finalModel, model_glm$bestTune$lambda)

imp_rf <- varImp(model_rf)$importance
imp_svm <- varImp(model_svm)$importance

genes_rf <- rownames(imp_rf)
genes_svm <- rownames(imp_svm)

common_rf_genes <- intersect(genes_rf, common_all_three)
common_svm_genes <- intersect(genes_svm, common_all_three)

common_rf_genes
common_svm_genes
