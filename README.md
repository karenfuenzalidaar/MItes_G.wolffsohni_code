# MItes_G.wolffsohni_code
Code used for presentation at CLA. This code include data cleaning, PCA analysis, MANOVA, ANOVA, Tukey's HDS, and Ranfom forest. The data set come from project from University of Concepcion, Chile. The dataset includes diferent morphometrics measures in mites Gigantolaelaps wolffsohni from Chile.

library(readr)
library(tidyverse)
library(patchwork)
library(ggsignif)
library(GGally)
library(randomForest)
library(caret)

mites <- readr::read_csv("https://drive.google.com/uc?export=download&id=1zZs2mBKOuersjemQOEtIaMjP9HlULo63")
mites$ZONA <- as.factor(mites$ZONA)
str(mites)
mites_num <- mites[, sapply(mites, is.numeric)]
mites_num <- subset(mites_num, select = -n)
mites_log <- log1p(mites_num)
mites_log_clean <- as.data.frame(lapply(mites_log, function(x) {
  x[is.na(x) | is.infinite(x)] <- median(x, na.rm = TRUE)
  return(x)
}))

outliers_mad <- function(df, threshold = 3.5) {
  numeric_cols <- sapply(df, is.numeric)
  df_numeric <- df[, numeric_cols]
  
  outlier_list <- lapply(df_numeric, function(col) {
    med <- median(col, na.rm = TRUE)
    mad_val <- mad(col, na.rm = TRUE)
    if(mad_val == 0) return(NULL)
    modified_z <- 0.6745 * (col - med) / mad_val
    which(abs(modified_z) > threshold)
  })
  outlier_list[sapply(outlier_list, length) > 0]
}
outliers_detected <- outliers_mad(mites_log_clean)
outliers_detected
mites_no_outliers <- mites_log_clean  # copia de tus datos transformados

for(col in names(outliers_detected)){
  mites_no_outliers[outliers_detected[[col]], col] <- NA
}

# Imputar mediana
mites_no_outliers <- as.data.frame(lapply(mites_no_outliers, function(x){
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}))
sum(is.na(mites_no_outliers))

#PCA
pca_mites <- prcomp(mites_no_outliers, center = TRUE, scale. = TRUE)
pc_scores <- as.data.frame(pca_mites$x)
pc_scores$ZONA <- mites$ZONA  # tu factor de interés
str(pc_scores$ZONA)
pc_scores$ZONA <- factor(pc_scores$ZONA)

#plot PCA
explained_var <- pca_mites$sdev^2
explained_var_percent <- explained_var / sum(explained_var) * 100
df <- data.frame(
  Dimension = 1:length(explained_var_percent),
  Variance = explained_var_percent
) %>% 
  dplyr::slice(1:10)  
df$Dimension <- factor(df$Dimension, levels = 1:10)
ggplot(df, aes(x = Dimension, y = Variance, group = 1)) +
  geom_bar(stat = "identity", fill = "steelblue") +       
  geom_line(aes(y = Variance), color = "black") +         
  geom_point(aes(y = Variance), color = "black") +       
  labs(title = "Scree plot",
       x = "Dimensions",
       y = "Percentage of explained variance") +
  theme_minimal()

var_cos2 <- get_pca_var(pca_mites)$cos2[, 1:2]  
var_cos2_sum <- rowSums(var_cos2)
df_cos2 <- data.frame(
  Variable = names(var_cos2_sum),
  Cos2 = var_cos2_sum
)
df_cos2$Variable <- factor(df_cos2$Variable, levels = df_cos2$Variable[order(df_cos2$Cos2, decreasing = TRUE)])
ggplot(df_cos2, aes(x = Variable, y = Cos2)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "", y = "Quality of representation (Cos²)", 
       title = "Cos² de las variables en Dim1 y Dim2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Manova and assumptions
manova_model <- manova(as.matrix(pc_scores[, 1:5]) ~ ZONA, data = pc_scores)
summary(manova_model, test = "Pillai")
resid_new <- residuals(manova_model)

# Shapiro-Wilk por componente
shapiro_results <- apply(resid_new, 2, function(x) shapiro.test(x)$p.value)
shapiro_results
library(biotools)
boxM(pc_scores[, 1:5], pc_scores$ZONA)
resid_long <- pivot_longer(as.data.frame(resid_new), everything(), names_to = "PC", values_to = "Residual")
ggplot(resid_long, aes(sample = Residual)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~PC, scales = "free") +
  theme_minimal()
resid_new <- residuals(manova_model)

# 3. Transformar a data frame
resid_df <- as.data.frame(resid_new)
resid_df$ZONA <- pc_scores$ZONA

# 4. Crear todos los pares de PCs
pc_pairs <- combn(1:5, 2, simplify = FALSE)

# 5. Graficar
plots <- lapply(pc_pairs, function(pair){
  x_pc <- paste0("PC", pair[1])
  y_pc <- paste0("PC", pair[2])
  
  ggplot(resid_df, aes_string(x = x_pc, y = y_pc, color = "ZONA")) +
    geom_point(size = 2, alpha = 0.7) +
    theme_minimal() +
    labs(title = paste(x_pc, "vs", y_pc, "- Residuals"),
         x = x_pc, y = y_pc)
})
wrap_plots(plots, ncol = 2)

set.seed(123)
adonis_result <- adonis2(pc_scores[,1:5] ~ ZONA, data = pc_scores, permutations = 999, method = "euclidean")
adonis_result

#Individual ANOVA
summary(aov(fit))

#ASSUMPTIONS graphically
residuals <- residuals(fit)
residuals_df <- as.data.frame(residuals)
residuals_long <- pivot_longer(residuals_df, everything(), names_to = "Component", values_to = "Residual")

# Homogeneity of variance
anova_results <- summary.aov(fit)
anova_results  

residuals_df <- as.data.frame(residuals(fit))
residuals_df$ZONA <- pc_scores$ZONA

ggpairs(residuals_df,
        columns = 1:5,        
        mapping = aes(color = ZONA, alpha = 0.7),
        upper = list(continuous = wrap("points", size = 2)),
        lower = list(continuous = wrap("points", size = 2)),
        diag = list(continuous = wrap("densityDiag"))) +
  theme_minimal() +
  labs(title = "Scatter plots de residuos por pares de PCs")


# Multicollinearity
cor(pc_scores[, 1:5])

#Assumptions met. Still, using a robust test in MANOVA (Pillai) 
#for not completely clear homogeneity of variance of PC5

# Tukey HSD and Relationships between variables in boxplots

comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))

# Función para crear boxplot con Tukey HSD 
plot_pc_tukey <- function(pc_name, pc_scores) {
  aov_pc <- aov(as.formula(paste(pc_name, "~ ZONA")), data = pc_scores)
  tukey <- TukeyHSD(aov_pc)
  tukey_pvals <- tukey$ZONA[, "p adj"]
  labels <- ifelse(tukey_pvals < 0.001, "p < 0.001",
                   paste0("p = ", formatC(tukey_pvals, format = "f", digits = 3)))
  max_val <- max(pc_scores[[pc_name]], na.rm = TRUE)
  y_positions <- seq(max_val * 1.05, max_val * 1.25, length.out = length(labels))
  p <- ggplot(pc_scores, aes_string(x = "ZONA", y = pc_name, fill = "ZONA")) +
    geom_boxplot() +
    geom_signif(comparisons = comparisons,
                annotations = labels,
                y_position = y_positions,
                tip_length = 0.03,
                textsize = 5,   
                vjust = 0.3)+
    scale_x_discrete(labels = c("1" = "Mediterráneo",
                                "2" = "Bosque Templado",
                                "3" = "Bosque Magallánico")) +
    scale_fill_manual(values = c("1" = "#FDD835",  
                                 "2" = "#7CB342",  
                                 "3" = "#3366CC")) + 
    theme_bw() +
    theme(legend.position = "none", 
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5), 
          axis.text.y = element_text(size = 15),                           
          axis.title.y = element_text(size = 17)                          
          ) +
    labs(x = NULL, y = pc_name)
  return(p)
}

# Aplicar función a las 5 PCs
plot_pc_tukey("PC1", pc_scores)
plot_pc_tukey("PC2", pc_scores)
plot_pc_tukey("PC3", pc_scores)
plot_pc_tukey("PC4", pc_scores)
plot_pc_tukey("PC5", pc_scores)


####################

# Modelo Random Forest 

set.seed(123)
train_indices <- sample(1:nrow(mites), 0.7 * nrow(mites))
train_data <- mites[train_indices, ]
test_data <- mites[-train_indices, ]

model_rf <- randomForest(ZONA ~ ., data = train_data, ntree = 500)
print(model_rf)

predictions <- predict(model_rf, newdata = test_data)
conf_matrix_rf <- table(Predicted = predictions, Actual = test_data$ZONA)
print(conf_matrix_rf)

accuracy_rf <- sum(diag(conf_matrix_rf)) / sum(conf_matrix_rf)
accuracy_rf

importance(model_rf)
varImpPlot(model_rf)


# Random Forest con hyperparameter tuning 
set.seed(123)
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

tune_grid <- expand.grid(mtry = c(2, 5, 10)) 
model_tuned <- train(
  ZONA ~ .,
  data = train_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  ntree = 500
)

print(model_tuned)

final_predictions <- predict(model_tuned, newdata = test_data)
conf_matrix_tuned <- confusionMatrix(final_predictions, test_data$ZONA)
print(conf_matrix_tuned)

importance_df <- varImp(model_tuned)$importance
importance_df$Variable <- rownames(importance_df)

ggplot(importance_df, aes(x = reorder(Variable, Overall), y = Overall)) +
  geom_col(fill = "skyblue4") +
  coord_flip() +
  geom_text(aes(label = round(Overall, 1)), hjust = 1, size = 3, color="white") +
  scale_y_continuous(expand = c(0, 0.5)) +
  labs(
    title = "Importancia de Variables (Random Forest Tunado)",
    x = "Variable",
    y = "Importancia"
  ) +
  theme_bw() +
  theme(axis.ticks = element_blank())

accuracy_tuned <- conf_matrix_tuned$overall['Accuracy']
accuracy_tuned
sensitivity_tuned <- conf_matrix_tuned$byClass[, "Sensitivity"]
sensitivity_tuned 
specificity_tuned <- conf_matrix_tuned$byClass[, "Specificity"]
specificity_tuned
mean_sensitivity <- mean(sensitivity_tuned)
mean_sensitivity
mean_specificity <- mean(specificity_tuned)
mean_specificity
precision <- conf_matrix_tuned$byClass[, "Precision"]
f1_score <- 2 * (sensitivity_tuned * precision) / (sensitivity_tuned + precision)
mean_f1 <- mean(f1_score)
mean_f1
