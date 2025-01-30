

# CANCER CELL DIAGNOSTIC MODEL


# Setting Working Directory

getwd()

wd <- "/Users/alicethesen/Desktop/R/Practice_Data_Sets"

setwd(wd)

getwd()

#install.packages("factoextra")

# Installing packages

library(randomForest)
library(corrplot)
library(GGally)
library(gridExtra)
library(FactoMineR)
library(factoextra)
library(caret)

# Importing Data

bc_data <- read.csv("bcdata.csv")
bc_data$diagnosis <- as.factor(bc_data$diagnosis)
bc_data <-bc_data[,1:32] # removing x col

# Scaling Data: (so that modelling is not biased)
bc_data[,3:32] <- scale(bc_data[,3:32])

# ------------------------------------------------------------------------------
#                               EDA
# ------------------------------------------------------------------------------

head(bc_data)

# (1) Outlier Detection

radius_box <- ggplot(bc_data, aes(x=radius_mean)) + geom_boxplot(fill = "mediumpurple4") + 
  coord_flip() + theme(legend.position = "none")
texture_box <- ggplot(bc_data, aes(x=texture_mean)) + geom_boxplot(fill = "mediumpurple3") + 
  coord_flip() + theme(legend.position = "none")
perimeter_box <- ggplot(bc_data, aes(x=perimeter_mean)) + geom_boxplot(fill = "mediumpurple2") + 
  coord_flip() + theme(legend.position = "none")
area_box <- ggplot(bc_data, aes(x=area_mean)) + geom_boxplot(fill = "mediumpurple1") + 
  coord_flip() + theme(legend.position = "none")
smoothness_box <- ggplot(bc_data, aes(x=smoothness_mean)) + geom_boxplot( fill = "mediumpurple") + 
  coord_flip() + theme(legend.position = "none")
compactness_box <- ggplot(bc_data, aes(x=compactness_mean)) + geom_boxplot(fill = "mediumorchid1") + 
  coord_flip() + theme(legend.position = "none")
concavity_box <- ggplot(bc_data, aes(x=concavity_mean)) + geom_boxplot(fill = "mediumorchid2") + 
  coord_flip() + theme(legend.position = "none")
concave.points_box <- ggplot(bc_data, aes(x=concave.points_mean)) + geom_boxplot(fill = "mediumorchid3") + 
  coord_flip() + theme(legend.position = "none")
symmetry_box <- ggplot(bc_data, aes(x=symmetry_mean)) + geom_boxplot(fill = "mediumorchid") + 
  coord_flip() + theme(legend.position = "none")
fractal_dimension_box <- ggplot(bc_data, aes(x=fractal_dimension_mean)) + geom_boxplot(fill = "mediumorchid4") + 
  coord_flip() + theme(legend.position = "none")

grid.arrange(radius_box, texture_box, perimeter_box, area_box, smoothness_box, nrow = 1)

grid.arrange(compactness_box, concavity_box, concave.points_box, symmetry_box, fractal_dimension_box, nrow = 1)




# (2) Pairs Plot
pairs_plot <- ggpairs(bc_data[,2:12], aes(color = diagnosis ,alpha = 0.4))
pairs_plot


# Based on the above plot: radius, perimeter, area and concavity appear to show
# the greatest degree of separation between Malignant and Benign cells.


# (3) Correlation Plot
bc_cor <- cor(bc_data[,3:32])

par(mfrow=c(1,1))

corrplot(bc_cor, order = 'hclust', addrect = 2) 
# clustering argument groups most correlated variables








# --------------------------------------------------------------------------------------------
# Setting seed and taking samples:
# --------------------------------------------------------------------------------------------


set.seed(2323)

# Partition data into training and test data
samps <- sample(1:nrow(bc_data), size = round(0.8*length(bc_data[,1])), replace = FALSE)

train <- bc_data[samps,]
test <-bc_data[-samps,]




# It is clear from the correlation plots and the pairs plot that there are a 
# large amount of highly-correlated variables. 
# There is also a large number of dimensions. This data set could benefit from 
# dimension reduction techniques.


# --------------------------------------------------------------------------------------------
# PCA
# --------------------------------------------------------------------------------------------

train_PCA <- PCA(train[,3:32], scale.unit = TRUE, ncp = 17, graph = FALSE)
test_PCA <- PCA(test[,3:32], scale.unit = TRUE, ncp = 17, graph = FALSE)

# Elbow plot - package
fviz_screeplot(train_PCA, ncp=17)
# Seems to be clear elbow at 3 PC's

# Elbow plot (Manual) - looks at percentage of variance explained by PC's
elbow_plot <- barplot(train_PCA$eig[, 2], names.arg=1:nrow(train_PCA$eig), 
        main = "Variance Explained by PCs",
        xlab = "Principal Components",
        ylab = "Percentage of Variance",
        col ="powderblue")
# Add connected line segments to the plot
lines(x = elbow_plot, y = train_PCA$eig[, 2], 
      type="b", pch=19, col = "brown3")
# 95% of variation explained
abline(h= (train_PCA$eig[, 2][1] - (train_PCA$eig[, 2][1])*0.95), col = "orange", lwd = 1.5)
grid()




# Elbow at 3 PC's, but 7PC's explain 95% of variation (orange line)



# Structuring PCA in preparation for Random Forest:
#----------------------------------------------------------------------------------

# Isolate first 7 principal components (use the first few PCs that capture most variance)
train_PCs <- as.data.frame(train_PCA$ind$coord[, 1:7])  
# Add diagnosis column end to the data
train_PCs$diagnosis <- train$diagnosis  

# Repeat PCA on test data
test_PCs <- as.data.frame(test_PCA$ind$coord[, 1:7]) 
# Add diagnosis column end to the data
test_PCs$diagnosis <- test$diagnosis  



# --------------------------------------------------------------------------------------------
# Random Forest
# --------------------------------------------------------------------------------------------

# Use PC training data set to train a Random Forest model
set.seed(2323)
rf_model <- randomForest(diagnosis ~ ., data = train_PCs, importance = TRUE)
rf_pred <- predict(rf_model, test_PCs)

# Confusion Matrix
confusionMatrix(rf_pred, test_PCs$diagnosis)  # if using the caret package for metrics
# +- 92% accuracy
# error rate: +- 8%
# p-value: 6.043e-07

# Scatter Plot to Visualise Efficacy:

# Construct DF of Classified Diagnosis, Actual Diagnosis, PCs 1 & 2 and a Correctness Variable

rf_data_frame <- data.frame(True = test_PCs$diagnosis,Predicted = rf_pred,PC1 = test_PCs$Dim.1,
                        PC2 = test_PCs$Dim.2)
rf_data_frame$Correct <- rf_data_frame$True == rf_data_frame$Predicted


# In the following section, I play with how I think the data is best represented. 
# The plot I chose is the last one, named FINAL PLOT

# Plot of RF (1)
ggplot(rf_data_frame, aes(x = PC1, y = PC2, color = True, shape = Predicted)) +
  # Colour gives TRUE grouping, shape gives PREDICTED grouping
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  scale_shape_manual(values = c(16, 17)) +  # Different shapes for predictions
  # Incorrectly classified points are circled
  geom_point(data = subset(rf_data_frame, Correct == FALSE),
             aes(x = PC1, y = PC2),
             size = 5, shape = 21, fill = NA, color = "red", stroke = 1) +
  # Ellipse gives TRUE grouping
  stat_ellipse(data = rf_data_frame, 
               aes(x = PC1, y = PC2, group = True, color = True), 
               type = "norm", level = 0.95, size = 0.7, linetype = "dashed")+
  labs(
    title = "Random Forest Predictions with Misclassifications Highlighted",
    x = "Principal Component 1", y = "Principal Component 2",
    color = "True Diagnosis", shape = "Predicted Diagnosis"
  ) +
  theme_minimal()


# Plot of RF (2)
ggplot(rf_data_frame, aes(x = PC1, y = PC2, color = Predicted, shape = True)) +
  # Colour gives PREDICTED diagnosis, shape gives TRUE diagnosis
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  scale_shape_manual(values = c(16, 17)) +  # Different shapes for correct/incorrect
  # incorrectly classified points are circled
  geom_point(data = subset(rf_data_frame, Correct == FALSE),
  aes(x = PC1, y = PC2),
  size = 5, shape = 21, fill = NA, color = "red", stroke = 1.5) +
  # Ellipse gives PREDICTED grouping
  stat_ellipse(aes(x = PC1, y = PC2, group = Predicted, color = True), type = "norm", level = 0.95, 
                   linetype = "dashed", size = 0.7) +
  labs(
    title = "Random Forest Predictions with Misclassifications Highlighted",
    x = "Principal Component 1", y = "Principal Component 2",
    color = "Predicted Diagnosis", shape = "True Diagnosis"
      ) +
  theme_minimal()


# Sames as RF (2) but with different aesthetics
ggplot(rf_data_frame, aes(x = PC1, y = PC2, color = True, shape = Predicted)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Predicted), type = "norm", level = 0.95, size = 1, linetype = "dashed") +
  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
  labs(
    title = "Random Forest Predictions with Ellipses Around True Diagnosis",
    x = "Principal Component 1", y = "Principal Component 2",
    color = "True Diagnosis", shape = "Predicted Diagnosis"
  ) +
  theme_minimal()


# FINAL PLOT
# Same as above but has removed dashed from True Diagnosis legend key
RF_test_ggplot <-ggplot(rf_data_frame, aes(x = PC1, y = PC2, color = True, shape = Predicted)) +
                  # Colour gives TRUE diagnosis, Shape gives PREDICTED diagnosis
                  geom_point(size = 3, alpha = 0.7) + # Plotting Data Points
                  geom_point(data = subset(rf_data_frame, Correct == FALSE), # Circling incorrectly classified points
                             aes(x = PC1, y = PC2),
                             size = 5, shape = 21, fill = NA, color = "red", stroke = 1) +
                  # Ellipse gives PREDICTED diagnosis
                  stat_ellipse(aes(x = PC1, y = PC2, group = Predicted), type = "norm", level = 0.95, 
                               linetype = "dashed", inherit.aes = FALSE, size = 0.75, color = "darkgrey") +
                  scale_color_manual(values = c("#2E9FDF", "#FC4E07")) +
                  labs(
                    title = "Random Forest Model on Test Set",
                    x = "Principal Component 1", y = "Principal Component 2",
                    color = "True Diagnosis", shape = "Predicted Diagnosis"
                  ) +
                  theme_minimal()

RF_test_ggplot



# --------------------------------------------------------------------------------------------
# K-Means Clustering
# --------------------------------------------------------------------------------------------


# Ensure data is scaled for k-means. (done in data import step)

set.seed(2323)
# k-means clustering with k = 2
# must remove diagnosis variable - unsupervised techniques
kmeans_train <- kmeans(train_PCs[1:7], centers = 2, nstart = 25)
# must use PCA data for suitable comparison to RF model


# Create contingency table to compare clusters with actual diagnosis
table(Predicted = kmeans_train$cluster, Actual = train_PCs$diagnosis)

# Training results data frame
train_results <- data.frame(Predicted = kmeans_train$cluster, Actual = train_PCs$diagnosis)

# Adding a character vector for prediction
train_results$predicted_letter <- vector(mode = 'character', length = length(train_results$Predicted))

# Iterate over prediction variable to assign
for (i in 1:length(train_results$Predicted)){
  
  if(train_results$Predicted[i] == 1){
    #print("M")
    train_results$predicted_letter[i] = "M" 
  } else if (train_results$Predicted[i] == 2) {
    #print("B")
    train_results$predicted_letter[i] = "B" 
  }
}

# Adding PC's for plotting purposes
train_results <- data.frame(train_results, PC1 = train_PCs$Dim.1, PC2 = train_PCs$Dim.2)

# Adding 'Correct' Col for highlighting misclassified values
train_results$Correct <- train_results$Actual == train_results$predicted_letter
# data frame is finished and ready for plotting

# Visualise clusters with the actual diagnosis labels for comparison
#library(factoextra)
fviz_cluster(kmeans_train, data = train_PCs[1:7], geom = "point", ellipse.type = "convex", 
             habillage = train$diagnosis, palette = c("#2E9FDF", "#FC4E07"),
             ggtheme = theme_minimal(), main = "K-Means Clustering with True Diagnosis Labels")


# Visualising clusters with centroids & misclassified points circled

kmeans_training_ggplot<- ggplot(train_results, aes(x = PC1, y = PC2, color = Actual, shape = predicted_letter)) +
                        # Colour gives TRUE diagnosis, Shape gives PREDICTED diagnosis
                        geom_point(size = 3, alpha = 0.7) + # Plotting Data Points
                         geom_point(data = subset(train_results, Correct == FALSE), # Circling incorrectly classified points
                                   aes(x = PC1, y = PC2),
                                  size = 5, shape = 21, fill = NA, color = "red", stroke = 1) +
                        # Ellipse gives PREDICTED diagnosis
                        stat_ellipse(aes(x = PC1, y = PC2, group = predicted_letter), type = "norm", level = 0.95, 
                                     linetype = "dashed", inherit.aes = FALSE, size = 0.75, color = "darkgrey") +
                        scale_color_manual(values = c("palegreen3", "violetred3")) +
                        # center of cluster 1
                        geom_point(aes(x = kmeans_train$centers[,1][1], y = kmeans_train$centers[,2][1]),
                                   color = "black", size = 4, shape = 8) +
                        # center of cluster 2
                        geom_point(aes(x = kmeans_train$centers[,1][2], y = kmeans_train$centers[,2][2]),
                                   color = "black", size = 4, shape = 8) +
                        labs(
                          title = "K-Means Clustering on Training Data",
                          x = "Principal Component 1", y = "Principal Component 2",
                          color = "True Diagnosis", shape = "Predicted Diagnosis"
                        ) +
                        theme_minimal()

kmeans_training_ggplot


# Predict clusters for the test set using the training cluster centers
# --------------------------------------------------------------------------------------

# Calculate distances from each test point to each cluster center
test_distances <- as.matrix(dist(rbind(kmeans_train$centers, test_PCs[1:7])))[1:2, -c(1:2)]
test_clusters <- apply(test_distances, 2, which.min)

# Convert to a dataframe for easy comparison
test_results <- data.frame(Predicted = factor(test_clusters, levels = c(1, 2)), Actual = test_PCs$diagnosis)

# Make a character variable to convey prediction clearly
test_results$predicted_letter <- vector(mode = 'character', length = length(test_results$Predicted))

# Iterate over prediction variable to assign
for (i in 1:length(test_results$Predicted)){
  
  if(test_results$Predicted[i] == 1){
    #print("M")
    test_results$predicted_letter[i] = "M" 
  } else if (test_results$Predicted[i] == 2) {
    #print("B")
    test_results$predicted_letter[i] = "B" 
  }
}

# Adding 'Correct' Col
test_results <- data.frame(
  test_results,
  PC1 = test_PCs[,1],
  PC2 = test_PCs[,2],
  Correct = test_results$Actual == test_results$predicted_letter
)


# Step 3: Evaluate clustering performance on the test set
table(Predicted = test_results$predicted_letter, Actual = test_results$Actual)
# 11 points were misclassified 
# 91% accuracy on test set

# Step 4: Visualise clusters on test set
#library(factoextra)
fviz_cluster(list(data = test[,3:32], cluster = test_clusters), geom = "point",
             ellipse.type = "convex", palette = c("#2E9FDF", "#FC4E07"),
             ggtheme = theme_minimal(), main = "Test Set Clustering based on Train Centers",
             habillage = test$diagnosis)


# Visualise with color for actual diagnosis and shape for predicted clusters

kmeans_test_ggplot<- ggplot(test_results, aes(x = PC1, y = PC2, color = Actual, shape = predicted_letter)) +
                          # Colour gives TRUE diagnosis, Shape gives PREDICTED diagnosis
                          geom_point(size = 3, alpha = 0.7) + # Plotting Data Points
                          geom_point(data = subset(test_results, Correct == FALSE), # Circling incorrectly classified points
                                     aes(x = PC1, y = PC2),
                                     size = 5, shape = 21, fill = NA, color = "red", stroke = 1) +
                          # Ellipse gives PREDICTED diagnosis
                          stat_ellipse(aes(x = PC1, y = PC2, group = predicted_letter), type = "norm", level = 0.95, 
                                       linetype = "dashed", inherit.aes = FALSE, size = 0.75, color = "darkgrey") +
                          scale_color_manual(values = c("palegreen3", "violetred3")) +
                          # center of cluster 1
                          geom_point(aes(x = kmeans_train$centers[,1][1], y = kmeans_train$centers[,2][1]),
                                     color = "black", size = 4, shape = 8) +
                          # center of cluster 2
                          geom_point(aes(x = kmeans_train$centers[,1][2], y = kmeans_train$centers[,2][2]),
                                     color = "black", size = 4, shape = 8) +
                          labs(
                            title = "K-Means Clustering on Test Set",
                            x = "Principal Component 1", y = "Principal Component 2",
                            color = "True Diagnosis", shape = "Predicted Diagnosis"
                          ) +
                          theme_minimal()

kmeans_test_ggplot
