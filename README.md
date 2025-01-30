# supervised_vs_unsupervised_bc_classification

For this project I looked at two grouping approaches, namely random forest classification and k-means clustering, and assessed how they performed on a data set consisting of a physical characteristics of breast cancer cells and an ultimate diagnosis as either malignant 'M' or benign 'B'.

The data is fairly interpretable, but here is a key for the first set of variables.

id: Cell ID number
diagnosis: he diagnosis of breast tissues (M = malignant, B = benign)
radius_mean: mean of distances from cell center to perimeter
texture_mean: standard deviation of gray-scale values
perimeter_mean: mean size of the core tumor
area_mean: mean area of tumor
smoothness_mean: mean of local variation in radius lengths
compactness_mean: mean of perimeter^2 / area - 1.0
concavity_mean: mean of severity of concave portions of the contour
concave_points: mean for number of concave portions of the contour

The code can be divided into three sections:
1) EDA: data exploration, pairs plot, correlation plot, bar graphs etc.

2) PCA: Conducting principal component analysis, Elbow plot to choose number of principal components

3) Random Forest Classification: producing random forest model, testing on test set, producing data frames for plotting, plotting results

4) K-means Classification: conducting K-means classification, testing on test set, producing data frames for plotting, plotting results
