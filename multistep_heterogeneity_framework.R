#loading libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, caret, twang, stats, MatchIt, tidyverse,labelled, tableone, rbounds, cobalt, marginaleffects, 
               boot, factoextra, cluster, grf, sandwich, lmtest, Hmisc, ggplot2, gridExtra, corrplot)
#load data
setwd("~/NSLM")
df_includingDescription <- read_excel("data_clean.xlsx", 
                                      sheet = "Sheet3")
#prepare data
df <- df_includingDescription
colnames(df)[2:13] <- c("Z", "Y", "S3", "C1", "C2", "C3", "XC", "X1", "X2", "X3", "X4", "X5")
df$schoolID <- as.factor(df$schoolID)
df$Z <- as.factor(df$Z)
df$S3 <- as.factor(df$S3)
df$C1 <- as.factor(df$C1)
df$C2 <- as.factor(df$C2)
df$C3 <- as.factor(df$C3)
df$XC <- as.factor(df$XC)

#one-hot encoding for the second part of the analysis
dummy <- dummyVars(~ S3 + C1 + C2 + C3 + XC, data=df)
dummies.df <- data.frame(predict(dummy, newdata=df))
df_one_hot <- cbind(df, dummies.df)
df_one_hot <- df_one_hot[,-c(4:8)]
df_one_hot$Z <- as.integer(df_one_hot$Z) #change treatment from type factor to integer
df_one_hot$Z <- ifelse(df_one_hot$Z == 1, 0, df_one_hot$Z) #change no-treatment to zero
df_one_hot$Z <- ifelse(df_one_hot$Z == 2, 1, df_one_hot$Z) #change treatment to 1

#exploratory analysis for non-random treatment assignment#
#perform logit PS score estimation for exploratory analysis including PS
m.out.logit <- matchit(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3.1 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 + C1.1 + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7
                       + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13 + C1.14 + C1.15 + C2.1 + C2.2 + C3.0 + C3.1 + XC.0 + XC.1 + XC.2 + XC.3 + XC.4,
                       data=df_one_hot,
                       distance='logit')
m.out.logit$distance

#add PS scores as a new column to the original df
dfPS <- df
dfPS["PS"] <- m.out.logit$distance

#make a boxplot, visualizing the dependence of PS on covariate
#create the boxplot
boxplot(dfPS$PS ~ dfPS$S3,
        outline = FALSE,
        xlim = c(0.5, 7.5),
        ylim = c(0.16, 0.48),
        xlab = "Student Expectation of Success",
        ylab = "Propensity Score")
#perform lowess smoothing
smoothed_line <- lowess(dfPS$PS ~ dfPS$S3)
#add smoothed line to the plot
lines(smoothed_line, col = "red")

##
##propensity score estimation##
##

#create df for gbm
df_gbm <- data.frame(df)
df_gbm$Z <- as.integer(df_gbm$Z) #change treatment from type factor to integer
df_gbm$Z <- ifelse(df_gbm$Z == 1, 0, df_gbm$Z) #change control to zero
df_gbm$Z <- ifelse(df_gbm$Z == 2, 1, df_gbm$Z) #change treatment to 1

#calculate propensity score using GBM, manual grid search
ps_gbm = ps(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3 + C1 + C2 + C3 + XC,  
            data = df_gbm,
            n.trees=20000,
            interaction.depth=2, 
            shrinkage=0.001,
            estimand = "ATT",
            stop.method=c("es.mean", "ks.max"),
            verbose=FALSE)
'optimal set of (n.trees, interaction.depth, shrinkage) = (20000, 2, 0.001), measured by minimizing the ASAM' 

#optimal number of iterations by minimizing the ASAM
ps_gbm$desc$es.mean.ATT$n.trees #17097 iterations
ps_gbm$desc$ks.max.ATT$n.trees #17916 iterations

#plot of the balance measures as a function of the number of iterations
plot(ps_gbm)
summary(ps_gbm$gbm.obj, plot=TRUE)

#calculate propensity scores using logistic model (benchmark)
ps.logit <- glm(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3 + C1 + C2 + C3 + XC,
                data= df_gbm,
                family=binomial)
ps.logit$fitted.values #propensity scores

##
##cluster on school-level variables##
##

#isolate the school-level variables
cluster.prep <- df_one_hot[,c(4:8,35:39)]
#choose best cluster
#Elbow Method
set.seed(111)
fviz_nbclust(cluster.prep, kmeans, method="wss", nstart = 18)

#Average Silhouette Method
set.seed(111)
fviz_nbclust(cluster.prep, kmeans, method = "silhouette", nstart = 18)
' --> 6 clusters seems to be the best outcome based on the Elbow- and Average Silhouette method'

#perform kmeans clustering using 6 clusters
set.seed(111)
k2 <- kmeans(cluster.prep, centers = 6, nstart = 25)
k2$cluster
print(k2)

#visualize clusters
fviz_cluster(k2, data = cluster.prep)
fviz_cluster(ellipse.type = "norm", k2, data = cluster.prep, geom = "point")
fviz_cluster(geom = "point", k2, data=cluster.prep)

#descriptives
cluster.prep %>%
  mutate(Cluster = k2$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")

#adding the cluster numbers to the rows
clusters <- df_one_hot %>%
  mutate(Cluster = k2$cluster)

#add PS to dataframe for es.mean
cluster_ps <- cbind(clusters, ps_gbm$ps$es.mean.ATT)

#do the same for ks.max
cluster_ps.ks <- cbind(clusters, ps_gbm$ps$ks.max.ATT) #add PS to dataframe

#do the same for benchmark model
cluster_ps.logit <- cbind(clusters, ps.logit$fitted.values)

#make a dataframe for each cluster
#es.mean
cluster_list <- lapply(1:6, function(i) cluster_ps %>% filter(Cluster == i))
cluster1 <- cluster_list[[1]]
cluster2 <- cluster_list[[2]]
cluster3 <- cluster_list[[3]]
cluster4 <- cluster_list[[4]]
cluster5 <- cluster_list[[5]]
cluster6 <- cluster_list[[6]]

#do the same for ks.max
cluster_list.ks <- lapply(1:6, function(i) cluster_ps.ks %>% filter(Cluster == i))
cluster1.ks <- cluster_list.ks[[1]]
cluster2.ks <- cluster_list.ks[[2]]
cluster3.ks <- cluster_list.ks[[3]]
cluster4.ks <- cluster_list.ks[[4]]
cluster5.ks <- cluster_list.ks[[5]]
cluster6.ks <- cluster_list.ks[[6]]

#do the same for benchmark model
cluster_list.logit <- lapply(1:6, function(i) cluster_ps.logit %>% filter(Cluster == i))
cluster1.logit <- cluster_list.logit[[1]]
cluster2.logit <- cluster_list.logit[[2]]
cluster3.logit <- cluster_list.logit[[3]]
cluster4.logit <- cluster_list.logit[[4]]
cluster5.logit <- cluster_list.logit[[5]]
cluster6.logit <- cluster_list.logit[[6]]

#do the same for the unadjusted sample
cluster_list.un <- lapply(1:6, function(i) clusters %>% filter(Cluster == i))
cluster1.un <- cluster_list.un[[1]]
cluster2.un <- cluster_list.un[[2]]
cluster3.un <- cluster_list.un[[3]]
cluster4.un <- cluster_list.un[[4]]
cluster5.un <- cluster_list.un[[5]]
cluster6.un <- cluster_list.un[[6]]

##
##propensity score matching##
##

#make list of clusters for each method
match_cluster_list <- list(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6)
match_cluster_list.ks <- list(cluster1.ks, cluster2.ks, cluster3.ks, cluster4.ks, cluster5.ks, cluster6.ks)
match_cluster_list.logit <- list(cluster1.logit, cluster2.logit, cluster3.logit, cluster4.logit, cluster5.logit, cluster6.logit)

#make empty lists for each model
match_gbm <- list()
match_gbm.ks <- list()
match_gbm.optimal <- list()
match_gbm.optimal.ks <- list()
match_gbm.logit <- list()

#use for loop to iterate over clusters
for (i in 1:6) {
  match_cluster <- match_cluster_list[[i]]
  match_cluster.ks <- match_cluster_list.ks[[i]]
  match_cluster.logit <- match_cluster_list.logit[[i]]
  
  #make matching models
  match_gbm[[i]] <- matchit(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3.1 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.1 + C1.2 + C1.3 + C1.4 + C1.5 + C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.1 + C2.2 + C3.0 + C3.1 + XC.0 + XC.1 + XC.2 + XC.3 + XC.4,
                            data = match_cluster, method = "nearest", caliper = 0.2,
                            distance = match_cluster$`ps_gbm$ps$es.mean.ATT`)
  
  match_gbm.ks[[i]] <- matchit(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3.1 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.1 + C1.2 + C1.3 + C1.4 + C1.5 + C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.1 + C2.2 + C3.0 + C3.1 + XC.0 + XC.1 + XC.2 + XC.3 + XC.4,
                               data = match_cluster.ks, method = "nearest", caliper = 0.2,
                               distance = match_cluster.ks$`ps_gbm$ps$ks.max.ATT`)
  
  match_gbm.optimal[[i]] <- matchit(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3.1 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                    + C1.1 + C1.2 + C1.3 + C1.4 + C1.5 + C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                    + C1.14 + C1.15 + C2.1 + C2.2 + C3.0 + C3.1 + XC.0 + XC.1 + XC.2 + XC.3 + XC.4,
                                    data = match_cluster, method = "optimal",
                                    distance = match_cluster$`ps_gbm$ps$es.mean.ATT`)
  
  match_gbm.optimal.ks[[i]] <- matchit(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3.1 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                       + C1.1 + C1.2 + C1.3 + C1.4 + C1.5 + C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                       + C1.14 + C1.15 + C2.1 + C2.2 + C3.0 + C3.1 + XC.0 + XC.1 + XC.2 + XC.3 + XC.4,
                                       data = match_cluster, method = "optimal",
                                       distance = match_cluster.ks$`ps_gbm$ps$ks.max.ATT`)
  
  match_gbm.logit[[i]] <- matchit(Z ~ schoolID + X1 + X2 + X3 + X4 + X5 + S3.1 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.1 + C1.2 + C1.3 + C1.4 + C1.5 + C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.1 + C2.2 + C3.0 + C3.1 + XC.0 + XC.1 + XC.2 + XC.3 + XC.4,
                                  data = match_cluster.logit, method = "nearest", caliper = 0.2,
                                  distance = match_cluster.logit$`ps.logit$fitted.values`)
}

#retrieve the matchit object for each cluster
#es.mean
match_gbm_cluster1 <- match_gbm[[1]]
match_gbm_cluster2 <- match_gbm[[2]]
match_gbm_cluster3 <- match_gbm[[3]]
match_gbm_cluster4 <- match_gbm[[4]]
match_gbm_cluster5 <- match_gbm[[5]]
match_gbm_cluster6 <- match_gbm[[6]]

#ks.max
match_gbm_cluster1.ks <- match_gbm.ks[[1]]
match_gbm_cluster2.ks <- match_gbm.ks[[2]]
match_gbm_cluster3.ks <- match_gbm.ks[[3]]
match_gbm_cluster4.ks <- match_gbm.ks[[4]]
match_gbm_cluster5.ks <- match_gbm.ks[[5]]
match_gbm_cluster6.ks <- match_gbm.ks[[6]]

#optimal with es.mean
match_gbm_cluster1.optimal <- match_gbm.optimal[[1]]
match_gbm_cluster2.optimal <- match_gbm.optimal[[2]]
match_gbm_cluster3.optimal <- match_gbm.optimal[[3]]
match_gbm_cluster4.optimal <- match_gbm.optimal[[4]]
match_gbm_cluster5.optimal <- match_gbm.optimal[[5]]
match_gbm_cluster6.optimal <- match_gbm.optimal[[6]]

#optimal with ks.max
match_gbm_cluster1.optimal.ks <- match_gbm.optimal.ks[[1]]
match_gbm_cluster2.optimal.ks <- match_gbm.optimal.ks[[2]]
match_gbm_cluster3.optimal.ks <- match_gbm.optimal.ks[[3]]
match_gbm_cluster4.optimal.ks <- match_gbm.optimal.ks[[4]]
match_gbm_cluster5.optimal.ks <- match_gbm.optimal.ks[[5]]
match_gbm_cluster6.optimal.ks <- match_gbm.optimal.ks[[6]]

#logit
match_gbm_cluster1.logit <- match_gbm.logit[[1]]
match_gbm_cluster2.logit <- match_gbm.logit[[2]]
match_gbm_cluster3.logit <- match_gbm.logit[[3]]
match_gbm_cluster4.logit <- match_gbm.logit[[4]]
match_gbm_cluster5.logit <- match_gbm.logit[[5]]
match_gbm_cluster6.logit <- match_gbm.logit[[6]]

#summary and plot of the matching procedure for cluster 1 (can be applied to each individual cluster)
summary(match_gbm_cluster1)
plot(match_gbm_cluster1, type="hist")
summary(match_gbm_cluster1.ks)
plot(match_gbm_cluster1.ks, type="hist")
summary(match_gbm_cluster1.optimal)
plot(match_gbm_cluster1.optimal, type="hist")
summary(match_gbm_cluster1.optimal.ks)
plot(match_gbm_cluster1.optimal.ks, type="hist")
summary(match_gbm_cluster1.logit)
plot(match_gbm_cluster1.logit, type="hist")

#check balance for each of the matching methods
#use bal.tab since CreateOneTable overestimates the SMD (specifically for discrete variables - Qin-Yu Zhao, 2021)
bal_table_cluster1 <- bal.tab(match_gbm_cluster1, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster2 <- bal.tab(match_gbm_cluster2, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster3 <- bal.tab(match_gbm_cluster3, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster4 <- bal.tab(match_gbm_cluster4, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster5 <- bal.tab(match_gbm_cluster5, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster6 <- bal.tab(match_gbm_cluster6, m.threshold = 0.1, un = TRUE, abs = TRUE)
#do the same for ks.max
bal_table_cluster1.ks <- bal.tab(match_gbm_cluster1.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster2.ks <- bal.tab(match_gbm_cluster2.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster3.ks <- bal.tab(match_gbm_cluster3.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster4.ks <- bal.tab(match_gbm_cluster4.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster5.ks <- bal.tab(match_gbm_cluster5.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster6.ks <- bal.tab(match_gbm_cluster6.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
#do the same for optimal matching
bal_table_cluster1.optimal <- bal.tab(match_gbm_cluster1.optimal, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster2.optimal <- bal.tab(match_gbm_cluster2.optimal, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster3.optimal <- bal.tab(match_gbm_cluster3.optimal, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster4.optimal <- bal.tab(match_gbm_cluster4.optimal, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster5.optimal <- bal.tab(match_gbm_cluster5.optimal, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster6.optimal <- bal.tab(match_gbm_cluster6.optimal, m.threshold = 0.1, un = TRUE, abs = TRUE)
#repeat for ks.max
bal_table_cluster1.optimal.ks <- bal.tab(match_gbm_cluster1.optimal.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster2.optimal.ks <- bal.tab(match_gbm_cluster2.optimal.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster3.optimal.ks <- bal.tab(match_gbm_cluster3.optimal.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster4.optimal.ks <- bal.tab(match_gbm_cluster4.optimal.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster5.optimal.ks <- bal.tab(match_gbm_cluster5.optimal.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster6.optimal.ks <- bal.tab(match_gbm_cluster6.optimal.ks, m.threshold = 0.1, un = TRUE, abs = TRUE)
#repeat for benchmark
bal_table_cluster1.logit <- bal.tab(match_gbm_cluster1.logit, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster2.logit <- bal.tab(match_gbm_cluster2.logit, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster3.logit <- bal.tab(match_gbm_cluster3.logit, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster4.logit <- bal.tab(match_gbm_cluster4.logit, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster5.logit <- bal.tab(match_gbm_cluster5.logit, m.threshold = 0.1, un = TRUE, abs = TRUE)
bal_table_cluster6.logit <- bal.tab(match_gbm_cluster6.logit, m.threshold = 0.1, un = TRUE, abs = TRUE)

#define a list of the match_gbm_clusters and their variations
cluster_objects <- list(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6)
clusters_list <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6")
variations <- c("", ".ks", ".optimal", ".optimal.ks")

#create an empty dataframe to store the results
df_bal_tables <- data.frame(Cluster = character(), Variation = character(), Mean_Diff_Adj = numeric())#, Matched_Sample_Size = integer())

#loop through clusters and variations
for (cluster in clusters_list) {
  for (variation in variations) {
    #construct the bal_table object name dynamically
    bal_table_name <- paste("bal_table_", cluster, variation, sep = "")
    #get the appropriate bal_table object
    bal_table_object <- get(bal_table_name)
    #calculate the mean standardized effect size for each balanced cluster per method
    mean_diff_adj <- mean(bal_table_object$Balance$Diff.Adj)
    #create a new row for the dataframe
    new_row <- data.frame(Cluster = cluster, Variation = variation, Mean_Diff_Adj = mean_diff_adj) #,Matched_Sample_Size = matched_sample_size)
    #append the new row to the df_bal_tables
    df_bal_tables <- rbind(df_bal_tables, new_row)
  }
}
#print the resulting dataframe for comparison between methods
print(df_bal_tables)

#bal_table for benchmark model
#define a list of the match_gbm_clusters and their variations
cluster_objects.logit <- list(cluster1.logit, cluster2.logit, cluster3.logit, cluster4.logit, cluster5.logit, cluster6.logit)
clusters_list.logit <- c("cluster1.logit", "cluster2.logit", "cluster3.logit", "cluster4.logit", "cluster5.logit", "cluster6.logit")
variations.logit <- c("")

#create an empty dataframe to store the results
df_bal_tables.logit <- data.frame(Cluster = character(), Variation = character(), Mean_Diff_Adj = numeric())#, Matched_Sample_Size = integer())

#loop through clusters and variations
for (cluster in clusters_list.logit) {
  for (variation in variations.logit) {
    #construct the bal_table object name dynamically
    bal_table_name <- paste("bal_table_", cluster, variation, sep = "")
    #get the appropriate bal_table object
    bal_table_object <- get(bal_table_name)
    #calculate the mean standardized effect size for each balanced cluster per method
    mean_diff_adj <- mean(bal_table_object$Balance$Diff.Adj)
    #create a new row for the dataframe
    new_row <- data.frame(Cluster = cluster, Variation = variation, Mean_Diff_Adj = mean_diff_adj) #,Matched_Sample_Size = matched_sample_size)
    #append the new row to the df_bal_tables
    df_bal_tables.logit <- rbind(df_bal_tables.logit, new_row)
  }
}
print(df_bal_tables.logit)

#balance visualizations
#make vector of cluster names
matched_cluster_names <- c("match_gbm_cluster1", "match_gbm_cluster2", "match_gbm_cluster3",
                           "match_gbm_cluster4", "match_gbm_cluster5", "match_gbm_cluster6")

#use a for loop to iterate through cluster names
for (matched_cluster_name in matched_cluster_names) {
  cluster <- get(matched_cluster_name)
  
  #use love.plot() for the current cluster
  print(love.plot(bal.tab(cluster, m.threshold = 0.1),
                  stat = "mean.diffs",
                  grid = TRUE,
                  stars = "raw",
                  abs = FALSE))
}

#choose colors
custom_colors <- c("black", "red2")

#make the love plot for match_gbm_cluster1
love_plot <- love.plot(bal.tab(match_gbm_cluster1, m.threshold = 0.1),
                       stat = "mean.diffs",
                       grid = TRUE,
                       stars = "raw",
                       abs = FALSE,
                       colors = c("black", "red2")) + labs(title = NULL)

#make the bal.plot for distances
distances_plot <- bal.plot(match_gbm_cluster1, var.name = "distance", which = "both",
                           type = "histogram", mirror = TRUE) + scale_fill_manual(values = custom_colors) + labs(title = NULL)
#make the bal.plot for X1
x1_plot <- bal.plot(match_gbm_cluster1, var.name = 'X1', which = 'both', grid = TRUE)
x1_plot <- x1_plot + 
  xlab("School-level fixed mindset") +
  scale_fill_manual(values = custom_colors) + 
  labs(title = NULL)

#create the bal.plot for XC.0
xc0_plot <- bal.plot(match_gbm_cluster1, var.name = 'XC.0', which = 'both', grid = TRUE)
xc0_plot <- xc0_plot + 
  xlab("Students' urbanicity category: 'Rural'") +
  scale_fill_manual(values = custom_colors) +
  labs(title = NULL)

#arrange the plots in a grid
grid.arrange(love_plot, distances_plot, x1_plot, xc0_plot, ncol = 2)

#create matched df for best matching method (es.mean)
mdata1 <- match.data(match_gbm_cluster1, data = cluster1)
mdata2 <- match.data(match_gbm_cluster2, data = cluster2)
mdata3 <- match.data(match_gbm_cluster3, data = cluster3)
mdata4 <- match.data(match_gbm_cluster4, data = cluster4)
mdata5 <- match.data(match_gbm_cluster5, data = cluster5)
mdata6 <- match.data(match_gbm_cluster6, data = cluster6)
#repeat for benchmark model
mdata1.logit <- match.data(match_gbm_cluster1.logit, data = cluster1.logit)
mdata2.logit <- match.data(match_gbm_cluster2.logit, data = cluster2.logit)
mdata3.logit <- match.data(match_gbm_cluster3.logit, data = cluster3.logit)
mdata4.logit <- match.data(match_gbm_cluster4.logit, data = cluster4.logit)
mdata5.logit <- match.data(match_gbm_cluster5.logit, data = cluster5.logit)
mdata6.logit <- match.data(match_gbm_cluster6.logit, data = cluster6.logit)

#check multicollinearity
#correlation matrix
mdata1$schoolID <- as.numeric(mdata1$schoolID)
mdata1.cor <- cor(mdata1[,1:39])
corrplot(mdata1.cor)

##
##G-computation##
##

#estimate ATT using marginaleffects package (https://cran.r-project.org/web/packages/MatchIt/vignettes/estimating-effects.html)
fit.cluster1 <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.2 + C3.1),
                   data = mdata1)
ATT.cluster1 <- avg_comparisons(fit.cluster1, variables = "Z",
                                vcov = ~subclass,
                                newdata = subset(mdata1, Z == 1))
fit.cluster2 <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.2 + C3.1),
                   data = mdata2)
ATT.cluster2 <- avg_comparisons(fit.cluster2, variables = "Z",
                                vcov = ~subclass,
                                newdata = subset(mdata2, Z == 1))
fit.cluster3 <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.2 + C3.1),
                   data = mdata3)
ATT.cluster3 <- avg_comparisons(fit.cluster3, variables = "Z",
                                vcov = ~subclass,
                                newdata = subset(mdata3, Z == 1))
fit.cluster4 <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.2 + C3.1),
                   data = mdata4)
ATT.cluster4 <- avg_comparisons(fit.cluster4, variables = "Z",
                                vcov = ~subclass,
                                newdata = subset(mdata4, Z == 1))
fit.cluster5 <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.2 + C3.1),
                   data = mdata5)
ATT.cluster5 <- avg_comparisons(fit.cluster5, variables = "Z",
                                vcov = ~subclass,
                                newdata = subset(mdata5, Z == 1))
fit.cluster6 <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                            + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                            + C1.14 + C1.15 + C2.2 + C3.1),
                   data = mdata6)
ATT.cluster6 <- avg_comparisons(fit.cluster6, variables = "Z",
                                vcov = ~subclass,
                                newdata = subset(mdata6, Z == 1))

#do the same for benchmark model
fit.cluster1.logit <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.2 + C3.1),
                         data = mdata1.logit)
ATT.cluster1.logit <- avg_comparisons(fit.cluster1.logit, variables = "Z",
                                      vcov = ~subclass,
                                      newdata = subset(mdata1.logit, Z == 1))
fit.cluster2.logit <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.2 + C3.1),
                         data = mdata2.logit)
ATT.cluster2.logit <- avg_comparisons(fit.cluster2.logit, variables = "Z",
                                      vcov = ~subclass,
                                      newdata = subset(mdata2.logit, Z == 1))
fit.cluster3.logit <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.2 + C3.1),
                         data = mdata3.logit)
ATT.cluster3.logit<- avg_comparisons(fit.cluster3.logit, variables = "Z",
                                     vcov = ~subclass,
                                     newdata = subset(mdata3.logit, Z == 1))
fit.cluster4.logit <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.2 + C3.1),
                         data = mdata4.logit)
ATT.cluster4.logit <- avg_comparisons(fit.cluster4.logit, variables = "Z",
                                      vcov = ~subclass,
                                      newdata = subset(mdata4.logit, Z == 1))
fit.cluster5.logit <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.2 + C3.1),
                         data = mdata5.logit)
ATT.cluster5.logit <- avg_comparisons(fit.cluster5.logit, variables = "Z",
                                      vcov = ~subclass,
                                      newdata = subset(mdata5.logit, Z == 1))
fit.cluster6.logit <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                                  + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                                  + C1.14 + C1.15 + C2.2 + C3.1),
                         data = mdata6.logit)
ATT.cluster6.logit <- avg_comparisons(fit.cluster6.logit, variables = "Z",
                                      vcov = ~subclass,
                                      newdata = subset(mdata6.logit, Z == 1))
#do the same for the unadjusted sample
fit.cluster1.un <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.2 + C3.1),
                      data = cluster1.un)
ATT.cluster1.un <- avg_comparisons(fit.cluster1.un, variables = "Z",
                                   newdata = subset(cluster1.un, Z == 1))
fit.cluster2.un <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.2 + C3.1),
                      data = cluster2.un)
ATT.cluster2.un <- avg_comparisons(fit.cluster2.un, variables = "Z",
                                   newdata = subset(cluster2.un, Z == 1))
fit.cluster3.un <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.2 + C3.1),
                      data = cluster3.un)
ATT.cluster3.un <- avg_comparisons(fit.cluster3.un, variables = "Z",
                                   newdata = subset(cluster3.un, Z == 1))
fit.cluster4.un <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.2 + C3.1),
                      data = cluster4.un)
ATT.cluster4.un <- avg_comparisons(fit.cluster4.un, variables = "Z",
                                   newdata = subset(cluster4.un, Z == 1))
fit.cluster5.un <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.2 + C3.1),
                      data = cluster5.un)
ATT.cluster5.un <- avg_comparisons(fit.cluster5.un, variables = "Z",
                                   newdata = subset(cluster5.un, Z == 1))
fit.cluster6.un <- lm(Y ~ Z * (X1 + X2 + X3 + X4 + X5 + S3.2 + S3.3 + S3.4 + S3.5 + S3.6 + S3.7 
                               + C1.2 + C1.3 + C1.4 + C1.5 +C1.6 + C1.7 + C1.8 + C1.9 + C1.10 + C1.11 + C1.12 + C1.13
                               + C1.14 + C1.15 + C2.2 + C3.1),
                      data = cluster6.un)
ATT.cluster6.un <- avg_comparisons(fit.cluster6.un, variables = "Z",
                                   newdata = subset(cluster6.un, Z == 1))

##
##sensitivity analysis##
##

#gbm approach

#define a function to compute psens for each cluster
psens_calculation <- function(cluster, match_gbm, Gamma = 2, GammaInc = 0.1) {
  m.pairs <- cbind(cluster[row.names(match_gbm$match.matrix), 'Y'], 
                   cluster[match_gbm$match.matrix, 'Y']) #first column is the treated outcome and the second column is the controlled outcome (Qin-Yu Zhao, 2021)
  
  m.pairs <- m.pairs[complete.cases(m.pairs), ]
  x <- as.vector(m.pairs[, 1])
  y <- as.vector(m.pairs[, 2])
  
  psens_outcome <- psens(x = x, y = y, Gamma = Gamma, GammaInc = GammaInc)
  return(psens_outcome)
}

#call the function for each cluster (cluster and logit)
psens_cluster1 <- psens_calculation(cluster1, match_gbm_cluster1)
psens_cluster2 <- psens_calculation(cluster2, match_gbm_cluster2)
psens_cluster3 <- psens_calculation(cluster3, match_gbm_cluster3)
psens_cluster4 <- psens_calculation(cluster4, match_gbm_cluster4)
psens_cluster5 <- psens_calculation(cluster5, match_gbm_cluster5)
psens_cluster6 <- psens_calculation(cluster6, match_gbm_cluster6)

#for logit approach
psens_cluster1.logit <- psens_calculation(cluster1.logit, match_gbm_cluster1.logit, 2, 0.1)
psens_cluster2.logit <- psens_calculation(cluster2.logit, match_gbm_cluster2.logit, 2, 0.1)
psens_cluster3.logit <- psens_calculation(cluster3.logit, match_gbm_cluster3.logit, 2.2, 0.1)
psens_cluster4.logit <- psens_calculation(cluster4.logit, match_gbm_cluster4.logit, 2, 0.1)
psens_cluster5.logit <- psens_calculation(cluster5.logit, match_gbm_cluster5.logit, 2, 0.1)
psens_cluster6.logit <- psens_calculation(cluster6.logit, match_gbm_cluster6.logit, 2.2, 0.1)

##
##causal forest##
##

'Athey, S., & Wager, S. (2019). Estimating treatment effects with causal forests: An application. Observational studies, 5(2), 37-51.'

#set seed
set.seed(111)
#load data
data.all <- df_includingDescription
colnames(data.all)[1:13] <- c("schoolid", "Z", "Y", "S3", "C1", "C2", "C3", "XC", "X1", "X2", "X3", "X4", "X5")

#prepare data
data.all$schoolid = factor(data.all$schoolid)

df.causal = data.all[,-1]
school.id = as.numeric(data.all$schoolid)

school.mat = model.matrix(~ schoolid + 0, data = data.all)
school.size = colSums(school.mat)

#prepare causal forest input
W = df.causal$Z
Y = df.causal$Y
X.raw = df.causal[,-(1:2)]

#one-hot encoding
C1.enc = model.matrix(~ factor(X.raw$C1) + 0)
XC.enc = model.matrix(~ factor(X.raw$XC) + 0)

X = cbind(X.raw[,-which(names(X.raw) %in% c("C1", "XC"))], C1.enc, XC.enc)

#apply causal forest
Y.forest = regression_forest(X, Y, clusters = school.id, equalize.cluster.weights = TRUE)
Y.hat = predict(Y.forest)$predictions
W.forest = regression_forest(X, W, clusters = school.id, equalize.cluster.weights = TRUE)
W.hat = predict(W.forest)$predictions

cf.pre = causal_forest(X, Y, W,
                       Y.hat = Y.hat, W.hat = W.hat,
                       clusters = school.id,
                       equalize.cluster.weights = TRUE)
varimp = variable_importance(cf.pre)
selected.idx = which(varimp > mean(varimp))

cf = causal_forest(X[,selected.idx], Y, W,
                   Y.hat = Y.hat, W.hat = W.hat,
                   clusters = school.id,
                   equalize.cluster.weights = TRUE,
                   tune.parameters = "all")
tau.hat = predict(cf)$predictions


#estimate ATE
ATE = average_treatment_effect(cf)
paste("95% CI for the ATE:", round(ATE[1], 3),
      "+/-", round(qnorm(0.975) * ATE[2], 3))

#test heterogeneity
#run test calibration
test_calibration(cf)

#compare subgroups with high- and low estimated out-of-bag CATEs
high_CATE = tau.hat > median(tau.hat)
ATE.high = average_treatment_effect(cf, subset = high_CATE)
ATE.low = average_treatment_effect(cf, subset = !high_CATE)
paste("95% CI for difference in ATE:",
      round(ATE.high[1] - ATE.low[1], 3), "+/-",
      round(qnorm(0.975) * sqrt(ATE.high[2]^2 + ATE.low[2]^2), 3))
