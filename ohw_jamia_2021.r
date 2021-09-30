
options(stringsAsFactors = FALSE);


# ===========================================================================
# parameters
# 

obs_hours <- 24;
obs_res <- 1.5;
n_lvs <- 10;
n_clusts <- 4;



# ===========================================================================
# Data preparation
# 

# ---------------------------------------------------------------------------
# development and validation set
# 

r.all <- 1:nrow(d.cohort);
set.seed(0);
r.dev <- sort(sample(r.all, round(nrow(d.cohort)*.8), replace=FALSE));
r.val <- setdiff(r.all, r.dev);

d.dev <- data.frame(
    MRN  = NA,
    rid  = r.dev,
    type = "dev"
);

d.val <- data.frame(
    MRN  = NA,
    rid  = r.val,
    type = "val"
);

d.eval <- rbind(d.dev, d.val);
d.eval <- dplyr::arrange(d.eval, rid);
d.eval <- dplyr::transmute(d.eval,
    MRN = d.cohort$MRN,
    type = factor(x=type, levels=c("dev", "val"))
);

table(d.eval$type);
## dev val 
## 829 207 



# ---------------------------------------------------------------------------
# sequence
# - Create sequences consisting of 16 non-overlapping interval slots
# 

seq_breaks <- seq(from=0, to=obs_hours, by=obs_res);
seq_breaks <- seq_breaks[-length(seq_breaks)];
seq_n_ticks <- length(seq_breaks) - 1;

# key: cartesian product of MRN and elapse.
d.tp.key <- data.frame(
    MRN    = rep(d.cohort$MRN, each=length(seq_breaks)),
    elapse = rep(seq_breaks, times=length(d.cohort$MRN))
);

# value: aggregate biomarkers and treatments.
d.tp.value <- dplyr::mutate(d.tp,
    elapse = as.numeric(as.character(cut(x=elapse, breaks=c(seq_breaks, obs_hours), labels=seq_breaks, right=FALSE)))
);
d.tp.value <- dplyr::summarise_all(dplyr::group_by(d.tp.value, MRN, elapse), max);
d.tp.value <- as.data.frame(d.tp.value);

# combine key and value
d.tp <- dplyr::left_join(
    x  = d.tp.key,
    y  = d.tp.value,
    by = c("MRN", "elapse")
);

d.tp <- dplyr::mutate_at(d.tp, 
    dplyr::vars(-MRN, -elapse), 
    function(x) { 
        ifelse(test=!is.na(x), yes=x*1L, no=0);
});
d.tp <- dplyr::arrange(d.tp, MRN, elapse);



# ---------------------------------------------------------------------------
# dimensionality reduction - entire
# 

# * Parameter tuning for for logistic PCA *
set.seed(0);
fit.tp.lpca.cv <- logisticPCA::cv.lpca(
    x = dplyr::select(d.tp, -MRN, -elapse), 
    ks = 1:(ncol(d.tp) - 2), 
    ms = seq(2, 20, by = 2), 
    folds = 5, 
    quiet = FALSE
);

png(file="fig/construct/dm/fig.tp.lpca.cv.png", width=6, height=6, units="in", res=600);
plot(fit.tp.lpca.cv);
dev.off();


# * Logistic PCA *
set.seed(0);
fit.tp.lpca <- logisticPCA::logisticPCA(
    x = dplyr::select(d.tp, -MRN, -elapse), 
    k = 10, 
    m = 8
);

d.tp.lpca <- predict(fit.tp.lpca, dplyr::select(d.tp, -MRN, -elapse));

d.tp.lpca <- data.frame(
    dplyr::select(d.tp, MRN, elapse),
    d.tp.lpca
);


# * Hierarchical clustering  *
# - Construct hierarchical clustering using logistic principal components
# - To make 10 levels.

set.seed(0);
nb.1days.all <- NbClust::NbClust(
    data     = dplyr::select(d.tp.lpca, -MRN, -elapse), 
    distance = "euclidean", 
    min.nc   = 2,
    max.nc   = 25, 
    method   = "ward.D2",
    index    = "all"
);

dist.tp.lpca = dist(dplyr::select(d.tp.lpca, -MRN, -elapse), method="euclidean");
set.seed(0);
fit.tp.hc.lpca = hclust(dist.tp.lpca, method="ward.D2");

png(file="fig/construct/dm/fig.tp.hc.lpca.png", width=6, height=6, units="in", res=600);
plot(fit.tp.hc.lpca, labels=FALSE);
dev.off();

# Average silhouette width
d.tp.hc.lpca.sil <- sapply(2:100, 
    function(k, hc, dist) {
        sils  <- cluster::silhouette(cutree(hc, k=k), dist);
        sil_m <- mean(sils[, "sil_width"]);
        
        return( sil_m );
}, hc=fit.tp.hc.lpca, dist=dist.tp.lpca);
d.tp.hc.lpca.sil <- c(0, d.tp.hc.lpca.sil);
d.tp.hc.lpca.sil <- data.frame(
    x = 1:length(d.tp.hc.lpca.sil),
    y = d.tp.hc.lpca.sil
);

g.tp.hc.lpca.sil <- ggplot(d.tp.hc.lpca.sil, aes(x=x, y=y)) +
    geom_line() +
    geom_point() +
    labs(x="Number of clusters k", y="Average silhouette width") +
    theme_classic()
ggsave(filename="fig/construct/dm/fig.tp.hc.lpca.silhouette.png", plot=g.tp.hc.lpca.sil, device="png", width=6, height=6, units="in", dpi=600);

# Elbow method
png(file="fig/construct/dm/fig.tp.hc.lpca.elbow.png", width=6, height=6, units="in", res=600);
factoextra::fviz_nbclust(dplyr::select(d.tp.lpca, -MRN, -elapse), FUN=factoextra::hcut, method="wss")
dev.off();


# * Final sequence *
d.tp.dr <- data.frame(
    dplyr::select(d.tp.lpca, MRN, elapse),
    cluster = cutree(fit.tp.hc.lpca, k=n_lvs)
);


# ---------------------------------------------------------------------------
# dimensionality reduction - development
# 

set.seed(0);
fit.tp.dev.lpca <- logisticPCA::logisticPCA(
    x = dplyr::select(d.tp.dev, -MRN, -elapse),
    k = 10,
    m = 8
);

d.tp.dev.lpca <- predict(fit.tp.dev.lpca, dplyr::select(d.tp.dev, -MRN, -elapse));

d.tp.dev.lpca <- data.frame(
    dplyr::select(d.tp.dev, MRN, elapse),
    d.tp.dev.lpca
);

dist.tp.dev.lpca = dist(dplyr::select(d.tp.dev.lpca, -MRN, -elapse), method="euclidean");
set.seed(0);
fit.tp.dev.hc.lpca = hclust(dist.tp.dev.lpca, method="ward.D2");

d.tp.dev.dr <- data.frame(
    dplyr::select(d.tp.dev.lpca, MRN, elapse),
    cluster = cutree(fit.tp.dev.hc.lpca, k=n_lvs)
);



# ---------------------------------------------------------------------------
# dimensionality reduction - validation
# 

set.seed(0);
fit.tp.val.lpca <- logisticPCA::logisticPCA(
    x = dplyr::select(d.tp.val, -MRN, -elapse),
    k = 10,
    m = 8
);

d.tp.val.lpca <- predict(fit.tp.val.lpca, dplyr::select(d.tp.val, -MRN, -elapse));

d.tp.val.lpca <- data.frame(
    dplyr::select(d.tp.val, MRN, elapse),
    d.tp.val.lpca
);

dist.tp.val.lpca = dist(dplyr::select(d.tp.val.lpca, -MRN, -elapse), method="euclidean");
set.seed(0);
fit.tp.val.hc.lpca = hclust(dist.tp.val.lpca, method="ward.D2");

d.tp.val.dr <- data.frame(
    dplyr::select(d.tp.val.lpca, MRN, elapse),
    cluster = cutree(fit.tp.val.hc.lpca, k=n_lvs)
);



# ---------------------------------------------------------------------------
# Convert the sequences into strings
# - LV takes strings.
# 

# Entire
d.string <- dplyr::group_by(d.tp.dr, MRN);
d.string <- dplyr::summarise(d.string,
    string = paste(LETTERS[cluster], collapse="")
);
d.string <- as.data.frame(d.string);

# Development
d.string.dev <- dplyr::group_by(d.tp.dev.dr, MRN);
d.string.dev <- dplyr::summarise(d.string.dev,
    string = paste(LETTERS[cluster], collapse="")
);
d.string.dev <- as.data.frame(d.string.dev);

# Validation
d.string.val <- dplyr::group_by(d.tp.val.dr, MRN);
d.string.val <- dplyr::summarise(d.string.val,
    string = paste(LETTERS[cluster], collapse="")
);
d.string.val <- as.data.frame(d.string.val);





# ===========================================================================
# Subphenotype derivation
# 

# ---------------------------------------------------------------------------
# Levenshtein distance
# 

edist <- function(data, method="lv", weight=c(d=1, i=1, s=1, t=1))
{
    dist <- stringdist::stringdistmatrix(a=data, b=data, method=method, weight=weight);
    dist <- as.dist(dist);
    
    return(dist);
}



# ---------------------------------------------------------------------------
# hierarchical clustering method
# 

# Cohort - entire
dist.string <- edist(d.string$string, method="lv");
set.seed(0);
fit.hclust.string = hclust(dist.string, method="ward.D2");


# Cohort - development
dist.string.dev <- edist(d.string.dev$string, method="lv");
set.seed(0);
fit.hclust.string.dev = hclust(dist.string.dev, method="ward.D2");


# Cohort - validation
dist.string.val <- edist(d.string.val$string, method="lv");
set.seed(0);
fit.hclust.string.val = hclust(dist.string.val, method="ward.D2");


# Elbow method
library(ggplot2);

set.seed(0);
g.hclust.string.elbow <- factoextra::fviz_nbclust(
        x             = as.matrix(dist.string), 
        FUNcluster    = factoextra::hcut, 
        method        = "wss", 
        diss          = dist.string,
        k.max         = 15,
        print.summary = FALSE) +
    labs(y="Total within-clusters sum of squares");

png(file="fig/cluster/hclust/fig.hclust.string.elbow.png", width=5, height=5, units="in", res=600);
g.hclust.string.elbow;
dev.off();


# Silhouette method
library(ggplot2);

set.seed(0);
g.hclust.string.silhouette <- factoextra::fviz_nbclust(
        x             = as.matrix(dist.string), 
        FUNcluster    = factoextra::hcut, 
        method        = "silhouette", 
        diss          = dist.string,
        k.max         = 15,
        print.summary = FALSE) +
    labs(y="Average silhouette width");

png(file="fig/cluster/hclust/fig.hclust.string.silhouette.png", width=5, height=5, units="in", res=600);
g.hclust.string.silhouette;
dev.off();


# Gap statistic
library(ggplot2);

set.seed(0);
g.hclust.string.gap <- factoextra::fviz_nbclust(
        x             = as.matrix(dist.string), 
        FUNcluster    = factoextra::hcut, 
        method        = "gap_stat", 
        diss          = dist.string,
        k.max         = 15,
        nstart        = 25,
        nboot         = 50, 
        print.summary = TRUE,
        maxSE         = list(method="Tibs2001SEmax", SE.factor=1)) +
    labs(y="Gap statistics (k)");

png(file="fig/cluster/hclust/fig.hclust.string.gap.png", width=5, height=5, units="in", res=600);
g.hclust.string.gap;
dev.off();


# Clest method
library(ggplot2);

set.seed(0);
g.hclust.string.clest <- factoextra::fviz_nbclust(
        x             = as.matrix(dist.string), 
        FUNcluster    = factoextra::hcut, 
        method        = "gap_stat", 
        diss          = dist.string,
        k.max         = 15,
        nstart        = 25,
        nboot         = 50, 
        print.summary = TRUE,
        maxSE         = list(method="globalSEmax", SE.factor=3)) +
    labs(y="Gap statistics (k)");

png(file="fig/cluster/hclust/fig.hclust.string.clest.png", width=5, height=5, units="in", res=600);
g.hclust.string.clest;
dev.off();



# ---------------------------------------------------------------------------
# k-means clustering method
# 

# Cohort - entire
dist.string <- edist(d.string$string, method="lv");
set.seed(0);
fit.kmeans.string <- kmeans(dist.string, centers=n_clusts);


# Cohort - development
dist.string.dev <- edist(d.string.dev$string, method="lv");
set.seed(0);
fit.kmeans.string.dev <- kmeans(dist.string.dev, centers=n_clusts);


# Cohort - validation
dist.string.val <- edist(d.string.val$string, method="lv");
set.seed(0);
fit.kmeans.string.val <- kmeans(dist.string.val, centers=n_clusts);





# ===========================================================================
# Evaluation
# 

# ---------------------------------------------------------------------------
# Hierarchical clustering
# 

k     <- n_clusts;
clust <- cutree(fit.hclust.string, k=k);
clust <- factor(x=clust, levels=1:k, labels=as.roman(1:k));

d.string.hclust <- data.frame(
    MRN   = d.string$MRN,
    clust = clust
);
rm(k, clust);



# ---------------------------------------------------------------------------
# K-means clustering
# 

k     <- n_clusts
clust <- fit.kmeans.string$cluster;
clust <- factor(x=clust, levels=1:k, labels=as.roman(1:k));

d.string.kmeans <- data.frame(
    MRN   = d.string$MRN,
    clust = clust
);
rm(k, clust);



# ---------------------------------------------------------------------------
# Dendrogram
# 

library(dendextend);

fit.dend.string <- as.dendrogram(fit.hclust.string);
fit.dend.string <- color_branches(fit.dend.string, k=n_clusts, col=colorspace::rainbow_hcl(n=4, c=90, l=50)[c(4, 2, 1, 3)]);

fit.dend.string <- set(
    dend  = fit.dend.string, 
    what  = "labels", 
    value = rep("", nrow(d.string))
);

png(file="fig/analysis/hclust/fig.hclust.string.dend.png", width=5, height=5, units="in", res=600);
plot(fit.dend.string, horiz=FALSE);
dev.off();



# ---------------------------------------------------------------------------
# Cumulative biomarker and treatment patterns
# 

# Common parameters

l.bx <- c(
    "sbp.90m", "hr.90p", "rr.20p", "paco2.32m", "wbc.12kp", 
    "p_f.150m", "glu.250p", "ph.7_25m", "ag.12p", "cr.past7_2xp"
);
l.tx <- c("rbt", "va", "di", "mv", "nb", "ins", "hd");

l.bx.full <- c(
    "SBP < 90 mm Hg", "HR > 90 BPS", "RR > 20 BPS", "PaCO2 < 32 mm Hg", "WBC > 12 K", 
    "P/F < 150", 
    "Glucose > 250 mg/dL", "pH < 7.25", "AG > 12 mmol/L", 
    "7 days Cr > 2X"
);
l.tx.full <- c(
    "Blood transfusion", "Vasoactive agents", 
    "Loop diuretics", 
    "Mechanical ventilation", "Neuromuscular Blockers", 
    "Insulin", 
    "Hemodialysis"
);

l.ft <- c(l.bx, l.tx);

seq_breaks <- seq(from=0, to=obs_hours, by=obs_res);
seq_breaks <- seq_breaks[-length(seq_breaks)];

library(ggplot2);

d.tp.bx <- dplyr::select_at(d.tp, .vars=c("MRN", "elapse", l.bx));
d.tp.bx <- lapply(seq_breaks, function(s, data, clust) {
    data <- dplyr::filter(data, elapse <= s);
    
    # 
    data <- dplyr::group_by(data, MRN);
    data <- dplyr::summarise_all(data, max);
    data <- as.data.frame(data);
    
    # 
    data <- dplyr::left_join(
        x = clust,
        y = data,
        by = c("MRN")
    );
    data <- dplyr::select(data, -MRN);
    data <- dplyr::group_by(data, clust);
    data <- dplyr::summarise_all(data, mean);
    data <- as.data.frame(data);
    data <- dplyr::arrange(data, clust);
    
    return(data);
}, data=d.tp.bx, clust=d.string.hclust);
d.tp.bx <- do.call(rbind, d.tp.bx);
d.tp.bx <- tidyr::gather_(d.tp.bx, "type", "value", l.bx);

d.tp.tx <- dplyr::select_at(d.tp, .vars=c("MRN", "elapse", l.tx));
d.tp.tx <- lapply(seq_breaks, function(s, data, clust) {
    data <- dplyr::filter(data, elapse <= s);
    
    # 
    data <- dplyr::group_by(data, MRN);
    data <- dplyr::summarise_all(data, max);
    data <- as.data.frame(data);
    
    # 
    data <- dplyr::left_join(
        x = clust,
        y = data,
        by = c("MRN")
    );
    data <- dplyr::select(data, -MRN);
    data <- dplyr::group_by(data, clust);
    data <- dplyr::summarise_all(data, mean);
    data <- as.data.frame(data);
    data <- dplyr::arrange(data, clust);
    
    return(data);
}, data=d.tp.tx, clust=d.string.hclust);
d.tp.tx <- do.call(rbind, d.tp.tx);
d.tp.tx <- tidyr::gather_(d.tp.tx, "type", "value", l.tx);

d.tp.btx <- rbind(d.tp.bx, d.tp.tx);
d.tp.btx <- dplyr::mutate(d.tp.btx,
    type = factor(x=type, levels=c(l.bx, l.tx), labels=c(l.bx.full, l.tx.full))
);

g.tp.btx.cum <- ggplot(d.tp.btx, aes(x=elapse, y=value, shape=clust, colour=clust)) +
    geom_point() + 
    geom_line() +
    facet_wrap(~ type, ncol=5) +
    scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
    labs(x="Hours", y="Cumulative prevalence", shape="Subphenotype", colour="Subphenotype") +
    theme_light() +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/hclust/fig.hclust.string.tp.btx.cum.png", plot=g.tp.btx.cum, device="png", width=8, height=7, units="in", dpi=600);



# ---------------------------------------------------------------------------
# 30 days in-hospital mortality
# 

library(ggplot2);
library(survival);
library(survminer);

fit.string.km.mort <- survfit(Surv(mort_tm, mort) ~ clust, data=data.frame(d.mort30, d.string.hclust));
png(file="fig/analysis/hclust/fig.hclust.string.km.mort.png", width=5, height=5, units="in", res=600);
ggsurvplot(
  fit          = fit.string.km.mort,
  conf.int     = TRUE, # Add confidence interval
  # pval         = TRUE, # Add p-value
  ggtheme      = theme_classic(), # Change ggplot2 theme
  legend       = "bottom",
  legend.title = "Subphenotype",
  legend.labs  = as.roman(1:n_clusts)
);
dev.off();



# ---------------------------------------------------------------------------
# 30 days discharge
# 

library(ggplot2);
library(survival);
library(survminer);

fit.string.km.disch <- survfit(Surv(disch_tm, disch) ~ clust, data=data.frame(d.disch30, d.string.hclust));
png(file="fig/analysis/hclust/fig.hclust.string.km.disch.png", width=5, height=5, units="in", res=600);
ggsurvplot(
  fit          = fit.string.km.disch,
  fun          = "event",
  conf.int     = TRUE, # Add confidence interval
  # pval         = TRUE, # Add p-value
  ggtheme      = theme_classic(), # Change ggplot2 theme
  ylim         = c(0, 1),
  legend       = "bottom",
  legend.title = "Subphenotype",
  legend.labs  = as.roman(1:n_clusts)
)
dev.off();



# ---------------------------------------------------------------------------
# Monthy trend
# 

library(ggplot2);

d.monthly <- data.frame(
    month = data.table::month(d.cohort$icu_admit_dt)
);
d.monthly <- dplyr::mutate(d.monthly,
    month = factor(x=month, levels=sort(unique(month)))
);
d.monthly <- dplyr::group_by(d.monthly, month);
d.monthly <- dplyr::summarise(d.monthly, 
    n = dplyr::n()
);
d.monthly <- as.data.frame(d.monthly);

g.monthly <- ggplot(d.monthly, aes(x=month, y=n)) +
    geom_bar(stat="identity") +
    labs(x=NULL, y="Frequency") +
    theme_classic() + 
    theme(axis.text.x=element_blank());
ggsave(filename="fig/analysis/hclust/fig.hclust.string.monthly.entire.png", plot=g.monthly, device="png", width=5, height=1.5, units="in", dpi=600);


d.monthly.clusters <- data.frame(
    month = data.table::month(d.cohort$icu_admit_dt), 
    class = d.string.hclust$clust
);
d.monthly.clusters <- dplyr::mutate(d.monthly.clusters,
    month = factor(x=month, levels=sort(unique(month)))
);

d.monthly.clusters <- dplyr::group_by(d.monthly.clusters, month);
d.monthly.clusters <- dplyr::mutate(d.monthly.clusters, 
    tot_month = dplyr::n()
);
d.monthly.clusters <- as.data.frame(d.monthly.clusters);

d.monthly.clusters <- dplyr::group_by(d.monthly.clusters, month, class);
d.monthly.clusters <- dplyr::summarise(d.monthly.clusters,
    perc = dplyr::n() / tot_month[1]
);
d.monthly.clusters <- as.data.frame(d.monthly.clusters);

g.monthly.clusters <- ggplot(d.monthly.clusters, aes(x=month, y=perc, fill=class)) +
    geom_bar(stat="identity") +
    labs(x="Month", y="Proportion", fill="Subphenotype") +
    theme_classic() +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/hclust/fig.hclust.string.monthly.clusters.png", plot=g.monthly.clusters, device="png", width=5, height=3.5, units="in", dpi=600);



# ---------------------------------------------------------------------------
# Sensitiveness of clustering algorithms
# 

library(ggplot2);

# K-means clustering
k     <- n_clusts;
clust <- fit.kmeans.string$cluster;
clust <- factor(x=clust, levels=c(4,2,1,3), labels=as.roman(1:k));

d.string.kmeans <- data.frame(
    MRN   = d.string$MRN,
    clust = clust
);
rm(k, clust);

d.clust <- table(hclust=d.string.hclust$clust, kmeans=d.string.kmeans$clust) / as.vector(table(d.string.hclust$clust))
d.clust <- as.data.frame(d.clust);

g.clust <- ggplot(d.clust, aes(x=kmeans, y=hclust, fill=Freq)) +
    geom_tile() +
    scale_fill_gradient(name="Proportion", low="white", high="red", limits=c(0,1)) +
    labs(x="Rederived subphenotypes", y="Subphenotypes") +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/hclust/fig.robustness.algorithms.png", plot=g.clust, device="png", width=5, height=5, units="in", dpi=600);



# ---------------------------------------------------------------------------
# Sensitiveness of sampling
# 

library(ggplot2);

# Development
k     <- n_clusts;
clust <- cutree(fit.hclust.string.dev, k=k);
clust <- factor(x=clust, levels=1:k, labels=as.roman(1:k));

d.string.hclust.dev <- data.frame(
    MRN   = d.string.dev$MRN,
    clust = clust
);
rm(k, clust);

# Validation
k     <- n_clusts;
clust <- cutree(fit.hclust.string.val, k=k);
clust <- factor(x=clust, levels=c(2,1,4,3), labels=as.roman(1:k));

d.string.hclust.val <- data.frame(
    MRN   = d.string.val$MRN,
    clust = clust
);
rm(k, clust);

# k-nearest neighbor
knn.predict <- function(x, y, newdata)
{
    classes <- sapply(newdata, function(n, x, y) {
        dist <- stringdist::stringdist(a=x, b=n, method="lv");
        class <- y[ which.min(dist) ];
        return(class);
    }, x=x, y=y);
    names(classes) <- NULL;
    return(classes);
}

yhat.val <- knn.predict(x=d.string.dev$string, y=d.string.hclust.dev$clust, newdata=d.string.val$string);

d.sampling <- table(y=d.string.hclust.val$clust, yhat=yhat.val) / as.vector(table(d.string.hclust.val$clust))
d.sampling <- as.data.frame(d.sampling);

g.sampling <- ggplot(d.sampling, aes(x=yhat, y=y, fill=Freq)) +
    geom_tile() +
    scale_fill_gradient(name="Proportion", low="white", high="red", limits=c(0,1)) +
    labs(x="Rederived subphenotypes", y="Subphenotypes") +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/hclust/fig.robustness.sampling.png", plot=g.sampling, device="png", width=5, height=5, units="in", dpi=600);





# ===========================================================================
# Appendix
# 

# ---------------------------------------------------------------------------
# Flatten data
# 

d.fl <- dplyr::group_by(d.tp, MRN);
d.fl <- dplyr::summarise_at(d.fl, dplyr::vars(sbp.90m:hd), max);
d.fl <- as.data.frame(d.fl);
d.fs <- data.frame(
    MRN = dplyr::select(d.fl, MRN),
    string = apply(dplyr::select(d.fl, -MRN), 1, paste, collapse="")
);



# ---------------------------------------------------------------------------
# Cosine distance
# 

edist <- function(data, method="lv", weight=c(d=1, i=1, s=1, t=1))
{
    dist <- stringdist::stringdistmatrix(a=data, b=data, method=method, weight=weight);
    dist <- as.dist(dist);
    
    return(dist);
}

dist.fl.cosine <- edist(d.fs$string, method="cosine");
set.seed(0);
fit.hclust.fl.cosine = hclust(dist.fl.cosine, method="ward.D2");

k     <- n_clusts;
clust <- cutree(fit.hclust.fl.cosine, k=k);
clust <- factor(x=clust, levels=1:k, labels=as.roman(1:k));

d.hclust.fl.cosine <- data.frame(
    MRN   = d.string$MRN,
    clust = clust
);
rm(k, clust);



# ---------------------------------------------------------------------------
# Jaccard distance
# 

dist.fl.jaccard <- philentropy::distance(d.fl[, -1], method="jaccard");
dist.fl.jaccard <- as.dist(dist.fl.jaccard);
set.seed(0);
fit.hclust.fl.jaccard = hclust(dist.fl.jaccard, method="ward.D2");

k     <- n_clusts;
clust <- cutree(fit.hclust.fl.jaccard, k=k);
clust <- factor(x=clust, levels=1:k, labels=as.roman(1:k));

d.hclust.fl.jaccard <- data.frame(
    MRN   = d.string$MRN,
    clust = clust
);
rm(k, clust);



# ---------------------------------------------------------------------------
# Dendrogram
# 

# Cosine distance
library(dendextend);

fit.dend.fl.cosine <- as.dendrogram(fit.hclust.fl.cosine);
fit.dend.fl.cosine <- color_branches(fit.dend.fl.cosine, k=n_clusts, col=colorspace::rainbow_hcl(n=4, c=90, l=50)[c(4, 3, 2, 1)]);

fit.dend.fl.cosine <- set(
    dend  = fit.dend.fl.cosine, 
    what  = "labels", 
    value = rep("", nrow(d.hclust.fl.cosine))
);

png(file="fig/analysis/flat/fig.hclust.fl.cosine.dend.png", width=5, height=5, units="in", res=600);
plot(fit.dend.fl.cosine, horiz=FALSE);
dev.off();


# Jaccard distance
library(dendextend);

fit.dend.fl.jaccard <- as.dendrogram(fit.hclust.fl.jaccard);
fit.dend.fl.jaccard <- color_branches(fit.dend.fl.jaccard, k=n_clusts, col=colorspace::rainbow_hcl(n=4, c=90, l=50)[c(1, 4, 3, 2)]);

fit.dend.fl.jaccard <- set(
    dend  = fit.dend.fl.jaccard, 
    what  = "labels", 
    value = rep("", nrow(d.hclust.fl.jaccard))
);

png(file="fig/analysis/flat/fig.hclust.fl.jaccard.dend.png", width=5, height=5, units="in", res=600);
plot(fit.dend.fl.jaccard, horiz=FALSE);
dev.off();


# ---------------------------------------------------------------------------
# Biomarkers and treatments patterns
# 

# Common parameters
l.bx <- c(
    "sbp.90m", "hr.90p", "rr.20p", "paco2.32m", "wbc.12kp", 
    "p_f.150m", "glu.250p", "ph.7_25m", "ag.12p", "cr.past7_2xp"
);
l.tx <- c("rbt", "va", "di", "mv", "nb", "ins", "hd");

l.bx.full <- c(
    "SBP < 90 mm Hg", "HR > 90 BPS", "RR > 20 BPS", "PaCO2 < 32 mm Hg", "WBC > 12 K", 
    "P/F < 150", 
    "Glucose > 250 mg/dL", "pH < 7.25", "AG > 12 mmol/L", 
    "7 days Cr > 2X"
);
l.tx.full <- c(
    "Blood transfusion", "Vasoactive agents", 
    "Loop diuretics", 
    "Mechanical ventilation", "Neuromuscular Blockers", 
    "Insulin", 
    "Hemodialysis"
);

l.ft <- c(l.bx, l.tx);
l.ft.full <- c(l.bx.full, l.tx.full);

seq_breaks <- seq(from=0, to=obs_hours, by=obs_res);
seq_breaks <- seq_breaks[-length(seq_breaks)];


# Cosine distance
library(ggplot2);

d.tp.bx <- dplyr::select_at(d.tp, .vars=c("MRN", "elapse", l.bx));
d.tp.bx <- lapply(seq_breaks, function(s, data, clust) {
    data <- dplyr::filter(data, elapse <= s);
    
    # 
    data <- dplyr::group_by(data, MRN);
    data <- dplyr::summarise_all(data, max);
    data <- as.data.frame(data);
    
    # 
    data <- dplyr::left_join(
        x = clust,
        y = data,
        by = c("MRN")
    );
    data <- dplyr::select(data, -MRN);
    data <- dplyr::group_by(data, clust);
    data <- dplyr::summarise_all(data, mean);
    data <- as.data.frame(data);
    data <- dplyr::arrange(data, clust);
    
    return(data);
}, data=d.tp.bx, clust=d.hclust.fl.cosine);
d.tp.bx <- do.call(rbind, d.tp.bx);
d.tp.bx <- tidyr::gather_(d.tp.bx, "type", "value", l.bx);

d.tp.tx <- dplyr::select_at(d.tp, .vars=c("MRN", "elapse", l.tx));
d.tp.tx <- lapply(seq_breaks, function(s, data, clust) {
    data <- dplyr::filter(data, elapse <= s);
    
    # 
    data <- dplyr::group_by(data, MRN);
    data <- dplyr::summarise_all(data, max);
    data <- as.data.frame(data);
    
    # 
    data <- dplyr::left_join(
        x = clust,
        y = data,
        by = c("MRN")
    );
    data <- dplyr::select(data, -MRN);
    data <- dplyr::group_by(data, clust);
    data <- dplyr::summarise_all(data, mean);
    data <- as.data.frame(data);
    data <- dplyr::arrange(data, clust);
    
    return(data);
}, data=d.tp.tx, clust=d.hclust.fl.cosine);
d.tp.tx <- do.call(rbind, d.tp.tx);
d.tp.tx <- tidyr::gather_(d.tp.tx, "type", "value", l.tx);

d.tp.btx <- rbind(d.tp.bx, d.tp.tx);
d.tp.btx <- dplyr::mutate(d.tp.btx,
    type = factor(x=type, levels=c(l.bx, l.tx), labels=c(l.bx.full, l.tx.full))
);

g.fl.cosine.tp.btx.cum <- ggplot(d.tp.btx, aes(x=elapse, y=value, shape=clust, colour=clust)) +
    geom_point() + 
    geom_line() +
    facet_wrap(~ type, ncol=5) +
    scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
    labs(x="Hours", y="Cumulative prevalence", shape="Subphenotype", colour="Subphenotype") +
    theme_light() +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/flat/fig.hclust.fl.cosine.tp.btx.cum.png", plot=g.fl.cosine.tp.btx.cum, device="png", width=8, height=7, units="in", dpi=600);


# Jaccard distance
library(ggplot2);

d.tp.bx <- dplyr::select_at(d.tp, .vars=c("MRN", "elapse", l.bx));
d.tp.bx <- lapply(seq_breaks, function(s, data, clust) {
    data <- dplyr::filter(data, elapse <= s);
    
    # 
    data <- dplyr::group_by(data, MRN);
    data <- dplyr::summarise_all(data, max);
    data <- as.data.frame(data);
    
    # 
    data <- dplyr::left_join(
        x = clust,
        y = data,
        by = c("MRN")
    );
    data <- dplyr::select(data, -MRN);
    data <- dplyr::group_by(data, clust);
    data <- dplyr::summarise_all(data, mean);
    data <- as.data.frame(data);
    data <- dplyr::arrange(data, clust);
    
    return(data);
}, data=d.tp.bx, clust=d.hclust.fl.jaccard);
d.tp.bx <- do.call(rbind, d.tp.bx);
d.tp.bx <- tidyr::gather_(d.tp.bx, "type", "value", l.bx);

d.tp.tx <- dplyr::select_at(d.tp, .vars=c("MRN", "elapse", l.tx));
d.tp.tx <- lapply(seq_breaks, function(s, data, clust) {
    data <- dplyr::filter(data, elapse <= s);
    
    # 
    data <- dplyr::group_by(data, MRN);
    data <- dplyr::summarise_all(data, max);
    data <- as.data.frame(data);
    
    # 
    data <- dplyr::left_join(
        x = clust,
        y = data,
        by = c("MRN")
    );
    data <- dplyr::select(data, -MRN);
    data <- dplyr::group_by(data, clust);
    data <- dplyr::summarise_all(data, mean);
    data <- as.data.frame(data);
    data <- dplyr::arrange(data, clust);
    
    return(data);
}, data=d.tp.tx, clust=d.hclust.fl.jaccard);
d.tp.tx <- do.call(rbind, d.tp.tx);
d.tp.tx <- tidyr::gather_(d.tp.tx, "type", "value", l.tx);

d.tp.btx <- rbind(d.tp.bx, d.tp.tx);
d.tp.btx <- dplyr::mutate(d.tp.btx,
    type = factor(x=type, levels=c(l.bx, l.tx), labels=c(l.bx.full, l.tx.full))
);

g.fl.jaccard.tp.btx.cum <- ggplot(d.tp.btx, aes(x=elapse, y=value, shape=clust, colour=clust)) +
    geom_point() + 
    geom_line() +
    facet_wrap(~ type, ncol=5) +
    scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
    labs(x="Hours", y="Cumulative prevalence", shape="Subphenotype", colour="Subphenotype") +
    theme_light() +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/flat/fig.hclust.fl.jaccard.tp.btx.cum.png", plot=g.fl.jaccard.tp.btx.cum, device="png", width=8, height=7, units="in", dpi=600);



# ---------------------------------------------------------------------------
# Comparison
# 

# Cosine distance
d.clust <- table(lv=d.string.hclust$clust, cosine=d.hclust.fl.cosine$clust) / as.vector(table(d.string.hclust$clust))
d.clust <- as.data.frame(d.clust);

g.clust <- ggplot(d.clust, aes(x=cosine, y=lv, fill=Freq)) +
    geom_tile() +
    scale_fill_gradient(name="Proportion", low="white", high="red", limits=c(0,1)) +
    labs(x="Cosine", y="Levenshtein") +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/flat/fig.comparison.lv_cosine.png", plot=g.clust, device="png", width=5, height=5, units="in", dpi=600);


# Jaccard distance
d.clust <- table(lv=d.string.hclust$clust, jaccard=d.hclust.fl.jaccard$clust) / as.vector(table(d.string.hclust$clust))
d.clust <- as.data.frame(d.clust);

g.clust <- ggplot(d.clust, aes(x=jaccard, y=lv, fill=Freq)) +
    geom_tile() +
    scale_fill_gradient(name="Proportion", low="white", high="red", limits=c(0,1)) +
    labs(x="Jaccard", y="Levenshtein") +
    theme(legend.position="bottom");
ggsave(filename="fig/analysis/flat/fig.comparison.lv_jaccard.png", plot=g.clust, device="png", width=5, height=5, units="in", dpi=600);









