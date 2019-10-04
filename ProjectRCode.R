###############################################
####    Advanced Lab 5: Link Prediction    ####
####              Rahul Shah               ####
###############################################

# Setting Working Directory 
dir_path <- getwd()
setwd(dir_path)

# Clearing everything out of memory
rm(list=ls()) 

# Input file
infile_nasdaq <- "wide_twitter_daily_lab5.csv"

## Load package
library(igraph)
nasdaq = read.csv(infile_nasdaq, header = TRUE, sep = ",")
class(nasdaq)
# ---
# [1] "data.frame"
# ---

# Describe the data frame
str(nasdaq)

# Converting 'datestart' column to rownames
rownames(nasdaq) <- nasdaq[,1]
nasdaq[,1] <- NULL

############################
#### From Regular Lab 5 ####
############################
# CHUNK 5
heatmap(scale(nasdaq), Rowv=NA)

# Code of CHUNK 8 THROUGH CHUNK 13 uses the full correlation matrix to predict edges in the network
# CHUNK 8
mycorr <- cor(nasdaq)

# CHUNK 9
# Fisher's transformation
z <- 0.5 * log((1 + mycorr) / (1 - mycorr))

# CHUNK 10
z.vec <- z[upper.tri(z)]
n <- dim(nasdaq)[1]
corr.pvals <- 2 * pnorm(abs(z.vec), 0, 
                        sqrt(1 / (n-3)), lower.tail=FALSE)

# CHUNK 11
length(corr.pvals)
# ---
## [1] 4186
# ---

# Number of edges predicted: using statistical significance at the p < 0.05 threshold
length(corr.pvals[corr.pvals < 0.05])
# ---
## [1] 483
# ---

# CHUNK 12
# Benjamini-Hochberg adjustment to control for the false discovery rate
corr.pvals.adj <- p.adjust(corr.pvals, "BH")

# CHUNK 13
# Number of edges predicted: using  the overall correlation and statistical significance at the p < 0.05 threshold
corr.edges <- (corr.pvals.adj < 0.05)
length(corr.pvals.adj[corr.edges])
# ---
## [1] 180
# ---

summary(corr.pvals)
# ---
##    Min.  st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.2315  0.5164  0.4973  0.7665  0.9999  
# ---

summary(corr.pvals.adj)
# ---
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.9253  0.9989  0.8658  0.9989  0.9999 
# ---

# Creating the graph predicted by the statistically significant overall correlations
corr.A <- matrix(0, 92, 92)
corr.A[lower.tri(corr.A)] <- as.numeric(corr.edges)
g.corr <- graph.adjacency(corr.A, "undirected")
# plot(g.corr)


# CHUNK 14
library(fdrtool)

# CHUNK 15
mycorr.vec <- mycorr[upper.tri(mycorr)]
fdr <- fdrtool(mycorr.vec, statistic="correlation")

# Note the code of CHUNK 16 through CHUNK 19 uses partial correlations to predict edges
# CHUNK 16
pcorr.pvals <- matrix(0, dim(mycorr)[1], 
                      dim(mycorr)[2])
for(i in seq(1, 92)){
  for(j in seq(1, 92)){
    rowi <- mycorr[i, -c(i, j)]
    rowj <- mycorr[j, -c(i, j)]
    tmp <- (mycorr[i, j] - 
              rowi*rowj)/sqrt((1-rowi^2) * (1-rowj^2))
    tmp.zvals <- (0.5) * log((1+tmp) / (1-tmp))
    tmp.s.zvals <- sqrt(n-4) * tmp.zvals
    tmp.pvals <- 2 * pnorm(abs(tmp.s.zvals), 
                           0, 1, lower.tail=FALSE)
    pcorr.pvals[i, j] <- max(tmp.pvals)
  }
}

# CHUNK 17
pcorr.pvals.vec <- pcorr.pvals[lower.tri(pcorr.pvals)]
# Benjamini-Hochberg adjustment to control for the false discovery rate
pcorr.pvals.adj <- p.adjust(pcorr.pvals.vec, "BH")

# CHUNK 18
# Number of edges predicted: using  the partial correlation and statistical significance at the p < 0.05 threshold
pcorr.edges <- (pcorr.pvals.adj < 0.05)
length(pcorr.pvals.adj[pcorr.edges])
# ---
## [1] 48
# ---

# CHUNK 19
# Creating the graph predicted by the statistically significant partial correlations (p < 0.05)
pcorr.A <- matrix(0, 92, 92)
# length(lower.tri(pcorr.A))
# length(as.numeric(pcorr.edges))
pcorr.A[lower.tri(pcorr.A)] <- as.numeric(pcorr.edges)
g.pcorr <- graph.adjacency(pcorr.A, "undirected")
# plot(g.pcorr)

# Number of edges predicted: using  the partial correlation and statistical significance at the p < 0.01 threshold
pcorr.edges.01 <- (pcorr.pvals.adj < 0.01)
length(pcorr.pvals.adj[pcorr.edges.01])
# ---
## [1] 32
# ---

# Creating the graph predicted by the statistically significant partial correlations (p < 0.01)
pcorr.A.01 <- matrix(0, 92, 92)
pcorr.A.01[lower.tri(pcorr.A.01)] <- as.numeric(pcorr.edges.01)
g.pcorr.01 <- graph.adjacency(pcorr.A.01, "undirected")
# plot(g.pcorr.01)

# Find overlap between the two graphs
graph.intersection(g.pcorr, g.pcorr.01, byname=FALSE)
# ---
## IGRAPH 915bccc U--- 92 32 -- 
## + edges from 915bccc:
##  [1] 87--92 72--88 65--85 62--77 61--72 57--89 52--71 46--54 43--77 43--62 34--88 31--52 28--32 26--81 23--29 22--30
## [17] 20--90 20--54 20--27 19--64 19--61 16--64 15--31 12--60 10--70  9--72  9--61  9--28  7--69  3--38  3-- 9  1--31
# ---

# CHUNK 21
# FDR tool can also be used to adjust for false discovery rate and predict new edges based on partial correlations
fdr <- fdrtool(pcorr.pvals.vec, statistic="pvalue", 
               plot=FALSE)
pcorr.edges.2 <- (fdr$qval < 0.05)
length(fdr$qval[pcorr.edges.2])
# ---
## [1] 47
# ---

# Creating the graph predicted by the statistically significant partial correlations of FDR tool
pcorr.A.2 <- matrix(0, 92, 92)
pcorr.A.2[lower.tri(pcorr.A.2)] <- as.numeric(pcorr.edges.2)
g.pcorr.2 <- graph.adjacency(pcorr.A.2, "undirected")
# plot(g.pcorr.2)

# Find overlap between the graph predicted by partial correlations of FDR tool with the partial correlations of BH
graph.intersection(g.pcorr, g.pcorr.2, byname=FALSE)
# ---
## IGRAPH b003926 U--- 92 47 -- 
## + edges from b003926:
##  [1] 87--92 80--83 78--79 72--88 65--85 62--77 61--72 57--89 52--71 50--63 46--54 43--77 43--62 39--80 34--88
## [16] 32--76 32--72 32--37 31--52 28--32 26--81 24--53 24--38 23--77 23--61 23--29 22--78 22--30 20--90 20--54
## [31] 20--27 19--64 19--61 16--64 15--31 12--60 10--70 10--20  9--72  9--61  9--47  9--28  7--69  7--59  3--38
## [46]  3-- 9  1--31
# ---

# HUGE (High-dimensional undirected graph estimation library) procedure for predicted links
# CHUNK 22
library(huge)
set.seed(42)
# Converting the data frame into a matrix for HUGE
nasdaq_mat <- data.matrix(nasdaq,rownames.force = NA)
huge.out <- huge(nasdaq_mat)

# CHUNK 23
huge.opt <- huge.select(huge.out, criterion="ric")
summary(huge.opt$refit)
# ---
## 92 x 92 sparse Matrix of class "dsCMatrix", with 0 entries 
## [1] i j x
## <0 rows> (or 0-length row.names)
# ---

# CHUNK 24
huge.opt <- huge.select(huge.out, criterion="stars")
g.huge <- graph.adjacency(huge.opt$refit, "undirected")
# plot(g.huge, vertex.size=3, vertex.label=NA)
summary(g.huge)
# ---
## IGRAPH 65dcc01 U--- 92 109 --
# ---

# CHUNK 25
# Find overlap between the graph produced by the partial adjusted correlations (BH) with graph produced by HUGE library
graph.intersection(g.pcorr, g.huge)
# ---
## IGRAPH 8a89f75 U--- 92 48 -- 
## + edges from 8a89f75:
##  [1] 87--92 80--83 78--79 72--88 65--85 62--77 61--72 57--89 52--71 50--63 46--54 43--77 43--62 39--80 34--88
## [16] 32--76 32--72 32--37 31--52 28--32 26--81 24--53 24--38 23--77 23--61 23--29 22--78 22--30 20--90 20--54
## [31] 20--27 19--64 19--61 16--64 15--31 12--60 10--70 10--20  9--72  9--64  9--61  9--47  9--28  7--69  7--59
## [46]  3--38  3-- 9  1--31
# ---

# CHUNK 26
# Find overlap between the graph produced by the partial adjusted correlations (FDR tools) with graph produced by HUGE library
graph.intersection(g.pcorr.2, g.huge, byname=FALSE)
# ---
## IGRAPH e6835ab U--- 92 47 -- 
## + edges from e6835ab:
##  [1] 87--92 80--83 78--79 72--88 65--85 62--77 61--72 57--89 52--71 50--63 46--54 43--77 43--62 39--80 34--88
## [16] 32--76 32--72 32--37 31--52 28--32 26--81 24--53 24--38 23--77 23--61 23--29 22--78 22--30 20--90 20--54
## [31] 20--27 19--64 19--61 16--64 15--31 12--60 10--70 10--20  9--72  9--61  9--47  9--28  7--69  7--59  3--38
## [46]  3-- 9  1--31
# ---


##################
### Viola Plot ###
##################

nv <- vcount(g.corr)
ncn <- numeric()
A<-get.adjacency(g.corr)

nv <- vcount(g.pcorr)
ncn <- numeric()
A<-get.adjacency(g.pcorr)

nv <- vcount(g.pcorr.01)
ncn <- numeric()
A<-get.adjacency(g.pcorr.01)

nv <- vcount(g.pcorr.2)
ncn <- numeric()
A<-get.adjacency(g.pcorr.2)

nv <- vcount(g.huge)
ncn <- numeric()
A <- huge.opt$refit

# Find the number of common neighbors for each pair of nodes in the g.huge network
for(i in (1:(nv-1))){
  ni <- neighborhood(g.huge, 1, i)
  nj <- neighborhood(g.huge, 1, (i+1):nv)
  nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
  temp <- unlist(lapply(nbhd.ij, length)) - 
    2*A[i, (i+1):nv]
  ncn <- c(ncn, temp)
}

library(vioplot)
Avec <- A[lower.tri(A)]
vioplot(ncn[Avec==0], ncn[Avec==1], 
        names=c("No Edge", "Edge"), col='purple')
title(ylab="Number of Common Neighbors")



#######################
### Network Plot: 1 ###
#######################
# Benjamini-Hochberg adjustment used to control for the false discovery rate;
# and a threshold of p < 0.05 used to identify statistically significant overall
# correlations. 

graph.corrL <- set.vertex.attribute(g.corr, "name", value=names(nasdaq))
graph.corr = delete.vertices(graph.corrL, which(degree(graph.corrL) < 1))

degree(graph.corr)
table(degree(graph.corr))
neighbors(graph.corr, v=c('EXPE'))

set.seed(9143)  # setting seed so as to make the layout reproducible
layout1 <- layout.fruchterman.reingold(graph.corr, niter=500) # Creating a layout object
# layout2 <- layout.kamada.kawai(graph.corr)

# Node Size and Color
# V(graph.corr)$size <- 12 

# V(graph.corr)$color <- "grey"
V(graph.corr)[degree(graph.corr) <= 2]$color <- "lightcyan"
V(graph.corr)[degree(graph.corr) == 3]$color <- "cyan"
V(graph.corr)[degree(graph.corr) == 4]$color <- "cyan3"
V(graph.corr)[degree(graph.corr) == 5]$color <- "mediumblue"
V(graph.corr)[degree(graph.corr) >= 6]$color <- "navyblue"

# V(graph.corr)$label.dist <- ifelse(V(graph.corr)$size >= 10, 0, 1.2)
V(graph.corr)$label.color <- 'black'
# V(graph.corr)$label.size <- 0.1

# Edge width and Color
E(graph.corr)$width <- 4
E(graph.corr)$color <- "gold3"

#Plotting 
plot(graph.corr, layout=layout1, vertex.size=12,
     vertex.label.cex = 0.9, vertex.label.degree = 2, vertex.label=names(nasdaq))
legend("left", inset=1.0, title="Node Color Scheme",
       c("Deg <= 2","Deg = 3","Deg = 4","Deg = 5","Deg >= 6"), 
       fill=c('lightcyan','cyan','cyan3','mediumblue', 'navyblue'), horiz=FALSE, box.lty=0, cex=0.7)
# title("Nasdaq100 Network Graph: Fruchterman Reingold Layout")


#######################
### Network Plot: 2 ###
#######################
# Benjamini-Hochberg adjustment used to control for the false discovery rate;
# and a threshold of p < 0.05 used to identify statistically significant partial
# correlations. 

graph.pcorrL <- set.vertex.attribute(g.pcorr, "name", value=names(nasdaq))
graph.pcorr = delete.vertices(graph.pcorrL, which(degree(graph.pcorrL) < 1))

degree(graph.pcorr)
table(degree(graph.pcorr))
neighbors(graph.pcorr, v=c('AMZN'))

set.seed(9143)  # setting seed so as to make the layout reproducible
layout1 <- layout.fruchterman.reingold(graph.pcorr, niter=500) # Creating a layout object
# layout2 <- layout.kamada.kawai(graph.pcorr)

# Node Size and Color
# V(graph.pcorr)$size <- 12 

# V(graph.pcorr)$color <- "grey"
V(graph.pcorr)[degree(graph.pcorr) == 1]$color <- "lightcyan"
V(graph.pcorr)[degree(graph.pcorr) == 2]$color <- "cyan"
V(graph.pcorr)[degree(graph.pcorr) == 3]$color <- "cyan3"
V(graph.pcorr)[degree(graph.pcorr) == 4]$color <- "mediumblue"
V(graph.pcorr)[degree(graph.pcorr) == 6]$color <- "navyblue"

# V(graph.pcorr)$label.dist <- ifelse(V(graph.pcorr)$size >= 10, 0, 1.2)
V(graph.pcorr)$label.color <- 'black'
# V(graph.pcorr)$label.size <- 0.1

# Edge width and Color
E(graph.pcorr)$width <- 4
E(graph.pcorr)$color <- "gold3"

#Plotting 
plot(graph.pcorr, layout=layout1, vertex.size=12,
     vertex.label.cex = 0.9, vertex.label.degree = 2)
legend("left", inset=1.0, title="Node Color Scheme",
       c("Deg = 1","Deg = 2","Deg = 3","Deg = 4","Deg = 6"), 
       fill=c('lightcyan','cyan','cyan3','mediumblue', 'navyblue'), horiz=FALSE, box.lty=0, cex=0.7)
# title("Nasdaq100 Network Graph: Fruchterman Reingold Layout")

#######################
### Network Plot: 3 ###
#######################
# Benjamini-Hochberg adjustment used to control for the false discovery rate;
# and a threshold of p < 0.01 used to identify statistically significant partial
# correlations. 

graph.pcorr.01L <- set.vertex.attribute(g.pcorr.01, "name", value=names(nasdaq))
graph.pcorr.01 = delete.vertices(graph.pcorr.01L, which(degree(graph.pcorr.01L) < 1))

degree(graph.pcorr.01)
table(degree(graph.pcorr.01))
neighbors(graph.pcorr.01, v=c('AMZN'))

set.seed(9143)  # setting seed so as to make the layout reproducible
layout1 <- layout.fruchterman.reingold(graph.pcorr.01, niter=500) # Creating a layout object
# layout2 <- layout.kamada.kawai(graph.pcorr.01)

# Node Size and Color
# V(graph.pcorr.01)$size <- 12 

# V(graph.pcorr.01)$color <- "grey"
V(graph.pcorr.01)[degree(graph.pcorr.01) == 1]$color <- "lightcyan"
V(graph.pcorr.01)[degree(graph.pcorr.01) == 2]$color <- "cyan"
V(graph.pcorr.01)[degree(graph.pcorr.01) == 3]$color <- "cyan3"
V(graph.pcorr.01)[degree(graph.pcorr.01) == 4]$color <- "mediumblue"
# V(graph.pcorr.01)[degree(graph.pcorr.01) == 6]$color <- "navyblue"

# V(graph.pcorr.01)$label.dist <- ifelse(V(graph.pcorr.01)$size >= 10, 0, 1.2)
V(graph.pcorr.01)$label.color <- 'black'
# V(graph.pcorr.01)$label.size <- 0.1

# Edge width and Color
E(graph.pcorr.01)$width <- 4
E(graph.pcorr.01)$color <- "gold3"

#Plotting 
plot(graph.pcorr.01, layout=layout1, vertex.size=12,
     vertex.label.cex = 0.9, vertex.label.degree = 2)
legend("left", inset=1.0, title="Node Color Scheme",
       c("Deg = 1","Deg = 2","Deg = 3","Deg = 4"), 
       fill=c('lightcyan','cyan','cyan3','mediumblue'), horiz=FALSE, box.lty=0, cex=0.7)
# title("Nasdaq100 Network Graph: Fruchterman Reingold Layout")


#######################
### Network Plot: 4 ###
#######################
# FDR tool library used to adjust for false discovery rate; and threshold of 
# p < 0.05 used to identify statistically significant partial correlations.

graph.pcorr.2L <- set.vertex.attribute(g.pcorr.2, "name", value=names(nasdaq))
graph.pcorr.2 = delete.vertices(graph.pcorr.2L, which(degree(graph.pcorr.2L) < 1))

degree(graph.pcorr.2)
table(degree(graph.pcorr.2))
neighbors(graph.pcorr.2, v=c('AMZN'))

set.seed(9143)  # setting seed so as to make the layout reproducible
layout1 <- layout.fruchterman.reingold(graph.pcorr.2, niter=500) # Creating a layout object
# layout2 <- layout.kamada.kawai(graph.pcorr.2)

# Node Size and Color
# V(graph.pcorr.2)$size <- 12 

# V(graph.pcorr.2)$color <- "grey"
V(graph.pcorr.2)[degree(graph.pcorr.2) == 1]$color <- "lightcyan"
V(graph.pcorr.2)[degree(graph.pcorr.2) == 2]$color <- "cyan"
V(graph.pcorr.2)[degree(graph.pcorr.2) == 3]$color <- "cyan3"
V(graph.pcorr.2)[degree(graph.pcorr.2) == 4]$color <- "mediumblue"
V(graph.pcorr.2)[degree(graph.pcorr.2) == 5]$color <- "navyblue"

# V(graph.pcorr.2)$label.dist <- ifelse(V(graph.pcorr.2)$size >= 10, 0, 1.2)
V(graph.pcorr.2)$label.color <- 'black'
# V(graph.pcorr.2)$label.size <- 0.1

# Edge width and Color
E(graph.pcorr.2)$width <- 4
E(graph.pcorr.2)$color <- "gold3"

#Plotting 
plot(graph.pcorr.2, layout=layout1, vertex.size=12,
     vertex.label.cex = 0.9, vertex.label.degree = 2)
legend("left", inset=1.0, title="Node Color Scheme",
       c("Deg = 1","Deg = 2","Deg = 3","Deg = 4","Deg = 5"), 
       fill=c('lightcyan','cyan','cyan3','mediumblue', 'navyblue'), horiz=FALSE, box.lty=0, cex=0.7)
# title("Nasdaq100 Network Graph: Fruchterman Reingold Layout")


#######################
### Network Plot: 5 ###
#######################
# HUGE procedure used for predicted links.

graph.hugeL <- set.vertex.attribute(g.huge, "name", value=names(nasdaq))
graph.huge = delete.vertices(graph.hugeL, which(degree(graph.hugeL) < 1))

degree(graph.huge)
table(degree(graph.huge))
neighbors(graph.huge, v=c('AMZN'))

set.seed(9143)  # setting seed so as to make the layout reproducible
layout1 <- layout.fruchterman.reingold(graph.huge, niter=500) # Creating a layout object
# layout2 <- layout.kamada.kawai(graph.huge)

# Node Size and Color
# V(graph.huge)$size <- 12 

# V(graph.huge)$color <- "grey"
V(graph.huge)[degree(graph.huge) == 1]$color <- "lightcyan"
V(graph.huge)[degree(graph.huge) == 2]$color <- "cyan"
V(graph.huge)[degree(graph.huge) == 3]$color <- "cyan3"
V(graph.huge)[degree(graph.huge) == 4]$color <- "mediumblue"
V(graph.huge)[degree(graph.huge) >= 5]$color <- "navyblue"

# V(graph.huge)$label.dist <- ifelse(V(graph.huge)$size >= 10, 0, 1.2)
V(graph.huge)$label.color <- 'black'
# V(graph.huge)$label.size <- 0.1

# Edge width and Color
E(graph.huge)$width <- 4
E(graph.huge)$color <- "gold3"

#Plotting 
plot(graph.huge, layout=layout1, vertex.size=12,
     vertex.label.cex = 0.9, vertex.label.degree = 2)
legend("left", inset=1.0, title="Node Color Scheme",
       c("Deg = 1","Deg = 2","Deg = 3","Deg = 4","Deg >= 5"), 
       fill=c('lightcyan','cyan','cyan3','mediumblue', 'navyblue'), horiz=FALSE, box.lty=0, cex=0.7)
# title("Nasdaq100 Network Graph: Fruchterman Reingold Layout")
