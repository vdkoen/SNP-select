library("ape")
library("ggfortify")
library('ggplot2')
#TODO is to add a argument parser which alles to run the script from the command line. 

#TODO add a input that says the set has to be atleast 



input_file <- "merged_5_snps.csv"
setwd("/home/kdi/Filezilla_download/")

df <- read.table(input_file, sep = "\t", na.strings = "-", header = TRUE,row.names = 1)
subset_df <- read.table("/home/kdi/Documents/SNP_sets/snps_5_freebayes_filtered.csv", sep = "\t", na.strings = "-", header = TRUE,row.names = 1)

gene_dist <- dist.gene(df[-1,], method = "pairwise", pairwise.deletion = TRUE, variance = FALSE)
gene_dist_set <- dist.gene(subset_df[-1,], method = "pairwise", pairwise.deletion = TRUE, variance = FALSE)

sum(gene_dist < 1)

#MDS test.
fit <- cmdscale(gene_dist, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y, pch = 16)
text(x, y, pos = 3, labels = rownames(subset_df[-1,]), cex = 0.7)

#Dendogram
hc <- hclust(as.dist(M), method = "average")                # apply hirarchical clustering
plot(hc)

plot(as.phylo(hc), type = "unrooted")
plot(as.phylo(hc), type = "fan", tip.color = hsv(runif(15, 0.65, 
0.95), 1, 1, 0.7), edge.color = hsv(runif(10, 0.65, 0.75), 1, 1, 0.7), edge.width = runif(20, 0.5, 3), use.edge.length = TRUE, col = "gray80")

colors = c("red", "blue", "green", "black")
clus4 = cutree(hc, 24)
plot(as.phylo(hc), type = "fan", tip.color = colors[clus4],
     label.offset = 1, cex = 0.7)






#fancy MDS plot
xy <- data.frame(cmdscale(d), factor(rownames(M)))
names(xy) <- c("x", "y", "cluster")
xy$model <- rownames(xy)
plot1 <- ggplot(xy, aes(x, y)) + geom_point(aes(colour=cluster), size=3) + geom_text(aes(label=cluster),hjust=1.5, vjust=0, size = 3.5) + theme(legend.position="none")

xy <- data.frame(cmdscale(gene_dist), factor(rownames(df[-1,])))
names(xy) <- c("x", "y", "cluster")
xy$model <- rownames(xy)
plot2 <- ggplot(xy, aes(x, y)) + geom_point(aes(colour=cluster), size=3) + geom_text(aes(label=cluster),hjust=1.2, vjust=0.2, size = 3.5) + theme(legend.position="none")

require(gridExtra)
grid.arrange(plot1, plot2, ncol=2)



#mantel test, could be used to test if performance of distance object is similar, if p < 0.005 it is otherwise it is not and the performance is much different p > 0.05
mantel.test(as.matrix(gene_dist), as.matrix(gene_dist_set), graph = TRUE, nperm = 999)



# some useless tree
result <- pvclust(as.matrix(gene_dist), method.dist="uncentered", method.hclust="average", nboot=10)
plot(result)
pvrect(result)


#----------------------------------------------------------------------------------------------------------

for(i in 1:ncol(all_samples[-1,])){
  all_samples[-1,][is.na(all_samples[-1,][,i]), i] <- mean(all_samples[-1,][,i], na.rm = TRUE)
}

test.pca <- prcomp(na.omit(all_samples[-1,]), center = TRUE, scale = TRUE)
plot(test.pca$x[,1:2])
plot(test.pca$x[,1],test.pca$x[,1])

autoplot(prcomp(all_samples[-1,]),label = TRUE)
