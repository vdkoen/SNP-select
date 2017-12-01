library("ape")

snps_out <- "subset_freebayes_5_E3_snps.csv"
matrix_out <- "filtered_pairs_framboos_only_full_matrix.csv"
input_file <- "filtered_new_pairs_framboos_only_freebayes.csv"

setwd("/home/kdi/Filezilla_download/NA")


input_file <- "filtered_framboos_subset_freebayes.csv"

df <- read.table("/home/kdi/Documents/SNP_sets/framboos_named/snps_5_framboos_name_filtered_samtools.csv", sep = "\t", na.strings = "-", header = TRUE,row.names = 1)

pairs_df <-read.delim("pairs.txt", sep=",", header=FALSE, fill = TRUE) #read pairs

gene_dist <- dist.gene(df[-1,], method = "pairwise", pairwise.deletion = TRUE, variance = FALSE)

new_df <- df
new_df <- df[ order(row.names(df)), ]

#df <- new_df


f <- function(x){
  nm <- (unlist(unname(na.omit((x)))))
  new_df <<- new_df[!rownames(new_df) %in% nm[-1], ]
  new_name <<- paste(nm, collapse = "")
  row_index <<- match(nm[1], rownames(new_df))
  row.names(new_df)[row_index] <<- new_name
}

apply(pairs_df,1, f)


fc <- file("pairs.txt")
mylist <- strsplit(readLines(fc), ",")



for (pair in mylist){
  nm <- (unlist(pair)) 
  new_df <- new_df[!rownames(new_df) %in% nm[-1], ]
  new_name <- paste(nm, collapse = "")
  row_index <- match(nm[1], rownames(new_df))
  row.names(new_df)[row_index] <- new_name
}


new_df <- new_df[!rownames(new_df) %in% nm[-1], ]
new_name <- paste(nm, collapse = "")
row_index <- match(nm[1], rownames(new_df))
row.names(new_df)[row_index] <- new_name

rownames(new_df)
mylist[[1]][1]


















# df <- df[c("A1","B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E3", "F2", "G2", "H2"),]

df <- read.table("/home/kdi/Documents/SNP_sets/snps_5_filtered_framboos_subset_freebayes.csv", sep = "\t", na.strings = "-", header = TRUE,row.names = 1)
new_data <- read.table("/home/kdi/Documents/SNP_sets/test_snp_set.csv", sep = "\t", na.strings = "-", header = TRUE,row.names = 1)

test_my_set <- rbind(df, new_data)

names(df)
names(new_data)

gene_dist <- dist.gene(test_my_set[-1,], method = "pairwise", pairwise.deletion = TRUE, variance = FALSE)

gene_dist <- as.data.frame(gene_dist)

a <- which(gene_dist==0,arr.ind=T)
row.names(gene_dist[a[,1],] )
names(gene_dist[a[,2]]) 

gene_dist[which.min(gene_dist ==0 ),]
names(gene_dist[which.min(gene_dist),])



gene_dist <- as.matrix(gene_dist)
gene_dist[upper.tri(gene_dist)] <- NA
write.table(gene_dist, file = "results.csv", sep = "\t", na = "", col.names = NA)

hc <- hclust(gene_dist)
plot(hc)

