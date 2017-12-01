library("ape")

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

snps_out <- "entropy_5_snps.csv"
matrix_out <- "entropy_5_matrix.csv"
input_file <- "filtered_framboos_subset_freebayes.csv"

setwd("/home/kdi/Filezilla_download/NA")
df <- read.table(input_file, sep = "\t", na.strings = "-", header = TRUE,row.names = 1)
df <- df[c("A1","B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E3", "F2", "G2", "H2"),]

pairs_df <-read.table("pairs.txt", sep=",", header=FALSE, fill = TRUE)



new_df <- df#[-1, ]
new_df <- new_df[ order(row.names(new_df)), ]
length(row.names(new_df))

df <- new_df

value = 1 # if one than it makes minimal snp set, if 3 their have to be atleast 3 snp's that seperate the groups.
# value <- arg_options[3]

output_files <- function(final_markers, gene_dist){
  write.table(final_markers, file = snps_out, sep = "\t", na = "-", row.names = TRUE, col.names = NA)
  gene_dist <- as.matrix(gene_dist)
  gene_dist[upper.tri(gene_dist)] <- NA
  write.table(gene_dist, file = matrix_out, sep = "\t", na = "", col.names = NA)
}

#make small set for testing the algorithm.
#df <- df[, 1:700]

#TODO improvement can be made by not just selecting the first two but rather implement 
#a method which selects the starting ones that contain lots of information already

snp_selection <- c(3)# snp to start with.

previous_value <- 0
previous_check <- 0
check = 1


#for (x in 1:length(df)){
  
while(check != value + 1){
  null_list <- c()
  snp_number <- c()

  for (i in 1:length(df)){
    if (any(snp_selection==i) == FALSE){

      selection <- df[-1 ,c(snp_selection,i)] #new dataframe with extra snp(i)
      gene_dist <- dist.gene(selection, method = "pairwise", pairwise.deletion = TRUE, variance = FALSE) #calc distance with new snp

      if (sum(gene_dist < value) == 0){ # if there are no 0's in the pairwise comparson they are all unique, meaning we found a snp set.
        cat("Minimal snp_set found! ")
        cat("Containing", length(snp_selection) + 1, "snps", "\n")
        snp_selection <<- c(snp_selection, i)
        final_markers <<- df[,c(snp_selection)]
        output_files(final_markers, gene_dist)
        stop_quietly()
      }
      
      null_list <- c(null_list, sum(gene_dist < check)) #save results of the additional snp(i), this is the amount of 0's
      snp_number <- c(snp_number, i) # save snp numbers, so we can relate it back to its performance
    }
  }
  
  if ( (sum(gene_dist < check) >= previous_value) & (check == previous_check) ){
    print("No optimal soluion found")
    cat("Closest solution contains", length(snp_selection) + 1, "snps", "\n")
    final_markers <<- df[,c(snp_selection)]
    output_files(final_markers, gene_dist)
    stop_quietly()
  }
  
  previous_value <- sum(gene_dist < check)
  previous_check <- check
  
  if (length(unique(null_list)) != 1){ #skip a non informative snp selection.
    best_hit <- which.min(null_list) # get the best hit with highest maf score
    snp_selection <- c(snp_selection,snp_number[best_hit]) # add the best hit to the snp set and iterate again. 
    
  }
  
  print(sum(gene_dist < check))
  
  if (sum(gene_dist < check) == 0){
    check <- check + 1
    # print(check)
  }
  
}


hc <- hclust(gene_dist)                # apply hirarchical clustering
plot(hc)

#gene_dist
#sum(gene_dist)

bed <- names(df[,c(snp_selection)])[[1]]
splitted_bed <- strsplit(bed, ".", fixed = TRUE)
as.integer(splitted_bed[[1]][2]) + 100
as.integer(splitted_bed[[1]][2]) - 100



bed_dat <- c()

write_bed <-function(snp_selection){
  
  for (name in names(df[,c(snp_selection)])){
    splitted_bed <- strsplit(name, ".", fixed = TRUE)
    start <- as.integer(splitted_bed[[1]][2]) - 100
    stop <- as.integer(splitted_bed[[1]][2]) + 100
    
    bed_dat <- c(bed_dat, splitted_bed[[1]][1], start, stop)
    
  } 
  bed_file <- matrix(bed_dat, ncol = 3, byrow = TRUE)
  write.table(bed_file, file = "test_bed_out.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FAL)
}



