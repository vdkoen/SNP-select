library("ape")

# setwd("/home/kdi/Filezilla_download/GRASSEN/algorithm")

arg_options <- commandArgs(trailingOnly = TRUE)
arg_options

dir.create(file.path(arg_options[2]))
file_name<- basename(arg_options[1])
# bed_name <- strsplit(file_name, ".", fixed = TRUE)

snps_file <- paste("snps_", arg_options[3],"_", file_name, sep = "") # Make output string
matrix_file <- paste("matrix_", arg_options[3],"_", file_name, sep = "")

# bed_file <- paste("bed_", arg_options[3], "_", bed_name[[1]][1],".bed", sep = "")
snps_out <- file.path(arg_options[2], snps_file) # Make output dir and name in string.
matrix_out <- file.path(arg_options[2], matrix_file)
input_file <- arg_options[1]

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# snps_out <- "/home/kdi/Filezilla_download/GRASSEN/algorithm/test_snps_3_test.csv" # Make output dir and name in string.
# matrix_out <- "/home/kdi/Filezilla_download/GRASSEN/algorithm/test_matrix_3_test.csv"
# input_file <- "/home/kdi/Pictures/allel_frequencies_grassen.csv"


df <- read.table(input_file, sep = "\t", na.strings = "-", header = TRUE, row.names = 1)


# if a pair file is specified merge the pairs. 
if (length(arg_options) > 3 ){
  fc <- file(arg_options[4])
  mylist <- strsplit(readLines(fc), ",")
  
  
  for (pair in mylist){
    nm <- (unlist(pair)) 
    df <- df[!rownames(df) %in% nm[-1], ]
    new_name <- paste(nm, collapse = "")
    row_index <- match(nm[1], rownames(df))
    row.names(df)[row_index] <- new_name
  }
  rownames(df)
}

# value = 1 # if one than it makes minimal snp set, if 3 their have to be atleast 3 snp's that seperate the groups.
# 
value <- as.numeric(arg_options[3])

output_files <- function(final_markers, gene_dist){
  write.table(final_markers, file = snps_out, sep = "\t", na = "-", row.names = TRUE, col.names = NA)
  gene_dist <- as.matrix(gene_dist)
  gene_dist[upper.tri(gene_dist)] <- NA
  write.table(gene_dist, file = matrix_out, sep = "\t", na = "", col.names = NA)
}

snp_selection <- c(1)# snp to start with. the one with the highest entropy

previous_value <- 0
previous_check <- 0
check = 0.2

while(check != value + 1){
  null_list <- c()
  snp_number <- c()
  
  for (i in 1:length(df)){
    if (any(snp_selection==i) == FALSE){
      
      selection <- df[-1 ,c(snp_selection,i)] #new dataframe with extra snp(i)
      gene_dist <- dist(selection) #calc distance with new snp
      gene_dist[is.na(gene_dist)] <- 0
      
      
      if (sum(gene_dist < value) == 0){ # if there are no 0's in the pairwise comparson they are all unique, meaning we found a snp set.
        cat("Minimal snp_set found! ")
        cat("Containing", length(snp_selection) + 1, "snps", "\n")
        snp_selection <<- c(snp_selection, i)
        final_markers <<- df[,c(snp_selection)]
        output_files(final_markers, gene_dist)
        # stop_quietly()
        quit()
      }
      
      
      null_list <- c(null_list, sum(gene_dist < check)) #save results of the additional snp(i), this is the amount of 0's
      snp_number <- c(snp_number, i) # save snp numbers, so we can relate it back to its performance
    }
  }
  
  # if ( (sum(gene_dist < check) >= previous_value) & (check == previous_check) ){
  #   print("No optimal soluion found")
  #   cat("Closest solution contains", length(snp_selection), "snps", "\n")
  #   cat("Lowest distance: ", min(gene_dist), "\n")
  #   final_markers <<- df[,c(snp_selection)]
  #   output_files(final_markers, gene_dist)
  #   quit()
  #   # stop_quietly()
  # }
  
  previous_value <- sum(gene_dist < check)
  previous_check <- check
  
  if (length(unique(null_list)) != 1){ #skip a non informative snp selection.
    best_hit <- which.min(null_list) # get the best hit with highest maf score
    snp_selection <- c(snp_selection,snp_number[best_hit]) # add the best hit to the snp set and iterate again. 
  }
  
  print(sum(gene_dist < check))
  
  if (sum(gene_dist < check) == 0){
    check <- check + 0.2
  }
  
}


# hc <- hclust(gene_dist)                # apply hirarchical clustering 
# plot(hc)

#gene_dist
#sum(gene_dist)


