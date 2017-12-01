library("ape")

arg_options <- commandArgs(trailingOnly = TRUE)
arg_options

# stop_quietly <- function() {
#   opt <- options(show.error.messages = FALSE)
#   on.exit(options(opt))
#   stop()
# }

dir.create(file.path(arg_options[2]))
file_name<- basename(arg_options[1])
bed_name <- strsplit(file_name, ".", fixed = TRUE)

snps_file <- paste("snps_", arg_options[3],"_", file_name, sep = "") # Make output string
matrix_file <- paste("matrix_", arg_options[3],"_", file_name, sep = "")
# bed_file <- paste("bed_", arg_options[3], "_", bed_name[[1]][1],".bed", sep = "")

# bed_out <- file.path(arg_options[2], bed_file)
snps_out <- file.path(arg_options[2], snps_file) # Make output dir and name in string.
matrix_out <- file.path(arg_options[2], matrix_file)
input_file <- arg_options[1]

# cat("It needs 3 inputs: input table, output folder, and indication of the minimal snp set (integer) where 1 is the smallest set.", "\n")
# cat("example: Rscript configure_snp_set.R input_file output_folder 1", '\n')
# stop_quietly()

df <- read.table(input_file, sep = "\t", na.strings = "-", header = TRUE,row.names = 1)

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
value <- as.numeric(arg_options[3])

output_files <- function(final_markers, gene_dist){
  write.table(final_markers, file = snps_out, sep = "\t", na = "-", row.names = TRUE, col.names = NA, quote = FALSE)
  gene_dist <- as.matrix(gene_dist)
  gene_dist[upper.tri(gene_dist)] <- NA
  write.table(gene_dist, file = matrix_out, sep = "\t", na = "", col.names = NA, quote = FALSE)
}

#
#
# write_bed <- function(snp_selection){
#
#   bed_dat <- c()
#
#   for (name in names(df[,c(snp_selection)])){
#     splitted_bed <- strsplit(name, ".", fixed = TRUE)
#     start <- as.integer(splitted_bed[[1]][2]) - 100
#     stop <- as.integer(splitted_bed[[1]][2]) + 100
#
#     bed_dat <- c(bed_dat, splitted_bed[[1]][1], start, stop)
#
#   }
#   bed_matrix <- matrix(bed_dat, ncol = 3, byrow = TRUE)
#   write.table(bed_matrix, file = bed_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# }

snp_selection <- c(1)# snp to start with. the one with the highest entropy

previous_value <- 0
previous_check <- 0
check = 1

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
        # stop_quietly()
        quit()
      }

      
      null_list <- c(null_list, sum(gene_dist < check)) #save results of the additional snp(i), this is the amount of 0's
      snp_number <- c(snp_number, i) # save snp numbers, so we can relate it back to its performance
    }
  }
  
  if ( (sum(gene_dist < check) >= previous_value) & (check == previous_check) ){
    print("No optimal soluion found")
    cat("Closest solution contains", length(snp_selection), "snps", "\n")
    final_markers <<- df[,c(snp_selection)]
    output_files(final_markers, gene_dist)
    quit()
    # stop_quietly()
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
  }
  
}


# hc <- hclust(gene_dist)                # apply hirarchical clustering 
# plot(hc)

#gene_dist
#sum(gene_dist)


