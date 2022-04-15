#Load data
msa <- read.csv("https://raw.githubusercontent.com/pine-bio-support/COVID-19-origin/main/SARS_CoV_2_1.csv", header=TRUE)

#Define selected positions (Here, we are defining position for spike protein)

startnt <- 21563
endnt <- 25384


#select specified rows from the letter data frame
msa_sel1 <- msa[startnt:endnt,] 

#Extract position column
pos <- msa_sel1[1]

msa_sel <- msa_sel1[,-1]  #remove the "POS" column
#Report on the number of mutations in each sample:
x <- ncol(msa_sel)

#create an empty data frame to store the new data
#df_new1 <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))

colnames(df_new1) <- colnames(msa_sel)

#define i value
i=2

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_new1[,i] <- ifelse(msa_sel[,1]== msa_sel[,i],"-", paste0(msa_sel[,1],">",msa_sel[,i]))
}


#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
df_new2 <- select(df_new1,-1, -ncol(df_new1))

#Add position column
final_df <- cbind(pos, df_new2)

#Check whether output is correct
head(final_df)

#Write into a file
write.table(final_df, file="variant_info.txt", row.names = F, sep ="\t")

dim(final_df)