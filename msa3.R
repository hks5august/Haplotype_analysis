################# MSA for variants ############
#load libraries
library(gplots)

#load and prepare data
#msa1 <- read.csv("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/SARS_CoV_2_1.csv", header=TRUE)
#msa1 <- read.table("input_data.txt", sep="\t",  header=TRUE)

msa1 <- read.table("42_key_haplotypes_variants2.txt", header=TRUE, sep="\t") 

dim(msa1)
#prepare a numeric matrix for msa1 data frame
msa1m <- data.matrix(msa1) 

#remove the "POS" column
msa1 <- msa1[,-1] 

startnt <- 1
endnt <- 81

#select from the letter data frame
msa1sel <- msa1[startnt:endnt,]
msa1sel
#select from the numeric matrix
msa1msel <- msa1m[startnt:endnt,]
msa1msel
#for numeric matrix, make sure to add the position to row names
row.names(msa1msel) <- msa1msel[,1] 
msa1msel
#now, remove the POS column from the input data
msa1msel <- msa1msel[,-1]



startnt <- 1
endnt <- 81
#create an empty data frame to store the new data
new_df <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))
colnames(new_df) <- colnames(msa1sel)
new_df
#define i value
i=2
x <- ncol(msa1msel)
#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  new_df[,i] <- ifelse(msa1sel[,1]== msa1sel[,i],"-", paste0(msa1sel[,i]))
}


new_df


#select only data where any sample has a variant to create a report
#new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9]),] 
new_dfT <- t(new_df)

#create an empty data frame to store the new data
new_matrix <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))
colnames(new_matrix) <- colnames(msa1sel)

#define i value
i=2

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  #new_matrix[,i] <- ifelse(msa1msel[,1]== msa1msel[,i],"-", paste0(msa1msel[,i]))
  new_matrix[,i] <- ifelse(msa1msel[,1]== msa1msel[,i],"5", paste0(msa1msel[,i]))
  
  }


new_matrix

#select only data where any sample has a variant in the matrix
#new_matrix <- msa1msel[!(msa1msel[,2] == msa1msel[,3] & msa1msel[,2] == msa1msel[,4] & msa1msel[,2] == msa1msel[,5] & msa1msel[,2] == msa1msel[,6] & msa1msel[,2] == msa1msel[,7] & msa1msel[,2] == msa1msel[,8] & msa1msel[,2] == msa1msel[,9]),]

#transpose
new_matrixT <- t(new_matrix)

#set colors and margins for plot
if("-" %in% new_matrixT){ #first, check if we have 4 or 5 characters (I am assuming that TCGA will be present, so only looking for "-")
  TCGAcolors <- c("white", "lightgreen", "pink", "lightblue", "yellow") 
  names(TCGAcolors) = c(1,2,3,4,5) #labels = c("-","A","T","C","G")
} else {TCGAcolors <- c("lightgreen", "pink", "lightblue", "yellow") }




#draw heat map for only variants
heatmap.2(new_matrixT, #data source
          
          #main settings
          cexRow = 0.7, #row name font size
          col = TCGAcolors, #set colors
          dendrogram = "none", #remove dendrogram
          Rowv = FALSE, #no reordering for rows
          Colv = FALSE, #no reordering for columns
          density.info="none", #remove density info
          trace="none", #remove row and column lines
          
          offsetRow=0.1, #change position of the row names
          offsetCol=0.1, #change position of the column names
          
          #add gray borders between cells
          sepwidth=c(0.05,0.05), #sets separation width and height
          sepcolor="gray", #color for border
          colsep=1:ncol(new_matrixT), #add separation for number of columns in source data
          rowsep=1:nrow(new_matrixT), #add separation for number of rows in source data
          
          #plot title
          main = cbind("No. of Variants:",ncol(new_matrixT)), #heat map title
          cex.main=0.8,
          #plot margins
          margins = c(5,10), #set margins
          lwid=c(0.2,4), 
          lhei=c(0.9,3),
          
          #adding letters inside the heatmap
          notecex=1.0, #size of font inside each cell
          cellnote = new_dfT, #data to use in cells
          notecol="black", #font color for cells
          
          #legend
          key = FALSE
)