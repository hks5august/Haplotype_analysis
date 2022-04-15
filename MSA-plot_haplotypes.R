#load libraries
library(gplots)

#load and prepare data

#for full gene
#
msa_sel1 <- read.table("42_key_haplotypes_variants2.txt", header=TRUE, sep="\t") 

# for domain1
msa_sel1 <- read.table("27_haplotypes_variants.txt", header=TRUE, sep="\t") 
dim(msa_sel1)

head(msa_sel1)

pos <- msa_sel1[1]
ref<- msa_sel1[2]

msa_sel <- msa_sel1[,-1]  #remove the "POS" column
#Report on the number of mutations in each sample:
x <- ncol(msa_sel)

msa_sel 

startnt <- 1
endnt <- 81
#create an empty data frame to store the new data
df_new1 <- data.frame(matrix(0, ncol = x+1, nrow = (endnt-startnt)+1))
colnames(df_new1) <- colnames(msa_sel)

#define i value
i=2

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_new1[,i] <- ifelse(msa_sel[,1]== msa_sel[,i],"-", paste0(msa_sel[,i]))
}


df_new1
#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
df_new2 <- select(df_new1,-1, -ncol(df_new1))

#Add position column
final_df <- cbind(pos,ref, df_new2)
final_df 
final_df1 <- cbind(ref, df_new2)
final_df1



#select only data where any sample has a variant to create a report
#new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9] & msa1sel[2] == msa1sel[10]),] 
#new_dfT <- t(new_df)
new_dfT <- t(final_df1)



#transform the data for ggplot
library(reshape2)
melted_mat <- melt(new_dfT)
head(melted_mat)

melted_mat1 <- melt(msa_sel1T)

#draw the plot
( plot1 <- ggplot(melted_mat) + #add data and fill cells
    geom_tile(data = melted_mat, aes(x=Var2, y=Var1, fill=value), colour = "black") + #add black border around cells
    geom_text(data = melted_mat, aes(x=Var2, y=Var1,label = value), size=2, family = "Avenir Next") + 
    scale_fill_manual(values = c("-" = "lightgray", "A" = "lightgreen", "C" = "pink", "G" = "lightblue", "T" = "yellow", "V"="red", "*"="white")) +
    coord_equal() +
    xlab("Nucleotide Position") +
    ylab("Haplotypes") +
    labs(fill = "NT") +
    labs(title="") +
    scale_x_continuous(breaks = seq(startnt, endnt, by = 2), expand = c(0, 0)) +
    theme(axis.text = element_text(size=7, family = "Avenir Next"),
          legend.title = element_text(color = "gray", size = 8),
          legend.text = element_text(color = "gray"),
          axis.title = element_text(size=7, vjust = 2, face="italic"),
          axis.text.x = element_text(angle = 90, vjust = 0.7, size = 7),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    ))

#make the first plot


png('AMA1_MSA.png', units="in", width=6, height=5, res=250)
ggplot(melted_mat) + #add data and fill cells
      geom_tile(data = melted_mat, aes(x=Var2, y=Var1, fill=value), colour = "black") + #add black border around cells
      geom_text(data = melted_mat, aes(x=Var2, y=Var1,label = value), size=1, family = "Avenir Next") + 
      scale_fill_manual(values = c("-" = "lightgray", "A" = "lightgreen", "C" = "pink", "G" = "lightblue", "T" = "yellow", "V"="red", "*"="white")) +
      coord_equal() +
      xlab("Nucleotide Position") +
      ylab("Haplotypes") +
      labs(fill = "NT") +
      labs(title="") +
      scale_x_continuous(breaks = seq(startnt, endnt, by = 2), expand = c(0, 0)) +
      theme(axis.text = element_text(size=3, family = "Avenir Next"),
            legend.title = element_text(color = "gray", size = 7),
            legend.text = element_text(color = "gray"),
            axis.title = element_text(size=4, vjust = 2, face="italic"),
            axis.text.x = element_text(angle = 90, vjust = 0.7, size = 4),
            axis.text.y = element_text(size = 4),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5)
      )  
 dev.off()

plot(plot1)













###############################

#load and prepare data
#msa1 <- read.csv("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/SARS_CoV_2_1.csv", header=TRUE)

msa1<- read.table("27_haplotypes_variants.txt", header=TRUE, sep="\t") 

#prepare a numeric matrix for msa1 data frame
msa1m <- data.matrix(msa1) 

#remove the "POS" column
msa1 <- msa1[,-1] 

#select from the letter data frame
msa1sel <- msa1

#select from the numeric matrix
msa1msel <- msa1m

#for numeric matrix, make sure to add the position to row names
row.names(msa1msel) <- msa1msel[,1] 

#now, remove the POS column from the input data
msa1msel <- msa1msel[,-1]

#select only data where any sample has a variant to create a report
#new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9] & msa1sel[2] == msa1sel[10] & msa1sel[2] == msa1sel[11] & msa1sel[2] == msa1sel[12] & msa1sel[2] == msa1sel[13] & msa1sel[2] == msa1sel[14]),] 

final_df <- cbind(pos,ref, df_new2)
final_df 
final_df1 <- cbind(ref, df_new2)
final_df1



#select only data where any sample has a variant to create a report
#new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9] & msa1sel[2] == msa1sel[10]),] 
#new_dfT <- t(new_df)
#new_dfT <- t(final_df1)
new_matrix <- data.matrix(final_df1)


new_dfT <- t(final_df1)

#select only data where any sample has a variant in the matrix
#new_matrix <- msa1msel[!(msa1msel[,2] == msa1msel[,3] & msa1msel[,2] == msa1msel[,4] & msa1msel[,2] == msa1msel[,5] & msa1msel[,2] == msa1msel[,6] & msa1msel[,2] == msa1msel[,7] & msa1msel[,2] == msa1msel[,8] & msa1msel[,2] == msa1msel[,9] & msa1msel[,2] == msa1msel[,10] & msa1msel[,2] == msa1msel[,11] & msa1msel[,2] == msa1msel[,12] & msa1msel[,2] == msa1msel[,13] & msa1msel[,2] == msa1msel[,14]),]

#transpose
new_matrixT <- t(new_matrix)

#set colors and margins for plot
if("-" %in% new_dfT){ #first, check if we have 4 or 5 characters (I am assuming that TCGA will be present, so only looking for "-")
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
          #main = cbind("MSA of Haplotypes"), #heat map title
          #plot title
          main = cbind("Number of Variants Found: ",ncol(new_matrixT)), #heat map title
          
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
