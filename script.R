

#load libraries
library(gplots)
library(dplyr) #required for modification of data
library(pryr) #required to save basic rplots as object
library(ggplot2)
library(gt)
library(ggpubr)
library(gridExtra)


setwd("../../AMA1-all_domains/Year_wise/")
msa_sel1 <- read.table("Domain1_yearwise_mutation_input.txt", header=TRUE, sep="\t") 
msa_sel1 <- read.table("Domain1_yearwise_mutation_input.txt", header=TRUE, sep="\t") 


dim(msa_sel1)

head(msa_sel1)

pos <- msa_sel1[1]
ref<- msa_sel1[2]

msa_sel <- msa_sel1  #remove the "POS" column
#Report on the number of mutations in each sample:
x <- ncol(msa_sel)
x

y <- nrow(msa_sel)
y

msa_sel 


#create an empty data frame to store the new data
df_new1 <- data.frame(matrix(0, ncol = x+1, nrow = y))
colnames(df_new1) <- colnames(msa_sel)

head(msa_sel)
#define i value
i=3

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_new1[,i] <- ifelse(msa_sel[,2]== msa_sel[,i],"-", paste0(msa_sel[,i]))
}


i=3

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 3:x) {
  
  df_new1[,i] <- ifelse(msa_sel[,2]== msa_sel[,i],"-", paste0(msa_sel[,2],">",msa_sel[,i]))
}



df_new1
#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
df_new2 <- select(df_new1,-1, -ncol(df_new1))

#Add position column
final_df <- cbind(df_new2)
final_df 

row.names(final_df) <- msa_sel1 [,1] 

head(final_df)

write.table(final_df, file= "AMA1_Domain1_Haplotypes_variant_summary.txt", row.names = T, sep ="\t")

final_df %>% gt()


#summary table

main.title <- paste("Summary of mutations variants in Haplotypes of AMA1-Domain1")
mutation_summary <- ggtexttable(final_df, rows = NULL, theme = ttheme("light", padding = unit(c(11, 2.5), "mm"), tbody.style = tbody_style(fill = "white", size = 9), colnames.style = colnames_style(fill = "white", size = 10, color = "Black", rot=90)))
mutation_summary1 <- mutation_summary %>% tab_add_title(text = main.title, face = "bold",  padding = unit(5, "line"))

mutation_summary1


Protein <- 'AMA1_haplotypes_Variants.txt'


report_name <- toString(sub(".txt", ".pdf", Protein))


grid <- grid.arrange(mutation_summary1, nrow=1, ncol=1)
ggsave(grid, file=report_name, height = 30 , width = 30, device = pdf,  dpi = 300)




#Generate MSA variant plot
#Transpose
new_dfT <- t(final_df)
new_dfT
write.table(new_dfT, file= "AMA1_Haplotypes_variant_info.txt", row.names = T, sep ="\t")


library(reshape2)
melted_mat <- melt(new_dfT)

head(melted_mat)


#draw MSA the plot
plot1 <- ggplot(melted_mat) + #add data and fill cells
  geom_tile(data = melted_mat, aes(x=factor(Var1), y=factor(Var2), fill=value), colour = "black") + #add black border around cells
  geom_text(data = melted_mat, aes(x=factor(Var1), y=factor(Var2),label = value), size=2, family = "Helvetica") + 
  scale_fill_manual(values = c("-" = "lightgray", "A" = "lightgreen", "C" = "pink", "G" = "lightblue", "T" = "yellow", "V"="red", "*"="white")) +
  #coord_equal() +
  ylab("Nucleotide Position") +
  xlab("Haplotypes") +
  labs(fill = "NT") +
  labs(title= paste("Number of variants in the Haplotypes for AMA1:" , ncol(new_dfT)  )) +
  theme(axis.text = element_text(size=6, family = "Helvetica"),
        legend.title = element_text(color = "black", size = 6),
        legend.text = element_text(color = "black"),
        axis.title = element_text(size=8, vjust = 2, face="italic"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

#make the first plot
plot(plot1)


Protein <- 'AMA1_haplotypes_MSA.txt'


report_name <- toString(sub(".txt", ".pdf", Protein))


grid <- grid.arrange(plot1, nrow=1, ncol=1)
ggsave(grid, file=report_name, height = 30 , width = 30, device = pdf,  dpi = 300)





########
head(new_dfT)



#Report on the number of mutations in each sample:
x <- ncol(msa_sel)
x

y <- nrow(msa_sel)
y




#create an empty data frame to store the new data
df_mut <- data.frame(matrix(0, ncol = x+1, nrow = y))
colnames(df_mut) <- colnames(msa_sel)

#define i value

i=3

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_mut[,i] <- ifelse(msa_sel[,2]== msa_sel[,i],0,1)
}



#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
df_mut <- select(df_mut, -1, -ncol(df_mut))
#df_mut%>% gt()
head(df_mut)
#Compute sum of mutations per samples
summ1 <- as.data.frame(colSums(df_mut))

#Samples
samples<- as.data.frame(row.names(summ1))

#Combine samples and sum
summ2 <- cbind(samples, summ1)

#Add column names
colnames(summ2) <- c("Sample", "Mutations_Number")

head(summ2)


#make a barplot
bb <- ggplot(summ2, aes(x = Sample, y = Mutations_Number)) +geom_col(aes(width = 0.8))
bb+ theme(legend.key.size = unit(0.6, "cm"), legend.text = element_text( color="Black", size=9), axis.text.x = element_text( color="Black", size=10, angle=90), axis.text.y = element_text( color="Black", size=10, angle=90))


bb <- ggplot(summ2, aes(x=Sample, y=Mutations_Number)) + geom_col(stat="identity", fill="steelblue", width = 0.5) +   
  labs(title= paste("No. of mutations per Haplotype" )) +
  theme(axis.text = element_text(size=10),
        legend.title = element_text(color = "black", size = 8),
        legend.text = element_text(color = "black"),
        axis.title = element_text(size=8, vjust = 2, face="italic"),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, size = 8),
        axis.text.y = element_text(color = "black",size = 8, angle = 90,vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

bb                                                                             



Protein <- 'AMA1_haplotypes_mutations.txt'


report_name <- toString(sub(".txt", ".pdf", Protein))


grid <- grid.arrange(bb, nrow=1, ncol=1)
ggsave(grid, file=report_name, height = 30 , width = 30, device = pdf,  dpi = 300)









