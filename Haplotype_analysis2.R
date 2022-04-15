library(ape)
library(pegas)

#load data (FASTA files)
seq<- read.FASTA("Popgendata_Begoro_with_ref_transl.fa", type="AA") 

seq<- read.dna("Popgendata_BegCape_combined.fa", format="fasta") 



#check loaded data
seq


#convert the DNA sequence data into haplotypes
AMA1_Haplotypes <- haplotype(seq)

#print haplotypes
AMA1_Haplotypes

all_haplotypes_diff <- as.data.frame(diffHaplo(AMA1_Haplotypes,1:194))

write.table(all_haplotypes_diff, file= "all_haplotypes_changes.txt")


haps1 <- subset(AMA1_Haplotypes, minfreq = 7)
haps1
######### Select frequency criteria to plot

haps1 <- subset(AMA1_Haplotypes, minfreq = 1)
haps2 <- subset(AMA1_Haplotypes, minfreq = 2)
haps3 <- subset(AMA1_Haplotypes, minfreq = 3)
(network <- haploNet(haps1, getProb = TRUE))

haps1

#Next to know which samples belong to which haplotype
ind.hap<- with(stack(setNames(attr(haps1, "index"), rownames(haps1))), table(hap=ind, pop=rownames(seq)[values]))



indiv_haplotypes <- as.matrix(ind.hap)
#barplot(indiv_haplotypes, las=2)

#write.table(indiv_haplotypes, file ="indiv_hapltype.txt")


#convert into dataframe
mydata <- as.data.frame(ind.hap)
hap_list<- mydata[mydata$Freq == 1,]
hap_list

# let create haplotypes based on locations, here, string split the names by underscores.
locations <- strsplit(as.character(hap_list$pop), "_")
locations

# Next extract first item in each list
locations1 <- sapply(locations, "[[", 1)

head(locations1)

#Now make a table with our new locations list and the corresponding haplotypes
new.hap <- table(hap_list$hap, locations1)
new.hap

write.table(new.hap, file="new_haplotypes.txt")

barplot(new.hap)

plot(network, size = attr(network, "freq"),show.mutation=2, labels=TRUE, cex = 0.5, main = "Haplotype network (freq >= 1) for Begro and Capecoast", pie=new.hap, lwd = 0.8, lty=0.4)
