library(ape)
library(pegas)

#load data (FASTA files)
seq<- read.dna("Popgendata_Begoro_copy.fa", format="fasta") 

#check loaded data
seq

#convert the DNA sequence data into haplotypes
AMA1_Haplotypes <- haplotype(seq)

haps2 <- subset(AMA1_Haplotypes, minfreq = 3)
(network <- haploNet(haps2))

plot(network, size = attr(network, "freq"),show.mutation=1,labels=T)

#print haplotypes
barplot(AMA1_Haplotypes, las=2)

write.table(AMA1_Haplotypes, file="AMA1_Haplotypes.txt")

#print all haplotypes
summary(AMA1_Haplotypes)

#Create haplotype network
AMA1_Net <- haploNet(AMA1_Haplotypes)
AMA1_Net

#Plot haplotype network
#Note: Use fast = TRUE, just to explore the data  
plot(AMA1_Net, size = attr(AMA1_Net, "freq"), fast = TRUE)

# Note: use fast = FALSE, because the structure of the resulting plot is less messy
plot(AMA1_Net, size = attr(AMA1_Net, "freq"), fast = FALSE)

#Next to know which samples belong to which haplotype
ind.hap<- with(stack(setNames(attr(AMA1_Haplotypes, "index"), rownames(AMA1_Haplotypes))), table(hap=ind, pop=rownames(seq)[values]))

#Print haplotype list
ind.hap



barplot(ind.hap, las=2)

#convert into dataframe
mydata <- as.data.frame(ind.hap)
mydata 
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

#Create haplotype network
plot(AMA1_Net, size=attr(AMA1_Net, "freq"), scale.ratio = 2, cex = 0.5, pie=new.hap)


plot(AMA1_Net, size=attr(AMA1_Net, "freq"), fast = TRUE, scale.ratio = 2, cex = 0.8, pie=new.hap)

legend("topright", c("Begro", "CapeCoast", "red"), text.col=2:5)