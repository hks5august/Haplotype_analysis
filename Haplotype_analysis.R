library(ape)
library(pegas)

#load data (FASTA files)
seq<- read.dna("Popgendata_BegCape_combined.fa", format="fasta") 

#check loaded data
seq

#convert the DNA sequence data into haplotypes
AMA1_Haplotypes <- haplotype(seq)

#print haplotypes
AMA1_Haplotypes

all_haplotypes_diff <- as.data.frame(diffHaplo(AMA1_Haplotypes,1:194))

write.table(all_haplotypes_diff, file= "all_haplotypes_changes.txt")

###
diffHaplo(AMA1_Haplotypes, a = "I", b = "V")

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


ind.hap1<-with(stack(setNames(attr(AMA1_Haplotypes, "index"), rownames(AMA1_Haplotypes))), table(hap=ind, pop=rownames(seq)[values]))

barplot(ind.hap1, las=2)

write.table(ind.hap1, file="individual_sample_haplotype.txt")



#convert into dataframe
mydata <- as.data.frame(ind.hap1)
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

write.table(new.hap, file="group_haplotype.txt")


#Create haplotype network
plot(AMA1_Net, size=attr(AMA1_Net, "freq"), scale.ratio = 2, cex = 0.5, pie=new.hap)

#Create haplotype network (quickly), but it woukd be messy
plot(AMA1_Net, size=attr(AMA1_Net, "freq"), fast = TRUE, scale.ratio = 2, cex = 0.8, pie=ind.hap)

legend("topright", c("Begro", "CapeCoast", "red"), text.col=2:5)

plot(AMA1_Net, size=attr(AMA1_Net, "freq"),  scale.ratio = 2, cex = 0.8, pie=ind.hap)
legend("topright", c("Begro", "CapeCoast", "red"), text.col=2:5)


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




plot(network, size = attr(network, "freq"),show.mutation=2, labels=TRUE, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.8, lty=0.4)

plot(network, size = attr(network, "freq"),show.mutation=2,labels=TRUE, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", threshold = c(1, 14), pie=new.hap, lwd = 0.5, lty=0.5)


legend(-130,10, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)






network1 <- haploNet(haps3, getProb = TRUE)
plot(network1, size = attr(network1, "freq"),show.mutation=2, labels=TRUE, cex = 0.5, threshold = c(1, 14), main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.8, lty=0.4)

legend(45,-50, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)




legend = c(-25, 30))


legend(50,50, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)

legend("bottomleft", colnames(new.hap), col=rainbow(ncol(new.hap)), pch=10)

legend(50,50, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)



legend("bottomleft", c("Begro", "CapeCoast"), text.col=2:6)

haps1

as.data.frame(diffHaplo(haps1,1:65))



diffHaplo(AMA1_Haplotypes, a = "I", b = "V")