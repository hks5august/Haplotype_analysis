library(ape)
library(pegas)

#load data (FASTA files)
seq<- read.dna("Popgendata_BegCape_combined.fa", format="fasta") 

#check loaded data
seq

#convert the DNA sequence data into haplotypes
AMA1_Haplotypes <- haplotype(seq)

#print all haplotypes
summary(AMA1_Haplotypes)

#Change rownames of haplotypes
row.names(AMA1_Haplotypes) <- 1 : 194 #look at no accordingly change it
#rownames(AMA1_Haplotypes) <- paste0('hap', rownames(AMA1_Haplotypes))

#print all haplotypes
summary(AMA1_Haplotypes)

#Extract mutations with position for all the haplotypes
all_haplotypes_diff <- as.data.frame(diffHaplo(AMA1_Haplotypes,1:194)) # no. can be changed after checking how many total haplotypes are there

#Write in a file
write.table(all_haplotypes_diff, file= "all_haplotypes_changes.txt")


#Create haplotype network
AMA1_Net <- haploNet(AMA1_Haplotypes)
AMA1_Net2 <- mjn(AMA1_Haplotypes)

#Plot haplotype network
#Note: Use fast = TRUE, just to explore the data  
plot(AMA1_Net, size = attr(AMA1_Net, "freq"), fast = TRUE, labels=F, show.mutation=2)

# Note: use fast = FALSE, because the structure of the resulting plot is less messy
plot(AMA1_Net, size = attr(AMA1_Net, "freq"), fast = FALSE)

#Next to know which samples belong to which haplotype
ind.hap<- with(stack(setNames(attr(AMA1_Haplotypes, "index"), rownames(AMA1_Haplotypes))), table(hap=ind, pop=rownames(seq)[values]))

#Print haplotype list
ind.hap

#Barplot
barplot(ind.hap, las=2)

#Write into a file
write.table(ind.hap, file="individual_sample_haplotype.txt")

#Next to convert based on region/loacation
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

#Draw a barplot
barplot(new.hap,  main= "No of.Haplotypes of AMA1 in two regions", cex.axis = 0.8, cex.lab=0.9, cex.names=0.8,  ylab="Number of Haplotypes", xlab="Location",  space=0.2)

#Write into a file
write.table(new.hap, file="group_haplotype.txt")

#Create haplotype network


#Create haplotype network (quickly), but it woukd be messy without labels
plot(AMA1_Net, size = attr(AMA1_Net, "freq"), fast = T, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

#Create haplotype network (quickly), but it woukd be messy with labels of haplotypes
plot(AMA1_Net, size = attr(AMA1_Net, "freq"), fast = T, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)




######### Select frequency criteria to plot

haps1 <- subset(AMA1_Haplotypes, minfreq = 1)
haps2 <- subset(AMA1_Haplotypes, minfreq = 2)
haps3 <- subset(AMA1_Haplotypes, minfreq = 3)


############################# Freuqncy 1 ###########
haps1 <- subset(AMA1_Haplotypes, minfreq = 1)


#Next to know which samples belong to which haplotype
# Note: Change the "haps" value according to the frequncy of haplotypes you want choose for network


haps_sort1 <- sort(haps1)


ind.hap<- with(stack(setNames(attr(haps1, "index"), rownames(haps1))), table(hap=ind, pop=rownames(seq)[values]))

#convert it into matrix
indiv_haplotypes <- as.matrix(ind.hap)
#barplot(indiv_haplotypes, las=2)

#Write in file
#write.table(indiv_haplotypes, file ="indiv_hapltype.txt")

#create group or location-wise table
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

#write.table(new.hap, file="new_haplotypes.txt")

#Create haplotype network 
haps_sort <- sort(haps1)
network <- haploNet(haps1, getProb = TRUE)



network11 <- haploNet(haps_sort)

#Plot without labels of samples, to quickly (fast=T), probably messy network
plot(network, size = attr(network, "freq"), fast = T, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 1) for Begro and Capecoast", pie=new.hap, lwd = 0.1, scale=1)

#Plot without labels of samples
plot(network1, size = attr(network1, "freq"), fast = T, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

#Plot proper network (fast=F) without labels
plot(network1, size = attr(network1, "freq"), fast = F, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

#Plot proper network (fast=F) with labels
plot(network1, size = attr(network1, "freq"), fast = F, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)



#Plot proper network (fast=F) with labels
plot(network22, size = attr(network22, "freq"), fast = F, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)


#Add legends
legend(180,40, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)




########################### Freuqncy 2 ################
haps2 <- subset(AMA1_Haplotypes, minfreq = 2)


#Next to know which samples belong to which haplotype
# Note: Change the "haps" value according to the frequncy of haplotypes you want choose for network


#haps_sort2 <- sort(haps2)
#haps = haps_sort2

ind.hap<- with(stack(setNames(attr(haps2, "index"), rownames(haps2))), table(hap=ind, pop=rownames(seq)[values]))

#convert it into matrix
indiv_haplotypes <- as.matrix(ind.hap)
#barplot(indiv_haplotypes, las=2)

#Write in file
#write.table(indiv_haplotypes, file ="indiv_hapltype.txt")

#create group or location-wise table
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

#write.table(new.hap, file="new_haplotypes.txt")

#Create haplotype network 
haps_sort2 <- sort(haps2)
network <- haploNet(haps, getProb = TRUE)


network1 <- haploNet(haps2)


network22 <- haploNet(haps_sort2)

#Plot without labels of samples, to quickly (fast=T), probably messy network
plot(network1, size = attr(network1, "freq"), fast = T, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.1, scale=1)

#Plot without labels of samples
plot(network1, size = attr(network1, "freq"), fast = T, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

#Plot proper network (fast=F) without labels
plot(network1, size = attr(network1, "freq"), fast = F, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

png(
  "test.png",
  width     = 10,
  height    = 10,
  units     = "in",
  res       = 150,
  pointsize = 1
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 10,
  cex.main = 20,
  cex.lab  = 10
)
plot(network1, size = attr(network1, "freq"), fast = F, show.mutation=2, labels=T, cex = 10, main = "", pie=new.hap, lwd = 0.5, scale=1)
 legend(-50,10, colnames(new.hap), col=rainbow(ncol(new.hap)), cex=10, pch=20)
dev.off()

#Plot proper network (fast=F) with labels
plot(network1, size = attr(network1, "freq"), fast = F, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 2) for Begoro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)
#Add legends
#Add legends
legend(-130,10, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)

dev.off()

legend(180,40, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)


#Plot proper network (fast=F) with labels
plot(network22, size = attr(network22, "freq"), fast = F, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 2) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)





############ for Frequency 3



haps3 <- subset(AMA1_Haplotypes, minfreq = 3)


#Next to know which samples belong to which haplotype
# Note: Change the "haps" value according to the frequncy of haplotypes you want choose for network
haps3

haps_sort3 <- sort(haps3)


ind.hap<- with(stack(setNames(attr(haps3, "index"), rownames(haps3))), table(hap=ind, pop=rownames(seq)[values]))

#convert it into matrix
indiv_haplotypes <- as.matrix(ind.hap)
#barplot(indiv_haplotypes, las=2)

#Write in file
#write.table(indiv_haplotypes, file ="indiv_hapltype.txt")

#create group or location-wise table
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

#write.table(new.hap, file="new_haplotypes.txt")

#Create haplotype network 
haps_sort <- sort(haps1)
network3 <- haploNet(haps3)
network4 <- mjn(haps3)

network33 <- haploNet(haps_sort3)


#Plot without labels of samples, to quickly (fast=T), probably messy network
plot(network3, size = attr(network3, "freq"), fast = T, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.1, scale=1)
plot(network4,  fast = T, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.1, scale=1)


#Plot without labels of samples
plot(network3, size = attr(network3, "freq"), fast = T, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

#Plot proper network (fast=F) without labels
plot(network3, size = attr(network3, "freq"), fast = F, show.mutation=2, labels=F, cex = 0.5, main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)

#Plot proper network (fast=F) with labels
plot(network3, size = attr(network3, "freq"), fast = F, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)



#Plot proper network (fast=F) with labels
plot(network33, size = attr(network33, "freq"), fast = F, show.mutation=2, labels=T, cex = 0.5, main = "Haplotype network (freq >= 3) for Begro and Capecoast", pie=new.hap, lwd = 0.5, scale=1)


#Add legends
legend(180,40, colnames(new.hap), col=rainbow(ncol(new.hap)), pch=20)

