#*********************
#LAB COURSE 2018
#Diversity indices
#Community ordination
#C. Karakoc
#*********************

#*****************************
#Install packages before class 
#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("reshape")
#*****************************

library(vegan)
library(ggplot2)
library(reshape)

# This is the place where R will look for files
# and where it will save files

# Set working directory (change with your working directory)
# run from the folder this script is in

# Haribo data 
haribo <- read.table("HariboMatrix.csv", sep=";", dec=",", header=T)

# Inspect data 
str(haribo)
head(haribo)

# Mesure diversity 
richness <- specnumber(haribo[,2:25]) #Richness

shannon <- diversity(haribo[,2:25], index="shannon") #Shannon
hill1   <- exp(shannon) #Exponential Shannon 

simpson <- diversity(haribo[,2:25], index="simpson") #Simpson's index
hill2 <- 1/simpson #Inverse Simpson

#Create a data set with all
haribo.data <- as.data.frame(cbind(richness,hill1,hill2))
haribo.data$sampleSize <- haribo$SampleSize
haribo.data$community <- haribo$Site

#Reshape the data for plotting 
haribo.data.long <- melt(haribo.data, id=c("sampleSize", "community"))

#Plot 
haribo.plot <- ggplot(haribo.data.long,aes(x=sampleSize, y=value, color=community))+
  geom_point()+
  geom_line()+
  facet_wrap(.~variable, scales="free")+
  theme_bw()
  
                        
#***********
#TRFLP Data
#***********

# Read data using read.delim
data <- read.delim("ClusBinMatrix_B_3_Praktikum_18S.txt") 

# Inspect data
str(data)
head(data) # shows first 6 rows

#Change the fragment names here
#This line produces rounded fragment
#Remove X from the column names 
colnames(data)<-gsub("X","",(colnames(data)))
dput(as.character((round((as.numeric(colnames(data[,-1]))), 1))))

colnames(data) <- c("samples", "50.3", "51.9", "53.7", "56.3", "60.6", "62.2", "64.3", "68.8", 
                    "70.4", "71.4", "78", "80", "85.8", "106.3", "110.3", "112.8", 
                    "137.8", "150.3", "157.1", "158.7", "163.9", "168.6", "171", 
                    "172.2", "173.2", "180.6", "182.5", "183.5", "184.8", "188.6", 
                    "191.2", "205.4", "207", "219.1", "220.5", "221.8", "233.7", 
                    "240.3", "255.1", "256.5", "260.9", "268.6", "269.8", "272.1", 
                    "301.7", "312.7", "320.4", "341", "369.9", "371.1", "373.6", 
                    "382.5", "384", "385.5", "386.9", "388", "393.1", "399.6", "401", 
                    "402.3", "404.5", "406.9", "408", "409.2", "410.8", "412.3", 
                    "414.5", "416", "417.5", "419.1", "420.3", "422.4", "424.4", 
                    "426.9", "430.7", "432.1", "433.5", "439.6", "446.9", "450", 
                    "471.2", "475.6", "490.5")

#*******************
#Diversity indices
#*******************

# Mesure diversity 
richness <- specnumber(data[,-1]) #species richness 

shannon <- diversity(data[,-1], index="shannon") #Shannon entropy
hill1   <- exp(shannon) #Exponential Shannon

simpson <- diversity(data[,-1], index="simpson") #Simpson's Index
hill2 <- 1/simpson #inverse Simpson

#Lets combine all and plot
diversity.data <- as.data.frame(cbind(richness, hill1, hill2))
diversity.data$samples <- data[,1]
diversity.data.long <- melt (diversity.data, id="samples")

#Plot indices
indices.plot <- ggplot(diversity.data.long, aes(x=samples, y=value, color=variable))+ 
  geom_point(stat="identity")+
  ylab("Index value")+
  xlab("Samples")+
  facet_wrap(.~variable, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position = "none")

#*******************
#Relative Abundances 
#*******************

#Longitidunal form
species <- melt(data, id="samples")
head(species)

species.plot <- ggplot(species, aes(x=samples, y=value, fill=variable))+ 
  geom_bar(stat="identity")+
  ylab("Relative abundances")+
  xlab("Samples")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

#****************************************
#Non-metric multidimensional scaling NMDS
#****************************************

mdsB     <- metaMDS(data[,-1], distance="bray")
stressplot(mdsB)

#NMDS Data
data.scores <- as.data.frame(scores(mdsB, display = "sites")) #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.mds <- cbind(samples=data[,1], data.scores) #Combine it with sample names
head(data.mds)

species.scores <- as.data.frame(scores(mdsB, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
head(species.scores)  


#Plot
plot_mds <- ggplot(data=data.mds,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(colour=samples),size=5)+
  #geom_text(data=species.scores, aes(x=NMDS1, y=NMDS2, label=rownames(species.scores)), 
  #                                   size=3, check_overlap = TRUE)+
  theme_bw()


#Bonus: Hierarchical Clustering
rownames(data)<- data$samples 
clusters <- hclust(vegdist(data[,-1]))
plot(clusters)

#STOP HERE!

#**********************************
#Why we need multivariate analysis
#**********************************

plot(0:10,0:10,type="n",axes=F,xlab="Abundance of Species 1",ylab="") 
axis(1)
points(1,0); text(1.2,0.5,labels="community A")
points(6,0); text(6.2,0.5,labels="community B")


plot(0:10,0:10,type="n",xlab="Abundance of Species 1",
     ylab="Abundance of Species 2")
points(1,5); text(1,4.5,labels="community A")
points(6,3); text(6,3.5,labels="community B")


#install.packages("scatterplot3d")
library(scatterplot3d)
d=scatterplot3d(0:10,0:10,0:10,type="n",xlab="Abundance of Species 1",
                ylab="Abundance of Species 2",zlab="Abundance of Species 3"); d
d$points3d(1,5,0); text(d$xyz.convert(5,5,0.5),labels="community A")
d$points3d(6,3,3); text(d$xyz.convert(3,3,3.5),labels="community B")








