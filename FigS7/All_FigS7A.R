#read in the file
data <- read.csv("MEMLg_Supp.csv")
head(data)

#Convert to a matrix and turn metabolites into rownames
dim(data)
data.mat <- as.matrix(data[,5:79])
row.names(data.mat) <- data$BIOCHEMICAL

#Remove rows where standard deviation = 0 
data.mat <- subset(data.mat, apply(data.mat,1, sd)!=0)

#Load strain info 
Sinfo <- read.csv("Sinfo_all.csv",stringsAsFactors = FALSE)

#heatmap of everything 
library(gplots)

##Create a PCA plot and get contributions data for specific genes##
H_PCA <- prcomp(t(data.mat), center = TRUE, scale = TRUE)
H_PCA

#PCA plots by different conditions
MEMtype <- Sinfo[['MEM']]
H_PCA_MEM <- ggbiplot::ggbiplot(pcobj = H_PCA, var.axes = FALSE, groups = MEMtype, circle.prob = TRUE, ellipse = TRUE, circle = TRUE, obs.scale = 1, var.scale = 1, size = 2, choices = c(1,2)) +
  geom_point(aes(colour=MEMtype), size = 2)+theme_grey(base_size = 22)
H_PCA_MEM

#Figure S7A
library(viridis)
H_PCA_MEM + theme_bw() + 
  theme(text = element_text(size = 20)) +
  scale_color_viridis(discrete = TRUE, 
                      name = "Condition", 
                      labels=c("MEM Control","MEM Supernatants","sMEM Control","sMEM Supernatants"))

