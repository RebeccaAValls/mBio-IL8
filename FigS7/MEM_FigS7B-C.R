#read in the file
data <- read.csv("MEMmin.csv")
head(data)

#Convert to a matrix and turn metabolites into rownames
dim(data)
data.mat <- as.matrix(data[,5:28])
row.names(data.mat) <- data$BIOCHEMICAL

#Remove rows where standard deviation = 0 
data.mat <- subset(data.mat, apply(data.mat,1, sd)!=0)

#Load strain info 
Sinfo <- read.csv("Sinfomin.csv",stringsAsFactors = FALSE)

#heatmap of everything 
library(gplots)
library(viridis)
library(scales)

StrainColors<-factor(Sinfo$Group)
levels(StrainColors)[levels(StrainColors)=="PLT_32B"]<-"#3B0F70FF"
levels(StrainColors)[levels(StrainColors)=="MEM"]<-"#000004FF"
levels(StrainColors)[levels(StrainColors)=="RH127O"]<-"#FE9F6DFF"
levels(StrainColors)[levels(StrainColors)=="SMC7758"]<-"#DE4968FF"
levels(StrainColors)[levels(StrainColors)=="VPI"]<-"#8C2981FF"

#Figure S7B
heatmap.2(log2(data.mat), trace="none",
          labRow = FALSE, labCol = FALSE,
          key.title = "",
          key.xlab = "Log 2 Relative Abundance",
          col = rev(redblue(100)),
          ColSideColors = as.character(StrainColors),
          sepcolor="white", sepwidth = c(0.05,0.05), colsep = c(5,10,15,19))

legend("center", title = "Strain",legend=c("MEM", "RH127O","CFPLTA003_2B", "SMC7757","SMC7758"), 
       fill=c("#000004FF","#FE9F6DFF","#3B0F70FF","#8C2981FF","#DE4968FF"), cex=0.8, box.lty=0)
dev.off()

##Create a PCA plot and get contributions data for specific genes##
H_PCA <- prcomp(t(data.mat), center = TRUE, scale = TRUE)
H_PCA

#PCA plots by different conditions
individual <- Sinfo[['Group']]
H_PCA_plot_individual <- ggbiplot::ggbiplot(pcobj = H_PCA, var.axes = FALSE, groups = individual, circle.prob = TRUE, ellipse = TRUE, circle = TRUE, obs.scale = 1, var.scale = 1, size = 2, choices = c(1,2)) +
  geom_point(aes(colour=individual), size = 2)+theme_grey(base_size = 22)

#Figure S7C (Left Panel)
H_PCA_plot_individual + 
  scale_color_manual(values=c("#000004FF","#3B0F70FF","#FE9F6DFF","#DE4968FF","#8C2981FF")) +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  labs(aesthetic='custom text') + 
  guides(color=guide_legend("Strain"))

#Figure S7C (Right panel)
H_PCA_plot_individual_2 <- ggbiplot::ggbiplot(pcobj = H_PCA, var.axes = FALSE, groups = individual, circle.prob = TRUE, ellipse = TRUE, circle = TRUE, obs.scale = 1, var.scale = 1, size = 2, choices = c(2,3)) +
  geom_point(aes(colour=individual), size = 2)+theme_grey(base_size = 22)

H_PCA_plot_individual_2 + 
  scale_color_manual(values=c("#000004FF","#3B0F70FF","#FE9F6DFF","#DE4968FF","#8C2981FF")) +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  labs(aesthetic='custom text') + 
  guides(color=guide_legend("Strain"))
