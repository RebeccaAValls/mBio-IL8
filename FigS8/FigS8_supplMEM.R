#read in the file
data <- read.csv("SupplMEM.csv")
head(data)

#Convert to a matrix and turn metabolites into rownames
dim(data)
data.mat <- as.matrix(data[,2:51])
row.names(data.mat) <- data$Metabolite

#Remove rows where standard deviation = 0 
data.mat <- subset(data.mat, apply(data.mat,1, sd)!=0)

#Load strain info 
Sinfo <- read.csv("Sinfo.csv",stringsAsFactors = FALSE)

#heatmap of everything 
library(gplots)
library(viridis)
library(scales)

show_col(viridis_pal(option = "magma")(10))

StrainColors<-factor(Sinfo$Group)
levels(StrainColors)[levels(StrainColors)=="MEM"]<-"#000004FF"
levels(StrainColors)[levels(StrainColors)=="PLT_21B"]<-"#FCFDBFFF"
levels(StrainColors)[levels(StrainColors)=="PLT_32B"]<-"#451077FF"
levels(StrainColors)[levels(StrainColors)=="PLT_41B"]<-"#721F81FF"
levels(StrainColors)[levels(StrainColors)=="AD126T3B"]<-"#180F3EFF"
levels(StrainColors)[levels(StrainColors)=="HAP2B"]<-"#F1605DFF"
levels(StrainColors)[levels(StrainColors)=="RH127O"]<-"#FD9567FF"
levels(StrainColors)[levels(StrainColors)=="SMC7758"]<-"#CD4071FF"
levels(StrainColors)[levels(StrainColors)=="TL139C3B"]<-"#FEC98DFF"
levels(StrainColors)[levels(StrainColors)=="VPI"]<-"#9F2F7FFF" 

##Create a PCA plot and get contributions data for specific genes##
H_PCA <- prcomp(t(data.mat), center = TRUE, scale = TRUE)
H_PCA

#Figure S8A
heatmap.2(log2(data.mat), trace="none",
          key.title = "",
          labRow = FALSE, labCol = FALSE,
          key.xlab = "Log 2 Relative Abundance",
          cexRow =.2, cexCol = .6,
          col = rev(redblue(100)),
          ColSideColors = as.character(StrainColors),
          sepcolor="white", sepwidth = c(0.05,0.05), colsep = c(4,7,9,12,15,18,23,28,32,35,39,40,46))

legend("center", title = "Strain",legend=c("sMEM","CFPLTA003_2B","CFPLTA004_1B",
                                           "AD126T_3B","HAP130N_2B","RH127O","ATCC_29741",
                                           "TL139C_3B","VPI","CFPLTA002_1B"), 
       fill=c("#000004FF","#451077FF","#721F81FF","#180F3EFF",
              "#F1605DFF","#FD9567FF","#CD4071FF","#FEC98DFF","#9F2F7FFF","#FCFDBFFF"), cex=0.8, box.lty=0)
dev.off()

#Figure S8B
Species <- Sinfo[['Species']]
H_PCA_plot_Species <- ggbiplot::ggbiplot(pcobj = H_PCA, var.axes = FALSE, groups = Species, circle.prob = TRUE, ellipse = TRUE, circle = TRUE, obs.scale = 1, var.scale = 1, size = 2, choices = c(1,2)) +
  geom_point(aes(colour=Species), size = 2)+theme_grey(base_size = 22)
H_PCA_plot_Species + theme_bw() + scale_color_viridis(discrete=TRUE, 
                                                      name="Species",
                                                      labels = c("B. fragilis", 
                                                                 "B. thetaiotaomicron", 
                                                                 "B. vulgatus", 
                                                                 "sMEM"))

