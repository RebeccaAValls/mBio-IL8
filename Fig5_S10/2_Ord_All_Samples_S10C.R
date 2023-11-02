library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")
library("vegan")
library("viridis")

###############################################################################
##merge metadata with phyloseq object##
#load phyloseq object with counts and ASVs
load("ASV_physeq.rda")
#load sample info metadata
sample_info <- read.csv("Minfo.csv", header=T, row.names=1)
#convert metadata to phyloseq format
sample_info_tab_phy <- sample_data(sample_info)
#add metadata to phyloseq object
physeq1 = merge_phyloseq(ASV_physeq, sample_info_tab_phy)
###############################################################################
##data cleanup##
#remove NAs
physeq2 <- subset_taxa(physeq1, !is.na(Phylum))

#any taxa with zero counts?
any(taxa_sums(physeq2) == 0)
sum(taxa_sums(physeq2) == 0)
#removes taxa with zero counts 
physeq2 = prune_taxa(taxa_sums(physeq2) > 0, physeq2)

###############################################################################
#remove control and baseline samples 

physeq2 <- subset_samples(physeq2, Condition=="Plus" 
                          | Condition=="Minus" 
                          | Condition=="Ctrl")
physeq2

###############################################################################
#some graphs of physeq2 object with logged counts-per-sample
qplot(log10(rowSums(otu_table(physeq2))),binwidth=0.2) +
  xlab("Logged counts") + ylab("Number of ASVs") + theme_bw()


# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf = apply(X = otu_table(physeq2),
               MARGIN = ifelse(taxa_are_rows(physeq2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq2),
                    tax_table(physeq2))

#a plot of ASV taxa abundances 
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq2, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq2),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 1% of total samples
prevalenceThreshold = 0.01 * nsamples(physeq2)
prevalenceThreshold

# Remove taxa with <1% prevalence
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
physeq3 = prune_taxa(keepTaxa, physeq2)

#exclude low count phyla
filterPhyla = c("Armatimonadota", "Chloroflexi", 
                "Deinococcota", "Euryarchaeota", 
                "Fusobacteriota", "Halobacterota",
                "Planctomycetota", "Spirochaetota",
                "Verrucomicrobiota")
physeq3 = subset_taxa(physeq3, !Phylum %in% filterPhyla)

#exclude contamination 
filterFam = c("Mitochondria")
physeq3 = subset_taxa(physeq3, !Family %in% filterFam)


#some graphs of ASV_physeq object with logged counts-per-sample
qplot(log10(rowSums(otu_table(physeq3))),binwidth=0.2) +
  xlab("Logged counts") + ylab("Number of ASVs") + theme_bw()

# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf = apply(X = otu_table(physeq3),
               MARGIN = ifelse(taxa_are_rows(physeq3), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq3),
                    tax_table(physeq3))

#a plot of ASV taxa abundances 
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq3, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq3),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme(legend.position="none") +
  theme_bw() + guides(color = FALSE, size = FALSE)

##############################################################################################################################################
#Ordination
#reduce dimensions by combining at genera level 
ASV_physeq_genus <- tax_glom(physeq3, taxrank="Genus")

#convert to a deseq object
deseq_counts <- phyloseq_to_deseq2(ASV_physeq_genus, ~OrdAll)

#variance stabilizing transformation; accounts for different depths
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

#convert back to phyloseq object#
# pull out transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
#making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
#threshold neg values up to 0
vst_count_phy[vst_count_phy < 0.0] <- 0.0
#make metadata sheet
sample_info_tab_phy <- sample_data(physeq3)
#make the phyloseq object
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

#######################
##distance and ordination##
#bray distance
dist = phyloseq::distance(vst_physeq, method="bray")

#ordination 
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="bray")

##plotting##
#with circles
OrdPMonth <- plot_ordination(vst_physeq, vst_pcoa, color = "OrdAll") + 
  geom_point(size=1)

Ord <- OrdPMonth + 
  stat_ellipse(type = "t") +
  theme_bw() + 
  #ggtitle("MDS Ordination") +
  labs(color="Samples") +
  scale_color_viridis(discrete=TRUE, labels = c("Baseline","Stool Pools", 
                                                "Early (-) Bacteroides",
                                                "Early (+) Bacteroides",
                                                "Late (-) Bacteroides", 
                                                "Late (+) Bacteroides")) +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14, face="bold"),
        plot.title = element_text(size =16, face = "bold"))

Ord
