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
physeq2 <- subset_samples(physeq2, Day=="Day2" | Day=="Day4")
#physeq2 <- subset_samples(physeq2, Day=="Day12" | Day=="Day13")
#physeq2 <- subset_samples(physeq2, Experiment=="2" | Experiment=="3")



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
filterPhyla = c("Armatimonadota","Campylobacterota", "Chloroflexi", 
                "Deferribacterota", "Deinococcota", "Desulfobacterota",
                "Euryarchaeota", "Fusobacteriota", "Halobacterota",
                "Patescibacteria", "Planctomycetota", "Spirochaetota",
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

##########################################################################################################
#Deseq2 to figure out which ASVs are causing significance

#This will model on the continuous variable Age_Years
ASV_deseq <- phyloseq_to_deseq2(physeq3, ~Condition)

#convert to deseq object
#workaround to deal with 0s
cts <- counts(ASV_deseq)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq <- estimateSizeFactors(ASV_deseq, geoMeans=geoMeans)

#deseq standard analysis
ASV_deseq <- DESeq(ASV_deseq)

#pull out results table
resultsNames(ASV_deseq)
results(ASV_deseq, name = "Condition_Plus_vs_Minus", alpha=0.01)
deseq_res <- results(ASV_deseq, name = "Condition_Plus_vs_Minus", alpha=0.01)

#look at results table
summary(deseq_res)

# subset table for significance
deseq_res <- deseq_res[which(deseq_res$padj < 0.01), ]

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_tax <- cbind(as(deseq_res, "data.frame"), as(tax_table(ASV_physeq)[row.names(deseq_res), ], "matrix"))

#sort by base mean 
sum_PlusvMinus <- deseq_res_tax[order(deseq_res_tax$baseMean, decreasing=T), ]

#visualization
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#make a graph
x = tapply(deseq_res_tax$log2FoldChange, deseq_res_tax$Phylum, function(x) max(x))
x = sort(x, TRUE)
deseq_res_tax$Phylum = factor(as.character(deseq_res_tax$Phylum), levels=names(x))

# Genus order
x = tapply(deseq_res_tax$log2FoldChange, deseq_res_tax$Genus, function(x) max(x))
x = sort(x, TRUE)

#Make graph
gen_Early <- ggplot(deseq_res_tax, aes(x=log2FoldChange, y=Genus, color=Phylum)) + 
  geom_point(size=6) +
  ggtitle("Early") +
  xlab("Log2 Fold Change") +
  ylab("Genus") +
  geom_vline(xintercept=0, color = "black", size=1.5) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        plot.title = element_text(size =16, face = "bold")) + 
  scale_color_viridis(discrete=TRUE)
gen_Early


