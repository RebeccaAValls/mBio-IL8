#packages
library("phyloseq"); packageVersion("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("gplots")
library("ggpubr")
library("dplyr")

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
#convert to relative abundance
physeq2 = transform_sample_counts(physeq2, function(x) (x / sum(x))*100 )

###############################################################################################################
#make tables
sample_info_tab <- as.data.frame(sample_data(physeq2))

count_tab <- as.data.frame(otu_table(physeq2))
count_tab <- t(count_tab)

phyla_counts_tab <- otu_table(tax_glom(physeq2, taxrank="Genus"))
phyla_counts_tab <- t(phyla_counts_tab)

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(physeq2, taxrank="Genus"))[,6]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#turn NAs into their own group and add this row to our phylum count table:
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(phyla_and_unidentified_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Day"=sample_info_tab$Day,
                                  "Sample"=row.names(sample_info_tab), 
                                  "Condition"=sample_info_tab$Condition, 
                                  "Experiment"=sample_info_tab$Experiment,
                                  stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#check the dimensions of this table at this point
dim(phyla_and_unidentified_counts_tab)

# here, we'll only keep rows (taxa) that make up greater than 10% in any
temp_filt_major_taxa_proportions_tab <- as.data.frame(phyla_and_unidentified_counts_tab[apply(phyla_and_unidentified_counts_tab, 1, max) > 10, ])

# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
rownames(temp_filt_major_taxa_proportions_tab)

# though each of the filtered taxa made up less than 10% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(phyla_and_unidentified_counts_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(filt_major_taxa_proportions_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
head(major_taxa_for_plot.g2)

###############################################################################
write.csv(major_taxa_for_plot.g2, "1_mouse_RelAbund.csv")
###############################################################################

