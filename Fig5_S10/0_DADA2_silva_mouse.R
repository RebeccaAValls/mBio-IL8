#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") 
# change the ref argument to get other versions
#Log into server, and copy this file and RDP file into the folder with all of the samples you want to look at 
#Open R by typing R into the terminal 
#install DADA2 & DECIPHER if they are not already installed 
#Quit R
#type: nohup Rscript DADA2_DECIPHER_RDP.R &
#This will run the script - it will take a while; I recommend testing on a subset of files before doing all of them

#if (!requireNamespace("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
#BiocManager::install("dada2")

library("dada2"); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")

#Set the filepath to your working directory
path <- getwd()
list.files(path)

#Sort forward file names
fnFs <- sort(list.files(path, pattern="out.R1.fastq", full.names = TRUE))
head(fnFs)
length(fnFs)

#This will return a single sample name for each pair of Forward and Reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_out"), `[`,1)
head(sample.names)
length(sample.names)

#Visualize quality profiles of forward & reverse  reads
pdf(file = "QualF.pdf",
    width = 10,
    height = 10)
plotQualityProfile(fnFs[1:10])
dev.off()


# Place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

head(filtFs)

#trim Forward reads at the 220th bp and rev reads at the 160th
out <- filterAndTrim(fnFs, filtFs, truncLen=c(220),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#This takes a little <10min each
#Option to increase nbases parameter
errF <- learnErrors(filtFs, multithread=TRUE)

#Check learnErrors function by plotting
#Both of these produce a warning message: Transformation introduced infinite values in continuous y-axis
pdf(file = "errF.pdf",
    width = 10,
    height = 10)
plotErrors(errF, nominalQ=TRUE)
dev.off()

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)

#Inspecting returned data-class object
dadaFs[[1]]
 
#Construct sequence table 
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
write.csv(seqtab, "seqtab.csv")

# Inspect distribution of sequence lengths
Tseq <- table(nchar(getSequences(seqtab)))
Tseq
write.csv(Tseq, "Tseqtab.csv")

#Remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track.csv")

##DADA2 ##Silva
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
taxa.species <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz", tryRC =TRUE)

taxa.print <- taxa.species # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#prepping files for phyloseq
ASV_physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               tax_table(taxa.species))

dna <- Biostrings::DNAStringSet(taxa_names(ASV_physeq))
names(dna) <- taxa_names(ASV_physeq)
ASV_physeq <- merge_phyloseq(ASV_physeq, dna)
taxa_names(ASV_physeq) <- paste0("ASV", seq(ntaxa(ASV_physeq)))
ASV_physeq

save(ASV_physeq, file = "ASV_physeq.rda")





