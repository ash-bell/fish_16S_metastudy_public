library(dada2)
library(phyloseq)
library(Biostrings)
library(tidyverse)

df <- read.csv("WKDIR/BIOPROJECT/BIOPROJECT_blastout.txt", sep="\t", header=F)

Mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}

ps <- readRDS(file = "WKDIR/BIOPROJECT/BIOPROJECT_PS_object.rds")

x <- df %>%
  group_by(V1) %>%
  summarise(start = Mode(V3), end = Mode(V4)) %>%
  full_join(data.frame(V1 = names(refseq(ps)), 
                       width = width(refseq(ps))), by = "V1") %>%
  mutate(start = replace_na(start, Mode(start)), 
         end = replace_na(end, Mode(end)),
         trim = end - start, 
         start = replace(start, trim < 248, Mode(start)),
         end = replace(end, trim < 248, Mode(end)))

x[x$end > x$width, "end"] <- x[x$end > x$width, "width"]
x[x$start > x$end, "end"] <- x[x$start > x$end, "width"]

trm_refseq <- subseq(refseq(ps)[x$V1], x$start, x$end)
trm_refseq <- trm_refseq[width(trm_refseq) %in% (248:256), ]
trm_otu_table <- otu_table(ps)[, colnames(otu_table(ps)) %in% names(trm_refseq)]

DNA = as.data.frame(trm_refseq)
taxa <- assignTaxonomy(unique(DNA$x), "WKDIR/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "WKDIR/silva_species_assignment_v138.1.fa.gz")

counts <- t(otu_table(ps))
x = merge(DNA, counts, by="row.names", all=T)
x$Row.names <- NULL
y = aggregate(. ~ x, data=x, sum, na.rm=F)
row.names(y) <- y$x
y$x <- NULL
ps <- merge_phyloseq(otu_table(y, taxa_are_row = T), tax_table(taxa))

ps <- subset_taxa(ps, Family != "Mitochondria")
ps <- subset_taxa(ps, Order != "Chloroplast")
ps <- subset_taxa(ps, Kingdom != "Eukaryota")

### prevalence filter ###
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps), MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps), tax_table(ps))

# Define prevalence threshold as present ASV in 2 samples or more
prevalenceThreshold <- 2
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps <- prune_taxa(keepTaxa, ps)

# Add metadata
metadata <- read.csv(Sys.glob("BIOPROJECT_metadata*"))
rownames(metadata) <- paste0(metadata$BioProject, "_", metadata$Run)
sample_names(ps) <- sapply(strsplit(basename(sample_names(ps)), "_2"), `[`, 1)
ps <- merge_phyloseq(sample_data(metadata), ps)

saveRDS(ps, file = "WKDIR/BIOPROJECT/BIOPROJECT_PS_new.rds")
