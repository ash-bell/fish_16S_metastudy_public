library(phyloseq)
library(Biostrings)
library(decontam)

path <- "WKDIR"

path <- "/lustre/home/bt273/BIOSSCOPE/metag/ashley/fish_16S"

bioprojects <- read.table("BIOPROJECT_list.txt")
bioprojects <- bioprojects$V1

list_ps <- sort(Sys.glob(paste0(path,"/*/*_PS_new.rds")))
all_ps <- lapply(list_ps, readRDS)

# if a project has less than 20 ASVs, remove it
for (i in 1:length(all_ps)) { if (nrow(otu_table(all_ps[[i]])) < 20) { all_ps[[i]] = NULL } }

# merge phyloseqs auto merges taxa that are the same if the are ASV sequences
ps <- merge_phyloseq(all_ps[[1]], all_ps[[2]])
for (i in 3:length(all_ps)) { 
  ps <- merge_phyloseq(all_ps[[i]], ps)
}


contam <- isContaminant(ps, neg="is.neg", threshold=0.5,
                        detailed=TRUE, normalize=TRUE,
                        method='prevalence')

ps <- prune_taxa(!contam$contaminant, ps)
ps <- prune_samples(sample_data(ps)$is.neg == FALSE, ps)

### Remove samples with low and high read counts as failed sequencing runs / low depth, too much ##
## find how many reads does each sample have and what quantile is 2.5% and 97.5 (keep middle 95%)
df <- colSums(otu_table(ps)) %>%
  as.data.frame() %>%
  dplyr::rename(count = ".")

# Remove samples with 0 reads
df <- df %>%
  filter(count > 0)

low_reads <- round(quantile(df$count, probs = 0.025))
cat("2.5 percentile is",low_reads[[1]])

high_reads <- round(quantile(df$count, probs = 0.975))
cat("97.5 percentile is",high_reads[[1]])

# remove samples with less 76 reads and more than per 60094 reads sample (95 percentile)
ps <- prune_samples(sample_sums(ps) >= low_reads[[1]], ps)
ps <- prune_samples(sample_sums(ps) <= high_reads[[1]], ps)

# also remove taxa that have no counts in PS obj after pruning these samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Take out the metadata so we can get query host lineages from NCBI plus curate in excel
# Make a column called "scientific_name" with scientific name and taxID of each fish
# Using the rownames from the exported metadata table, for each paper find all 
#the environmental parameters and re-join that to the metadata table to be added 
#in as extra columns. 
write.csv(sample_data(ps), "data/raw_metadata.csv")

# read in lineage of host taxid, split columns by taxonomic level
host_lineage <- read.csv("data/host_summarised_lineage.csv", col.names = c("host_taxid", "summarised_lineage"))
lineage <- c("superkingdom","phylum","class","order","family","genus","species")
lineage <- paste0("host_", lineage)
host_lineage <- host_lineage %>% 
  separate(summarised_lineage, lineage, ";") %>%
  select(all_of(lineage), host_taxid)

# Using the rownames from the exported metadata table, for each paper manually find all 
#the environmental parameters and re-join that to the metadata table to be added 
#in as extra columns. 
curated.metadata <- read_csv("data/curated_metadata.csv")

# join new host tax levels to ps object
metadata <- left_join(metadata, host_lineage, by = "host_taxid")
rownames(metadata) <- paste0(metadata$BioProject,"_",metadata$Run)
sample_data(ps) <- metadata

# remove any samples with no host_species 
ps <- prune_samples(!is.na(sample_data(ps)$host_species), ps)

# next filter for the top ASVs in each fish species. Only retain ASVs per species 
# across studies if ASV is present in top 95 percentile of studies.
# example: if 100 samples of fish species X, only keep ASV is present in at least 5 samples

# make first PS obj of first fish species to merge into
ps.new <- subset_samples(ps, host_species == unique(sample_data(ps)$host_species)[1])
wh0 <- genefilter_sample(ps.new, filterfun_sample(function(x) x > 1), 
                         A=0.05*nsamples(ps.new))
ps.new <- prune_taxa(wh0, ps.new)

for (spp in unique(sample_data(ps)$host_species)[-1]) {
  print(spp)
  tmp <- subset_samples(ps, host_species == spp)
  wh0 <- genefilter_sample(tmp, filterfun_sample(function(x) x > 1), 
                           A=0.05*nsamples(tmp))
  tmp <- prune_taxa(wh0, tmp)
  ps.new <- merge_phyloseq(ps.new, tmp)
}

# remove any taxa no longer present in any samples and samples with no ASVs
ps.new <- prune_taxa(rowSums(otu_table(ps.new)) > 0, ps.new)
ps.new <- prune_samples(sample_sums(ps) > 0, ps.new)

# get phylip tree of fish host taxonomy to determine phylogeny of fish
host_tree <- read.tree("data/summarised_lineage_phylip_tree.phy")
species_order <- gsub("\'", "", get_taxa_name(ggtree(host_tree)))

# order fish species by the order they are in taxnomic tree (ordered factor), so closer related fish are
# next to each other in plots
sample_data(ps.new)$host_species <- factor(sample_data(ps.new)$host_species, levels = species_order)
ordered <- with(sample_data(ps.new), sample_data(ps.new)[order(host_species), ])
ordered$host_order <- factor(ordered$host_order, levels = unique(ordered$host_order))
sample_data(ps.new) <- ordered

### how many samples per study made it though QC
metadata <- sample_data(ps.new) %>%
  data.frame()

table(metadata$BioProject)

#### export the refseq data so we can make a tree
dna <- Biostrings::DNAStringSet(taxa_names(ps))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(refseq(dna), ps)

saveRDS(ps, "all_studies_phyloseq.rds")
Biostrings::writeXStringSet(refseq(ps), "all_studies_ASVs.fna", 
                            append=FALSE, compress = F, compression_level = NA, 
                            format = "fasta")