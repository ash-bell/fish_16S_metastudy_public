library(dada2)
library(phyloseq)

path <- "WKDIR/BIOPROJECT/reads"

fnFs <- sort(list.files(path, pattern = "^BIOPROJECT_.*_1.QC.trm.fq.gz", full.names = T))
fnRs <- sort(list.files(path, pattern = "^BIOPROJECT_.*_2.QC.trm.fq.gz", full.names = T))

sample.names <- sapply(strsplit(basename(fnFs), "_1"), `[`, 1)
filtFs <- file.path(path, paste0(sample.names, "_1.filt.QC.trm.fq.gz"))
filtRs <- file.path(path, paste0(sample.names, "_2.filt.QC.trm.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=T, compress=T, multithread=T)

filtFs <- Sys.glob(paste0(path, "/BIOPROJECT_*_1.filt.QC.trm.fq.gz"))
filtRs <- Sys.glob(paste0(path, "/BIOPROJECT_*_2.filt.QC.trm.fq.gz"))

names(filtFs) <- sapply(strsplit(basename(filtRs), "_1"), `[`, 1)
names(filtRs) <- sapply(strsplit(basename(filtRs), "_2"), `[`, 1)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(248,999)]

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

ps <- otu_table(seqtab.nochim, taxa_are_rows=F)
sample_names(ps) <- sapply(strsplit(basename(sample_names(ps)), "_2"), `[`, 1)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("BIOPROJECT_ASV", seq(ntaxa(ps)))

Biostrings::writeXStringSet(refseq(ps), "BIOPROJECT_ASVs.fna", append=FALSE, compress = F, compression_level = NA, format = "fasta")
saveRDS(ps, file = "BIOPROJECT_PS_object.rds")
