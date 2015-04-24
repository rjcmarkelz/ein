#############
# Cody Markelz
# markelz@gmail.com
# voom/limma analysis pipeline for brassica ein mutants
#############

# go to data directory
setwd("/Users/Cody_2/Box Sync/Brassica (rjmarkelz@ucdavis.edu)")
B3_40_9_counts <- read.delim("Ein_v1.5_B3_ein9_ein40.tsv",
                     header = TRUE, sep = "\t")
dim(B3_40_9_counts)
#[1] 41612    58

B3_194_counts <- read.delim("Ein_v1.5_B3_ein194_renamed.tsv",
                     header = TRUE, sep = "\t")
dim(B3_194_counts)
# [1] 36138    13


#replace all NA values with 0 
B3_40_9_counts[is.na(B3_40_9_counts)] <- 0
head(B3_40_9_counts)
tail(B3_40_9_counts)

B3_194_counts[is.na(B3_194_counts)] <- 0
head(B3_194_counts)
tail(B3_194_counts)

#remove first row
B3_40_9_counts <- B3_40_9_counts[-1,]
head(B3_40_9_counts)[,1:10]

B3_194_counts <- B3_194_counts[-1,]
head(B3_194_counts)[,1:10]

# make gene names rownames
rownames(B3_40_9_counts) <- B3_40_9_counts[,1]
B3_40_9_counts <- B3_40_9_counts[,-1]
head(B3_40_9_counts)[,1:10]
dim(B3_40_9_counts)

rownames(B3_194_counts) <- B3_194_counts[,1]
B3_194_counts <- B3_194_counts[,-1]
head(B3_194_counts)
dim(B3_194_counts)

# clean up column names for subsetting and model building

B3_40_coln <- as.data.frame(colnames(B3_40_9_counts))
colnames(B3_40_coln) <- paste("names")
B3_40_coln

B3_40_coln$names <- sub("(JD015_)(.+)(.sorted)(.+)", "\\2", B3_40_coln$names)
B3_40_coln
B3_40_coln$names <- sub("(.+)(ein)(_)(.+)(.so)(.+)", "\\2\\4", B3_40_coln$names)
B3_40_coln
B3_40_coln$names <- sub("(B2)(.+)", "B3\\2", B3_40_coln$names)
B3_40_coln
str(B3_40_coln)
colnames(B3_40_9_counts) <- paste(as.character(B3_40_coln$names))
head(B3_40_9_counts)


B3_194_coln <- as.data.frame(colnames(B3_194_counts))
colnames(B3_194_coln) <- paste("names")
B3_194_coln

B3_194_coln$names <- sub("(B2)(.+)", "B3\\2", B3_194_coln$names)
B3_194_coln$names <- sub("E", "e", B3_194_coln$names)
B3_194_coln$names <- sub("INT", "Internode", B3_194_coln$names)
B3_194_coln$names <- sub("(.+)(.sorted.bam)", "\\1", B3_194_coln$names)
B3_194_coln
str(B3_194_coln)
colnames(B3_194_counts) <- paste(as.character(B3_194_coln$names))
head(B3_194_counts)

# use these cleaned tables for analysis
# note missing col.name for gene name column in these output tables
setwd("/Users/Cody_2/Box Sync/Brassica (rjmarkelz@ucdavis.edu)")
write.table(B3_194_counts, file="B3_194_counts.csv", sep=",", col.names = NA) 
write.table(B3_40_9_counts, file="B3_40_9_counts.csv", sep=",", col.names = NA) 

# merge dataframes to do joint analysis on all internode data
dim(B3_194_counts)
dim(B3_40_9_counts)


B3_194_40_9 <- merge(B3_194_counts, B3_40_9_counts, by = "row.names", all = T)
dim(B3_194_40_9)
head(B3_194_40_9)
tail(B3_194_40_9)

B3_194_40_9[is.na(B3_194_40_9)] <- 0
head(B3_194_40_9)
tail(B3_194_40_9)

rownames(B3_194_40_9) <- B3_194_40_9[,1]
# remove row name column
B3_194_40_9 <- B3_194_40_9[,-1]
head(B3_194_40_9)[,1:10]
dim(B3_194_40_9)


setwd("/Users/Cody_2/Box Sync/brassica")
save.image("brassica_ein.RData")


# load libs
library(edgeR)
library(limma)

dim(B3_40_9_counts)
head(B3_40_9_counts)
Br_group <- factor(sub("(.+)(_)(.+)(_)(S.+)(_)(\\d)",
                       "\\1\\2\\3\\4\\5", colnames(B3_40_9_counts)))
head(Br_group)
Br_group
Br_trt <- factor(sub("(.+)(_)(.+)(_)(S.+)(_)(\\d)",
                       "\\5", colnames(B3_40_9_counts)))
Br_geno <- factor(sub("(.+)(_)(.+)(_)(S.+)(_)(\\d)",
                       "\\1", colnames(B3_40_9_counts)))
Br_tissue <- factor(sub("(.+)(_)(.+)(_)(S.+)(_)(\\d)",
                       "\\3", colnames(B3_40_9_counts)))


group_design <- model.matrix(~ 0 + Br_group)


Br_trt <- relevel(Br_trt, ref = "Sun")
Br_trt

Br_geno <- relevel(Br_geno, ref = "B3")
Br_geno

Br_tissue
Br_tissue <- relevel(Br_tissue, ref = "Leaf")

ref_design <- model.matrix(~ Br_geno*Br_trt*Br_tissue)

brassica_DE <- DGEList(counts = B3_40_9_counts, group = Br_group)
brassica_DE$samples
dim(brassica_DE)
# [1] 41611    57

# keep genes with at least 1 count per million in at least 20 samples
brassica_DE <- brassica_DE[rowSums(cpm(brassica_DE) > 1 ) >= 20,]
dim(brassica_DE)
# [1] 25351    57

brassica_DE <- calcNormFactors(brassica_DE)

# output plots, they are much to large to fit into memory
setwd("/Users/Cody_2/Box Sync/brassica_ein")
png(file="Brassica_MDS.png",width=1000,height=1000,res=100)
plotMDS(brassica_DE)
dev.off()

system.time(brass_voom <- voom(brassica_DE, group_design, plot = TRUE))
head(brass_voom)


system.time(fit1 <-lmFit(brass_voom, ref_design))

fit1 <- eBayes(fit1)
toptable(fit1)
head(fit1)











