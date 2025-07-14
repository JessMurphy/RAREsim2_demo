
library(dplyr)

setwd("C:/Users/murphjes/Documents/Repositories/RAREsim2_demo")

pop = "AFR"
n = 8128
maf = 0.01

genes = c("ADGRE5", "DDX39A", "PKN1", "PTGER1", "GPIC1", "DNAJB1", "TECR", "NDUFB7", "CLEC17A", "ADGRE3", "ZNF333", "ADGRE2")
start = c(14491313, 14519631, 14543865, 14583278, 14588572, 14625582, 14627897, 14676890, 14693896, 14729929, 14800613, 14843205)
end = c(14519537, 14530192, 14582679, 14586174, 14606944, 14640582, 14676792, 14682874, 14721969, 14800839, 14844558, 14889353)

gnomad = read.table(paste0("./input/gnomad/gnomad.exomes.r2.1.1.sites.19.block37_", pop, ".vcf.INFO"), sep='\t', header=T) %>% rename(position = POS)
names(gnomad)[5:8] = sapply(strsplit(names(gnomad)[5:8], "_", fixed=T), head, 1)

master = read.table(paste0("./input/chr19.block37.", pop, ".master.legend"), sep='\t')
colnames(master) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
master$alleles = paste0(master$a0, '/', master$a1)

gnomad2 = gnomad %>% group_by(position) %>% slice_max(order_by=AC, n=1, with_ties=F) %>% 
  filter(position %in% master$position, AC!=0)

### SUBSET MASTER ###

# subset the file according to the number of alleles at each position
singles = master %>% filter(prob==1)
dups = master %>% filter(prob==0.5)
trips = master %>% filter(prob==".")

# get a list of all the  positions in trips
trips.po = levels(droplevels(as.factor(trips$position)))

# create a table of transitions/transversions
trips.po1  =  as.data.frame(matrix(NA, nrow=length(trips.po), ncol=3))
colnames(trips.po1) = c('position', 'draw', 'TiTv')
trips.po1$position = trips.po
trips.po1$draw = runif(nrow(trips.po1))
trips.po1$TiTv[which(trips.po1$draw < 0.7396)] = 'transition'
trips.po1$TiTv[which(trips.po1$draw >= 0.7396)] = 'transversion'
#head(trips.po1)

# subset the transitions
ti = trips.po1 %>% filter(TiTv  == 'transition')
ann.ti = trips %>% filter(position %in% ti$position, alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C')

# remove all of the transitions
ann.tv = trips %>% filter(!(position %in% ti$position), !(alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C'))

# randomly pick an allele from the transversions
ann.tv2 = ann.tv %>% group_by(position) %>% sample_n(size=1)

# merge transitions and transversions
trips2 = union(ann.ti, ann.tv2)

# randomly pick an allele from the duplicates
dups2 = dups %>% group_by(position) %>% sample_n(size=1)

# merge all positions
master2 = union(singles, union(dups2, trips2)) %>% arrange(position) %>% select(id, position, a0, a1, AC, prob, exonic, gene, fun)

master2$fun = ifelse(master2$fun=="synonymous SNV", "syn", "fun")
master2$exonic[grepl("exonic", master2$exonic)] = "exonic"
master2$gene[grepl("ZNF333", master2$gene)] = "ZNF333"

master2.filtered = master2 %>% filter(AC!=".") %>% mutate(AC=as.numeric(AC)) %>% filter(AC>0) %>%
  mutate(common=ifelse(AC/(2*n) > maf & AC/(2*n) < 1-maf, "common", "rare"))

sum.vars.master = master2.filtered %>% group_by(gene, fun, common) %>% summarize(n=n()) %>%
  filter(common=="rare", fun=="fun")


gnomad2.merged = merge(gnomad2, master2, by="position") %>%
  mutate(common=ifelse(AC.x/(2*n) > maf & AC.x/(2*n) < 1-maf, "common", "rare"))

sum.vars.gnomad = gnomad2.merged %>% group_by(gene, fun, common) %>% summarize(n=n()) %>%
  filter(common=="rare", fun=="fun")


