setwd('/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/')
library(dplyr)
library(ggplot2)
library(readxl)
library(here)

nmf = read.csv(here('code', 'revision','gc_nmf_gene-weights.csv', row.names=1))
nmf = as.matrix(nmf)

topN = lapply(c("nmf5","nmf14","nmf26"), function(x) names(sort(nmf[,x], decreasing=T))[1:300])
names(topN) = c("nmf5", "nmf14","nmf26")
sapply(topN, length) #150 in all


infant = read_excel('Supplementary_table_S3.xlsx', sheet =1)
teen = read_excel('Supplementary_table_S3.xlsx', sheet =2)
adult = read_excel('Supplementary_table_S3.xlsx', sheet =3)
elderly = read_excel('Supplementary_table_S3.xlsx', sheet =4)


colnames(infant)[6] <- "log2.fold.change"
colnames(teen)[6] <- "log2.fold.change"
colnames(adult)[6] <- "log2.fold.change"
colnames(elderly)[6] <- "log2.fold.change"


#many of the top weighted genes are shared
UpSetR::upset(UpSetR::fromList(topN))
#many of the top weighted genes are shared (without needing to install UpSetR)
t1 = table(unlist(topN))
common.weighted.genes = names(t1)[t1==3]
#the shared top weighted genes are also significant infant markers
nrow(filter(infant, gene_name %in% common.weighted.genes)) #13
nrow(filter(teen, gene_name %in% common.weighted.genes)) #2
nrow(filter(adult, gene_name %in% common.weighted.genes)) #0
nrow(filter(elderly, gene_name %in% common.weighted.genes)) #6

#restrict to only the top N genes that are unique to each NMF
topUnique = lapply(c("nmf5","nmf14","nmf26"), function(x) 
  setdiff(topN[[x]], unlist(topN[setdiff(names(topN),x)]))
)
names(topUnique) <- c("nmf5","nmf14","nmf26")
sapply(topUnique, length) #about 30 each for top 150
#40-45 each for top 300


infant.df = do.call(rbind, lapply(c("nmf5","nmf14","nmf26"), function(x) {
  tmp = infant[infant$gene_name %in% topUnique[[x]],]
  mutate(tmp[,c("gene_name", "log2.fold.change")], nmf=x)
})) %>% mutate(sign.logFC = sign(log2.fold.change)) %>%
  group_by(nmf, sign.logFC) %>% tally() %>% mutate(age="infant")


teen.df = do.call(rbind, lapply(c("nmf5","nmf14","nmf26"), function(x) {
  tmp = teen[teen$gene_name %in% topUnique[[x]],]
  mutate(tmp[,c("gene_name", "log2.fold.change")], nmf=x)
})) %>% mutate(sign.logFC = sign(log2.fold.change)) %>%
  group_by(nmf, sign.logFC) %>% tally() %>% mutate(age="teen")


adult.df = do.call(rbind, lapply(c("nmf5","nmf14","nmf26"), function(x) {
  tmp = adult[adult$gene_name %in% topUnique[[x]],]
  mutate(tmp[,c("gene_name", "log2.fold.change")], nmf=x)
})) %>% mutate(sign.logFC = sign(log2.fold.change)) %>%
  group_by(nmf, sign.logFC) %>% tally() %>% mutate(age="adult")

elderly.df = do.call(rbind, lapply(c("nmf5","nmf14","nmf26"), function(x) {
  tmp = elderly[elderly$gene_name %in% topUnique[[x]],]
  mutate(tmp[,c("gene_name", "log2.fold.change")], nmf=x)
})) %>% mutate(sign.logFC = sign(log2.fold.change)) %>%
  group_by(nmf, sign.logFC) %>% tally() %>% mutate(age="elderly")

ggplot(bind_rows(infant.df, teen.df, adult.df, elderly.df) %>%
         mutate(age=factor(age, levels=c("infant","teen","adult","elderly"))), 
       aes(x=factor(nmf, levels=c("nmf26", "nmf5", "nmf14")), y=n, fill=factor(sign.logFC)))+
  geom_bar(stat="identity", position=position_dodge(preserve="single"), color="black")+
  facet_grid(rows=NULL, vars(age))+
  labs(fill="sign(logFC)", x= "nmf pattern", y="# genes",
       title="Unique genes top-weighted to each NMF that are also age DEGs")+
  scale_fill_manual(values=c("#67a9cf","#ef8a62"))+theme_bw()

#heatmap

gene_list_unique = read.csv('gc_nmf_gene-weights.csv')

ggplot (nmf, aes(x=factor(nmf, level=c('nmf26', 'nmf5', 'nmf14')), y=factor(c(gene_list)))) + 
  geom_tile()



col_fun2 = colorRamp2(c(0, 0.002), c("white", "black"))

