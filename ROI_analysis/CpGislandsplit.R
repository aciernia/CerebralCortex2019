cpgs <- read.delim("/Users/annieciernia/Desktop/hg38CpG.bed",header=F)

islands <- cpgs %>% filter(V4 == "hg38_cpg_islands")
shelves <- cpgs %>% filter(V4 == "hg38_cpg_shelves")
shores <- cpgs %>% filter(V4 == "hg38_cpg_shores")
inter <- cpgs %>% filter(V4 == "hg38_cpg_inter")

write.table(islands,"/Users/annieciernia/Desktop/hg38_cpg_islands.bed",col.names=F,row.names=F,quote=F, sep="\t")
write.table(shelves,"/Users/annieciernia/Desktop/hg38_cpg_shelves.bed",col.names=F,row.names=F,quote=F, sep="\t")
write.table(shores,"/Users/annieciernia/Desktop/hg38_cpg_shores.bed",col.names=F,row.names=F,quote=F, sep="\t")
write.table(inter,"/Users/annieciernia/Desktop/hg38_cpg_inter.bed",col.names=F,row.names=F,quote=F, sep="\t")
