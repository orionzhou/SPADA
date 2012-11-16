f_eval = '/project/youngn/zhoup/Data/misc3/spfinder/Athaliana/41_eval/10_eval.tbl'
f_ref_gtb = '/project/youngn/zhoup/Data/misc3/spfinder/Athaliana/41_eval/01_model.gtb'
f_stat = '/project/youngn/zhoup/Data/misc3/spfinder/Athaliana/31_model/41_stat.tbl'
fo = '/project/youngn/zhoup/Data/misc3/spfinder/Athaliana/41_eval/21_stat.tbl'

tv1 = read.table(f_eval, sep="\t", header=T, as.is=T, na.strings=c(NA, ''))
t_stat = read.table(f_stat, sep="\t", header=T, as.is=T, na.strings=c(NA, ''))
genes = read.table(f_ref_gtb, sep="\t", header=T, as.is=T, na.strings=c(NA, ''))$id
genes_missed = genes[!genes %in% tv1$gene]

tv2 = rbind(tv1, data.frame(id=NA, gene=genes_missed, tag=10, lenO=0, len1=0, len2=100))
tv2$tag[(tv2$tag == 2 & tv2$len1+tv2$len2 > 30) | tv2$tag == 7 | tv2$tag == 8] = 5

tv3 = tv2[,c("id", "gene", "tag")]
tv4 = merge(tv3, t_stat, by="id", all.x=T)

write.table(tv4, file=fo, col.names=T, row.names=F, sep="\t", quote=F) 


