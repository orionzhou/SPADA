require(ggplot2)
require(RColorBrewer)

orgs = c('Athaliana', 'Mtruncatula_3.5', 'Osativa')
dirO = file.path(DIR_Data, "misc4/stat")

for (i in 1:length(orgs)) {
	org = orgs[i]
	fi = sprintf("%s/misc4/spada.crp.%s/41_perf_eval/51_stat.tbl", DIR_Data, org)
    ti = read.table(fi, sep="\t", header=T, as.is=T)
    ti = cbind(org=org, ti)
    if(i == 1) {
      t = ti
    } else {
      t = rbind(t, ti)
    }
}

#f05 = file.path(dirO, "05_stat.tbl")
#write.table(t, f05, sep="\t", quote=F, row.names=F, col.names=T)

te1 = t[t$soft == 'SPADA', c(-2, -4)]
te2 = reshape(te1, idvar=c("org", "e"), varying=list(3:6), timevar="type", v.names='value', times=colnames(te1)[3:6], direction='long')
te2 = cbind(te2, loge=log10(te2$e))
te2$type = factor(te2$type, levels=c('sp_nt','sp_exon','sn_nt','sn_exon'))
p <- ggplot(te2) + 
  geom_line(mapping=aes(loge, value, group=type), size=0.2) +
  geom_point(mapping=aes(loge, value, colour=type), size=1.5) +
  scale_colour_brewer(palette="Set1", labels=c("sn_nt"="Nucleotide Sensitivity", "sn_exon"='Exon Sensitivity', "sp_nt"='Nucleotide Specifity', "sp_exon"='Exon Specifity')) + 
  scale_x_continuous(name='E-value cutoff', breaks=seq(-8,0,1), labels=10^seq(-8,0,1)) +
  scale_y_continuous(name='Sensitivity / Specificity') +
  facet_grid(. ~ org) +
  theme(axis.text.x=element_text(size=7, hjust=1, vjust=1, angle=45), strip.text.x=element_text(face="italic")) +
  labs(shape='') +
  theme(legend.title=element_blank(), legend.position="top", legend.direction="horizontal", legend.text=element_text(size=7))
ggsave(file.path(dirO, "performance_e.pdf"), p, width=5, height=3.5)

ts1 = t[t$e == 0.001, c(-3, -4)]
ts2 = reshape(ts1, idvar=c("org", "soft"), varying=list(3:6), timevar="type", v.names='value', times=colnames(ts1)[3:6], direction='long')
ts2$type = factor(ts2$type, levels=c('sp_nt','sp_exon','sn_nt','sn_exon'))
ts2$soft = factor(ts2$soft, levels=c('GeneID', 'Augustus_de_novo', 'GlimmerHMM', 'GeneMark', 'GeneWise_SplicePredictor', 'Augustus_evidence', 'SPADA', 'All'))
p <- ggplot(ts2) + 
  geom_point(mapping=aes(soft, value, colour=type), size=1.5) +
#  geom_line(mapping=aes(soft, value, group=type), size=0.2) +
  scale_colour_brewer(palette="Set1", labels=c("sn_nt"="Nucleotide Sensitivity", "sn_exon"='Exon Sensitivity', "sp_nt"='Nucleotide Specifity', "sp_exon"='Exon Specifity')) + 
  scale_x_discrete(name='Organism') +
  scale_y_continuous(name='Sensitivity / Specificity') + 
  facet_grid(. ~ org) +
  theme(axis.text.x = element_text(size=7, hjust=1, vjust=1, angle=30), strip.text.x=element_text(face="italic")) +
  labs(fill='', colour='') +
  theme(legend.title=element_blank(), legend.position="top", legend.direction="horizontal", legend.text=element_text(size=7))
ggsave(file.path(dirO, "performance_soft.pdf"), p, width=5, height=3.5)

t1 = t[, -4]
t2 = reshape(t1, idvar=c("org", "soft", "e"), varying=list(4:7), timevar="type", v.names='value', times=colnames(t1)[4:7], direction='long')
t2 = cbind(t2, loge=log10(t2$e))
t2$type = factor(t2$type, levels=c('sp_nt','sp_exon','sn_nt','sn_exon'))
ts2$soft = factor(ts2$soft, levels=c('GeneID', 'Augustus_de_novo', 'GlimmerHMM', 'GeneMark', 'GeneWise_SplicePredictor', 'Augustus_evidence', 'SPADA', 'All'))
p <- ggplot(t2) + 
  geom_line(mapping=aes(loge, value, group=type), size=0.2) +
  geom_point(mapping=aes(loge, value, colour=type), size=1.5) +
  scale_colour_brewer(palette="Set1", labels=c("sn_nt"="Nucleotide Sensitivity", "sn_exon"='Exon Sensitivity', "sp_nt"='Nucleotide Specifity', "sp_exon"='Exon Specifity')) + 
  scale_x_continuous(name='E-value cutoff', breaks=seq(-8,0,1), labels=10^seq(-8,0,1)) +
  scale_y_continuous(name='Sensitivity / Specificity') +
  facet_grid(org ~ soft) +
  theme(axis.text.x = element_text(size=7, hjust=1, vjust=1, angle=45), strip.text.y=element_text(face="italic"), strip.text.x=element_text(size=7)) +
  theme(legend.title=element_blank(), legend.position="top", legend.direction="horizontal", legend.text=element_text(size=8)) +
  labs(fill='', colour='')
ggsave(file.path(dirO, "performance_all.pdf"), p, width=10, height=7)
