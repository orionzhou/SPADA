orgs = c("Athaliana", "Mtruncatula", "Osativa")
org = orgs[3] 

fi = file.path(DIR_Misc3, "spada", org, "41_eval/21_stat_SPADA_0.001.tbl")
t = read.table(fi, header=T, sep="\t", as.is=T)

table(t$tag)

ts = t[t$tag==9,1:9]
ts[order(ts$e),]

get_crp_cat = function(x) {if(x <= "CRP1030") "DEFL" else if(x <= "CRP1120") "CCP-like" else if(x <= "CRP1530") "CCP" else "other CRP"}
for (org in orgs) {
  fi = file.path(DIR_Misc3, "spada", org, "31_model_SPADA", "61_final.tbl")
  t = read.table(fi, header=T, sep="\t", as.is=T)
  t = cbind(t, group=sapply(t$family, get_crp_cat))
  table(t$group)
}


org = "Zmays_masked"
fi = file.path(DIR_Misc3, "spada", org, "31_model_SPADA/91_compare.tbl")
