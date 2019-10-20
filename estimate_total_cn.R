library(Battenberg)
# library(ggplot2)
# library(gridExtra)
library(GenomicRanges)
library(coin)

phasekmin = 25
phasegamma = 100
rect_height_padding = 0.2
num_bootstrap_iters = 100 # number of tests to do checking for subclonality
num_bootstrap_samples = 1000 # number of samples to base the test on

args = commandArgs(T)
samplename = args[1]
workdir = args[2]
purity = as.numeric(args[3])
psi_t = as.numeric(args[4])
sex = args[5]

# samplename = "NASCR-0016"
# samplename = "T9-2017"

setwd(workdir)

logr2tumcn = function(cellularity, total_ploidy, logR, normal_total_cn=2) {
  return(((total_ploidy*(2^logR)) - normal_total_cn*(1-cellularity)) / cellularity)
}

# load logr
logr = Battenberg::read_table_generic(paste0(samplename, "_mutantLogR.tab"))
colnames(logr)[3] = "raw_logr"
# logr$logr_smoothed = Battenberg:::runmed_data(logr$Chromosome, logr$raw_logr, 101)
logr$logr_smoothed = logr$raw_logr

# load purity/ploidy
# rho_psi = read.table(paste0(samplename, "_rho_and_psi.txt"), header=T, stringsAsFactors=F)
# purity = rho_psi["FRAC_GENOME", "rho"]
# psi = rho_psi["FRAC_GENOME", "psi"]
# # psi = Battenberg:::psit2psi(psi_t=psi, rho=purity)

# read subclones
subclones = readr::read_tsv(paste0(samplename, "_subclones.txt"))
subclones$total_major = Battenberg:::calc_total_cn_major(subclones)
subclones$total_minor = Battenberg:::calc_total_cn_minor(subclones)
subclones$total_cn = subclones$total_minor + subclones$total_major

#est_psi_t = Battenberg:::calc_ploidy(subclones)
est_psi = Battenberg:::psit2psi(purity, psi_t)

# calc total cn
logr$total_cn = NA
logr$total_cn_psi = NA
for (i in (1:nrow(subclones))) {
  print(i)
  sel = which(logr$Chromosome == subclones$chr[i] & logr$Position >= subclones$startpos[i] & logr$Position <= subclones$endpos[i])
  # setting the expected normal allele to 1 for X and Y
  if (sex=="male") {
    normal_total_cn = rep(ifelse(subclones$chr[i]=="X" | subclones$chr[i]=="Y", 1, 2), length(sel))
  } else {
   normal_total_cn = rep(2, length(sel))
  }
  #logr$total_cn[sel] = logr2tumcn(purity, est_psi, logr$logr_smoothed[sel], normal_total_cn=normal_total_cn)
  logr$total_cn[sel] = Battenberg:::logr2tumcn(purity, est_psi, logr$logr_smoothed[sel])
}
if (sex=="male") {
  logr$total_cn[logr$Chromosome=="X" | logr$Chromosome=="Y"] = logr$total_cn[logr$Chromosome=="X" | logr$Chromosome=="Y"] / 2
}
# load bafsegmented + determine breakpoints

# segment logr (in future with baf breakpoints?)
logroutput = data.frame()
for (chr in unique(logr$Chromosome)) {
  print(chr)
  
  logr_chr = as.data.frame(logr[logr$Chromosome==chr, ])
  items_not_na = !is.na(logr_chr$total_cn)
  # sdev <- Battenberg:::getMad(logr_chr$total_cn[items_not_na], k=25)
  # 
  # # Standard deviation is not defined for a single value
  # if (is.na(sdev)) {
  #   sdev = 0
  # }
  
  # #DCW 250314
  # #for cell lines, sdev goes to zero in regions of LOH, which causes problems.
  # #0.09 is around the value expected for a binomial distribution around 0.5 with depth 30
  # if(sdev<0.09){
  #   sdev = 0.09
  # }

  res = Battenberg:::selectFastPcf(logr_chr$total_cn[items_not_na], phasekmin, phasegamma,T)
  total_cn_segm = res$yhat
  # r = rle(total_cn_segm)
  # r$lengths

  logroutput = rbind(logroutput, data.frame(Chromosome=logr_chr$Chromosome[items_not_na],
                                            Position=logr_chr$Position[items_not_na],
                                            total_cn_segm=total_cn_segm))
}

# create a copy number profile
d = data.frame()
for (chrom in unique(logroutput$Chromosome)) {
  print(chrom)
  dat = subset(logroutput, logroutput$Chromosome==chrom)
  r = rle(dat$total_cn_segm)
  cs_index = cumsum(r$lengths)
  d = rbind(d, data.frame(chr=dat$Chromosome[cs_index],
               startpos=c(dat$Position[1], dat$Position[cs_index][1:(length(cs_index)-1)]),
               endpos=c(dat$Position[cs_index]),
               raw_total_cn=dat$total_cn_segm[cs_index]))
}

d$total_cn = NA
d$logr_mean = NA
d$logr_median = NA
d$logr_sd = NA
# now test segments for clonality
for(i in 1:nrow(d)) {
  print(i)
  
  # get raw_total_cn
  raw_total = d$raw_total_cn[i]
  raw_total_min = floor(raw_total)
  raw_total_max = ceiling(raw_total)
  
  # get total_cn of SNPs in this segment
  sel = logr$Chromosome==d$chr[i] & logr$Position >= d$startpos[i] & logr$Position <= d$endpos[i]
  num_snps_segment = sum(sel)
  snps_sd = sd(logr$total_cn[sel], na.rm=T)
  d$logr_median[i] = mean(logr$raw_logr[sel], na.rm=T)
  d$logr_mean[i] = median(logr$raw_logr[sel], na.rm=T)
  d$logr_sd[i] = sd(logr$raw_logr[sel], na.rm=T)
  if (length(num_snps_segment)==1) {
    snps_sd = 0.7
  }
  
  min_pvals = rep(NA, num_bootstrap_iters)
  max_pvals = rep(NA, num_bootstrap_iters)
  # perform t-test for min and max
  for (j in 1:num_bootstrap_iters) {
    obs_data = rnorm(n = num_bootstrap_samples, mean = raw_total, sd = snps_sd)
    # obs_data = sample(logr$total_cn_psi[sel], num_bootstrap_samples, replace = T)
    ref_data_min = rnorm(n = num_bootstrap_samples, mean = raw_total_min, sd = snps_sd)
    ref_data_max = rnorm(n = num_bootstrap_samples, mean = raw_total_max, sd = snps_sd)
  
    # min_pvals[j] = t.test(obs_data, ref_data_min)$p.value
    # max_pvals[j] = t.test(obs_data, ref_data_max)$p.value
    
    # min_pvals[j] = ks.test(obs_data, ref_data_min)$p.value
    # max_pvals[j] = ks.test(obs_data, ref_data_max)$p.value
    
    DV <- c(obs_data, ref_data_min)
    IV <- factor(rep(c("A", "B"), c(num_bootstrap_samples, num_bootstrap_samples)))
    min_pvals[j] = pvalue(oneway_test(DV ~ IV, alternative="greater", distribution=approximate(B=9999)))
    DV <- c(obs_data, ref_data_max)
    max_pvals[j] = pvalue(oneway_test(DV ~ IV, alternative="less", distribution=approximate(B=9999)))
  }
    
  # obs_data = rnorm(n = num_snps_segment, mean = raw_total, sd = snps_sd)
  # ref_data_min = rnorm(n = num_snps_segment, mean = raw_total_min, sd = snps_sd)
  # ref_data_max = rnorm(n = num_snps_segment, mean = raw_total_max, sd = snps_sd)
  # 
  # DV <- c(obs_data, ref_data_min)
  # IV <- factor(rep(c("A", "B"), c(num_snps_segment, num_snps_segment)))
  # min_pval = pvalue(oneway_test(DV ~ IV, alternative="greater", distribution=approximate(B=9999)))
  # DV <- c(obs_data, ref_data_max)
  # max_pval = pvalue(oneway_test(DV ~ IV, alternative="less", distribution=approximate(B=9999)))
  
  min_signif = all(min_pvals < (0.05/num_bootstrap_iters))
  max_signif = all(max_pvals < (0.05/num_bootstrap_iters))
  
  # case 1: minimum not significantly different -> round down to obtain clonal copy number
  if (!min_signif & max_signif) {
    d$total_cn[i] = raw_total_min
  # case 2: maximum not significantly different -> round up to obtain clonal copy number  
  } else if (min_signif & !max_signif) {
    d$total_cn[i] = raw_total_max
  # case 3: both significantly different -> keep current estimate as subclonal copy number
  } else if (min_signif & max_signif) {
    d$total_cn[i] = raw_total
  # case 4: both not significant -> we cannot distinguish between the states
  #         round conservatively to the closest clonal state
  } else {
    d$total_cn[i] = round(raw_total)
  }
}

d$chr = factor(d$chr, levels=gtools::mixedsort(unique(d$chr)))
subclones$chr = factor(subclones$chr, levels=gtools::mixedsort(unique(subclones$chr)))

# p1 = ggplot() +
#   geom_rect(data=d, mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E69F00") + 
#   facet_grid(~chr, scales="free_x") +
#   ylim(1,10)
# 
# p2 = ggplot() +
#   geom_rect(data=subclones, mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E69F00") + 
#   facet_grid(~chr, scales="free_x") +
#   ylim(1,10)
# 
# png("test.png", width=1200, height=400)
# grid.arrange(p1, p2, ncol=1)
# dev.off()

d.gr = makeGRangesFromDataFrame(d, start.field="startpos", end.field="endpos", keep.extra.columns=T)
subclones.gr = makeGRangesFromDataFrame(subclones, start.field="startpos", end.field="endpos", keep.extra.columns=T)
overlap = findOverlaps(d.gr, subclones.gr)
multi_hits = queryHits(overlap)[duplicated(queryHits(overlap))]
is_duplicated = duplicated(queryHits(overlap))

d$baf_total_cn = unlist(lapply(1:length(overlap), function(i) { 
  x = queryHits(overlap)[i]
  if (is_duplicated[i]) {
    NULL 
  } else if (x %in% multi_hits & !is_duplicated[i]) {
    NA
  } else {
    subclones$total_cn[subjectHits(overlap)[queryHits(overlap)==x]]
  }
}))

len = d$endpos/1000 - d$startpos/1000
d$diff_with_baf_estimate = abs(d$total_cn-d$baf_total_cn)*len / sum(len, na.rm=T)
d$rho = purity
d$psi = est_psi
d$psi_t = psi_t
write.table(d, file=paste0(samplename, "_logr_total_cn.txt"), quote=F, sep="\t", row.names=F)
save.image(file=paste0(samplename, "_logr_total_cn.RData"))

# 
# diff = sum(abs(d$total_cn-d$baf_total_cn)*len, na.rm=T) / sum(len, na.rm=T)
# 
# p = ggplot(d[!is.na(d$baf_total_cn),]) + aes(x=total_cn, y=baf_total_cn) + geom_abline(intercept=0, slope=1, linetype=2) + geom_point()
# png("compare.png", height=700, width=700)
# print(p)
# dev.off()

