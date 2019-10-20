library(ggplot2)
library(gridExtra)
library(grid)
library(gtools)
library(GenomicRanges)
# copy number maximum value to plot
ymax = 6

args = commandArgs(T)
samplename = args[1]
logr_file = args[2] # logr estimated profile
baf_file = args[3] # baf estimated profile
baf_purity_file = args[4] # rho_and_psi
outdir = args[5]


# samplename = "P9"
# outdir = ""
# 
# # logr_dir = "truth_rho_psi_all_v3/" # run on true purity and psi calculated from true CNA profile from all segments - v3 pipeline with adapted psi
# logr_dir = "battenberg_rho_psi_v3/" # run on battenberg guess of rho and psi - v3 pipeline
# 
# logr_file = file.path(logr_dir, paste0(samplename, "-2017_logr_total_cn.txt"))
# baf_file = paste0("baf//", samplename, "-2017/", samplename, "-2017_subclones.txt")
# baf_purity_file = paste0("baf//", samplename, "-2017/", samplename, "-2017_rho_and_psi.txt")


#' Calc total copy number per segment from a subclones data.frame
#' @noRd
calculate_bb_total_cn = function(bb) {
  return((bb$nMaj1_A+bb$nMin1_A)*bb$frac1_A + ifelse(!is.na(bb$frac2_A), (bb$nMaj2_A+bb$nMin2_A)*bb$frac2_A, 0))
}

calc_ploidy = function(startpos, endpos, total_cn) {
  len = endpos/1000-startpos/1000
  ploidy = sum(total_cn*len, na.rm=T) / sum(len, na.rm=T)
  return(ploidy)
}


logr = read.table(logr_file, header=T, stringsAsFactors=F)
# to prevent overlapping segments
logr$startpos = logr$startpos+1
baf = read.table(baf_file, header=T, stringsAsFactors=F)
purity = read.table(baf_purity_file, header=T, stringsAsFactors=F)["FRAC_GENOME", "rho"]

breakpoints_logr = data.frame(chrom=c(logr$chr,logr$chr), pos=c(logr$start, logr$end), stringsAsFactors=F)
breakpoints_logr$chrom = factor(breakpoints_logr$chrom, levels=gtools::mixedsort(unique(breakpoints_logr$chrom)))
breakpoints_baf = data.frame(chrom=c(baf$chr, baf$chr), pos=c(baf$startpos, baf$endpos), stringsAsFactors=F)
breakpoints_baf$chrom = factor(breakpoints_baf$chrom, levels=gtools::mixedsort(unique(breakpoints_baf$chrom)))

p1 = ggplot() +
  # geom_rect(data=truth_s, mapping=aes(xmin=startpos, xmax=endpos, ymin=total_cn-rect_height_padding, ymax=total_cn+rect_height_padding), fill="#E69F00") +
  geom_vline(data=breakpoints_baf, mapping=aes(xintercept=pos), linetype=1, colour="red") +
  geom_vline(data=breakpoints_logr, mapping=aes(xintercept=pos), linetype=2, colour="blue") +
  facet_grid(~chrom, scales="free_x") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x=element_blank(),
        panel.spacing = unit(1, "points"))


# synchronise the breakpoints
breakpoints_all = unique(rbind(breakpoints_baf, breakpoints_logr))
breakpoints_all = breakpoints_all[order(paste0(breakpoints_all$chrom, "_", breakpoints_all$pos)),]

combine_segmentation = data.frame()
for (chrom in unique(breakpoints_all$chrom)) {
  bp_chrom = breakpoints_all[breakpoints_all$chrom==chrom,]
  bp_chrom = bp_chrom[order(bp_chrom$pos, decreasing=F),]
  
  for (i in 1:(nrow(bp_chrom)-1)) {
    combine_segmentation = rbind(combine_segmentation, 
                                 data.frame(chrom=chrom,
                                            startpos=bp_chrom$pos[i],
                                            endpos=bp_chrom$pos[i+1],
                                            stringsAsFactors=F))
  }
}
# remove breakpoints that are 1bp away
combine_segmentation = combine_segmentation[(combine_segmentation$endpos-combine_segmentation$startpos) > 1,]

logr_gr = makeGRangesFromDataFrame(logr, start.field="startpos", end.field="endpos", keep.extra.columns=T)
baf_gr = makeGRangesFromDataFrame(baf, start.field="startpos", end.field="endpos", keep.extra.columns=T)
combine_segmentation_gr = makeGRangesFromDataFrame(combine_segmentation, start.field="startpos", end.field="endpos", keep.extra.columns=T)

overlap = findOverlaps(combine_segmentation_gr, logr_gr)
combine_segmentation$logr_total_cn = NA
combine_segmentation$logr_total_cn[queryHits(overlap)] = logr$total_cn[subjectHits(overlap)]

overlap = findOverlaps(combine_segmentation_gr, baf_gr)
combine_segmentation[,c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")] = NA
combine_segmentation[queryHits(overlap), c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")] = baf[subjectHits(overlap), c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]

combine_segmentation$nMaj = ifelse(combine_segmentation$frac1_A==1, combine_segmentation$nMaj1_A, combine_segmentation$nMaj1_A*combine_segmentation$frac1_A + combine_segmentation$nMaj2_A*combine_segmentation$frac2_A)
combine_segmentation$nMin = ifelse(combine_segmentation$frac1_A==1, combine_segmentation$nMin1_A, combine_segmentation$nMin1_A*combine_segmentation$frac1_A + combine_segmentation$nMin2_A*combine_segmentation$frac2_A)
combine_segmentation$nTot = combine_segmentation$nMaj + combine_segmentation$nMin

combine_segmentation$chrom = factor(combine_segmentation$chrom, levels=gtools::mixedsort(unique(combine_segmentation$chrom)))

rect_height_padding = 0.2
p = ggplot() +
  geom_rect(data=combine_segmentation, mapping=aes(xmin=startpos, xmax=endpos, ymin=nMin-rect_height_padding, ymax=nMin+rect_height_padding), fill="#2f4f4f") +
  geom_rect(data=combine_segmentation, mapping=aes(xmin=startpos, xmax=endpos, ymin=nTot-rect_height_padding, ymax=nTot+rect_height_padding), fill="#E69F00") +
  geom_rect(data=combine_segmentation, mapping=aes(xmin=startpos, xmax=endpos, ymin=logr_total_cn-rect_height_padding, ymax=logr_total_cn+rect_height_padding), fill="green") +
  facet_grid(~chrom, scales="free_x") +
  scale_y_continuous(breaks=0:ymax) +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(1, "points")) + 
  coord_cartesian(ylim=c(0, ymax))

# baf purity / ploidy / num segments
# logr ploidy and possibly purity / num segments
# proportion genome agree ?
combine_segmentation$baf_total_cn = calculate_bb_total_cn(combine_segmentation)
baf_ploidy = calc_ploidy(combine_segmentation$startpos, combine_segmentation$endpos, combine_segmentation$baf_total_cn)
logr_ploidy = calc_ploidy(combine_segmentation$startpos, combine_segmentation$endpos, combine_segmentation$logr_total_cn)


combine_segmentation$len = combine_segmentation$endpos/1000 - combine_segmentation$startpos/1000
percent_genome_agree = round((sum(combine_segmentation$len[combine_segmentation$baf_total_cn==combine_segmentation$logr_total_cn & combine_segmentation$chrom!="X" & combine_segmentation$chrom!="Y"], na.rm=T) / sum(combine_segmentation$len[combine_segmentation$chrom!="X" & combine_segmentation$chrom!="Y"], na.rm=T)) * 100)
percent_genome_agree_rounded = round((sum(combine_segmentation$len[round(combine_segmentation$baf_total_cn)==round(combine_segmentation$logr_total_cn) & combine_segmentation$chrom!="X" & combine_segmentation$chrom!="Y"], na.rm=T) / sum(combine_segmentation$len[combine_segmentation$chrom!="X" & combine_segmentation$chrom!="Y"], na.rm=T)) * 100)


if (any(combine_segmentation$logr_total_cn==1)) {
  single_copy_index = which(combine_segmentation$logr_total_cn==1 & combine_segmentation$chrom!="X" & combine_segmentation$chrom!="Y")
  selected_index = single_copy_index[which.max(combine_segmentation$len[single_copy_index])]
  nMaj = 1
  nMin = 0
  overlap = findOverlaps(combine_segmentation_gr[selected_index,], baf_gr)
  segment_baf = baf$BAF[subjectHits(overlap)[1]]
  logr_purity = round((2*segment_baf-1)/(2*segment_baf-segment_baf*(nMaj+nMin)-1+nMaj), 2)
} else {
  logr_purity = NA
}


plot_title = paste0(samplename, "\n BAF|logR purity=", round(purity, 2), "|", logr_purity,
                    " ploidy=", round(baf_ploidy, 2), "|", round(logr_ploidy, 2),
                    " segs=", nrow(baf), "|", nrow(logr),
                    " agree=", percent_genome_agree, "%",
                    " agree rounded=", percent_genome_agree_rounded, "%")

png(file.path(outdir, paste0(samplename, "_compare.png")), height=350, width=1200)
grid.arrange(p1, p, heights=c(1.5/5,3.5/5),
             top=plot_title)
dev.off()

output = data.frame(samplename=samplename,
                    baf_purity=purity,
                    baf_ploidy=baf_ploidy,
                    logr_purity=logr_purity,
                    logr_ploidy=logr_ploidy,
                    perc_agree=percent_genome_agree,
                    perc_agree_round=percent_genome_agree_rounded)
write.table(output, file=file.path(outdir, paste0(samplename, "_compare.txt")), quote=F, sep="\t", row.names=F)





# get_ccf = function(vcf, mt_output, purity) {
#   num_muts = nrow(mt_output)
#   output = data.frame(chromosome=as.character(seqnames(vcf)), position=as.numeric(start(vcf)), type=rep(NA, num_muts),
#                       ccf=rep(NA, num_muts), major_cn=mt_output$MajCN, minor_cn=mt_output$MinCN, mcn=rep(NA, num_muts), mult=mt_output$MutCN,
#                       chromosome2=rep(NA, num_muts), position2=rep(NA, num_muts),
#                       altcount=mt_output$altCount, wtcount=mt_output$wtCount, vaf=mt_output$altCount/(mt_output$altCount+mt_output$wtCount))
#   output$mcn = mutationBurdenToMutationCopyNumber(burden=mt_output$altCount / (mt_output$altCount + mt_output$wtCount), cellularity=purity, normalCopyNumber=rep(2, num_muts), totalCopyNumber=output$major_cn + output$minor_cn)
#   output$ccf = output$mcn / output$mult
#   return(output)
# }
# 
# res = get_ccf(vcf_snv, MCN$D, purity)
# res$yulia = (2+ purity* ((res$major_cn+res$minor_cn)-2))*res$vaf / purity
# 
# dat$mcn_cp = (2+ purity* ((dat$major_cn+dat$minor_cn)-2))*dat$vaf / purity
