#Rscript to convert start and end from genomic positions to centimorgans
#Code adapted from interpolate.R (https://github.com/halasadi/ibd_data_pipeline/blob/master/interpolate.R) by Hussein Al-Asadi

for(chr in 1:22) {

mapdir <- "/share/hennlab/reference/recombination_maps/hapmap_GRCh37_plinkFormat_genetic_map/plink."
map <- read.table(paste0(mapdir, "chr", chr, ".GRCh37.map"), as.is=TRUE, header= F, stringsAsFactors = FALSE);
names(map) <- c('chr', 'dot', 'Genetic_Map.cM.', 'position')
position = map$position
genetic_map_cm = map$Genetic_Map.cM.

ibd.merged <- read.table(paste0("/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/IBDOut/Pagani_Scheinfeldt_ibd_chr", chr, "_merged.ibd"))
#ibd <- data.frame(read.table(ibd.merged, header = F, stringsAsFactors = F));
names(ibd.merged) <- c("id1", "hap1", "id2", "hap2", "chr", "start", "end", "lod", "lod2");

start <- approx(position, genetic_map_cm, xout = ibd.merged$start,
                yleft = min(genetic_map_cm), yright = max(genetic_map_cm))

end <- approx(position, genetic_map_cm, xout = ibd.merged$end,
                yleft = min(genetic_map_cm), yright = max(genetic_map_cm))

ibd.interp <- data.frame(id1=ibd.merged$id1, hap1 = ibd.merged$hap1, id2=ibd.merged$id2, hap2 = ibd.merged$hap2, chr = ibd.merged$chr,
                       start=start$y, end=end$y, lod = ibd.merged$lod, stringsAsFactors = F)

write.table(ibd.interp, paste0("/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/IBDOut/Pagani_Scheinfeldt_ibd_chr", chr, "_merged_cM.ibd"), row.names=FALSE,sep="\t", quote = FALSE)

}
