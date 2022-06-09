#Script for phasing genetic data and calling identical-by-descent tracts in preparation for running MAPS
#Paths provided for Augrabies users

module load shapeit/2.r904
module load vcftools

export DATA=/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/data
export PHASED_OUT=/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/PhasedOut
export PHASED_NoRel=/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/PhasedOut_NoRelatives
export IBD_OUT=/share/hennlab/users/espless/PhasingAndIBD/hapibd_Pagani_Scheinfeldt/IBDOut
export MAP=/share/hennlab/reference/recombination_maps/genetic_map_b37
export MAP2=/share/hennlab/reference/recombination_maps/hapmap_GRCh37_plinkFormat_genetic_map
export REF=/share/hennlab/reference/1000G_Phase3_haps-sample-legend/1000GP_Phase3
export PROGS=/share/hennlab/progs
export UNRELATED=/share/hennlab/users/espless/MergedDataEastAfrica/Merge3/KingFiles_Pagani2015_scheinfeldt_geno0.05_mind0.1

#1: Phase with SHAPEIT2

for chr in `seq 22`

do

#Run SHAPEIT checks
shapeit -check \
    -B $DATA/Pagani2015_scheinfeldt_geno0.05_mind0.1_${chr} \
    -M $MAP/genetic_map_chr${chr}_combined_b37.txt \
    --input-ref $REF/1000GP_Phase3_chr${chr}.hap.gz $REF/1000GP_Phase3_chr${chr}.legend.gz $REF/1000GP_Phase3.sample \
    --output-log $PHASED_OUT/gwas.alignments.chr${chr}

shapeit -B $DATA/Pagani2015_scheinfeldt_geno0.05_mind0.1_${chr} \
    -M $MAP/genetic_map_chr${chr}_combined_b37.txt \
    --input-ref $REF/1000GP_Phase3_chr${chr}.hap.gz $REF/1000GP_Phase3_chr${chr}.legend.gz $REF/1000GP_Phase3.sample \
    --exclude-snp $PHASED_OUT/gwas.alignments.chr${chr}.snp.strand.exclude \
    --duohmm \
    -W 5 \
    -O $PHASED_OUT/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5 \
    -T 6 \
    --output-log $PHASED_OUT/chr${chr} \
    --force

#2: Convert to VCF

shapeit -convert --input-haps $PHASED_OUT/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5 --output-vcf $PHASED_OUT/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5.vcf

#3: Exclude 1st and 2nd degree relatives and any individuals not in the geographic region of interest

vcftools --vcf $PHASED_OUT/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5.vcf --keep $UNRELATED/kingunrelated_EAfr_vcftools.txt  --recode --recode-INFO-all --out $PHASED_NoRel/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5_noRel.vcf

#3: Call IBD segments with hap-ibd v1.0

nice java -Xmx10g -jar $PROGS/hap-ibd.jar gt=$PHASED_NoRel/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5_noRel.vcf.recode.vcf map=$MAP2/plink.chr${chr}.GRCh37.map out=$IBD_OUT/Pagani_Scheinfeldt_ibd_chr${chr}.out

gunzip  $IBD_OUT/Pagani_Scheinfeldt_ibd_chr${chr}.out.ibd.gz

#4 Repair IBD segments with merge-ibd-segments
cat $IBD_OUT/Pagani_Scheinfeldt_ibd_chr${chr}.out.ibd | java -jar $PROGS/merge-ibd-segments.16May19.ad5.jar $PHASED_NoRel/Pagani_Scheinfeldt_${chr}_SHAPEITphased_duoHMM_W5_noRel.vcf.recode.vcf  $MAP2/plink.chr${chr}.GRCh37.map 0.6 1 > $IBD_OUT/Pagani_Scheinfeldt_ibd_chr${chr}_merged.ibd

done
