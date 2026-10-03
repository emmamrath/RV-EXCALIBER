#!/bin/bash
#=======================================================================================================================
# build_gnomad211_amrevel_control.sh
#
# Rebuilds the RV-EXCALIBER gnomAD 2.1.1 (hg19) control dataset from ANNOVAR hg19_gnomad211_exome, annotated with the
# same table_annovar command used for internal datasets:
#   -protocol refGeneWithVer,gnomad211_exome,dbnsfp47a -operation g,f,f -arg -exonicsplicing,,
# Missense variants are flagged in column 'pathogenic_missense':
#   1 if AlphaMissense_pred == "P" or REVEL_score >= 0.644, otherwise 0.
# Only sites that are PASS in gnomAD 2.1.1 (gnomAD_211_exomes_hg19_filter) are kept.
#
# Output: ${outdir}/${chr}_gnomAD_pruned_annotation_for_R_input_{rcc,hcc}.txt.gz
#
# Run in two stages:
#   qsub -v MODE=prep build_gnomad211_amrevel_control.sh                # once
#   qsub -J 1-22 -v MODE=annotate build_gnomad211_amrevel_control.sh    # per chromosome, after prep
#=======================================================================================================================
#PBS -N build_gnomad211_amrevel
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -l walltime=24:00:00
#PBS -j oe

set -euo pipefail

#module load bedtools
export PATH=/srv/scratch/z3531501/software/bedtools/bedtools2/bin:${PATH}
module load r/4.5.1

#-----------------------------------------------------------------------------------------------------------------------
# Paths (edit if needed)
#-----------------------------------------------------------------------------------------------------------------------

annovar=/srv/scratch/z3531501/software/annovar/annovar
humandb=${annovar}/humandb
db=${humandb}/hg19_gnomad211_exome.txt
refgene=${humandb}/hg19_refGeneWithVer.txt
ref_fa=/srv/scratch/z3531501/references_and_databases/hs37d5/hs37d5.fa          # uncompressed or bgzip + .fai/.gzi
rvx=/srv/scratch/z3531501/software/RV-EXCALIBER-fork
filter_file=${rvx}/gnomAD_211_exomes_hg19_filter/ALL_CHROM_gnomAD_filter.txt.gz
hcc_bed=${rvx}/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed
outdir=${rvx}/gnomAD_211_exomes_hg19_AMREVEL
work=${outdir}/work
formatter=${rvx}/scripts/Rscripts/format_gnomad211_amrevel_pruned.R

MODE=${MODE:-prep}
mkdir -p ${work}/vcf ${work}/annovar

#-----------------------------------------------------------------------------------------------------------------------
# Stage 1: prep (run once)
#-----------------------------------------------------------------------------------------------------------------------

if [[ ${MODE} == "prep" ]]; then

    # 1a. coding exons (CDS) +/- 10 bp from refGeneWithVer, autosomes, no "chr" prefix
    awk 'BEGIN {OFS="\t"} $7 < $8 {
           split($10, s, ","); split($11, e, ",")
           for (i = 1; i <= $9; i++) {
             st = (s[i] > $7 ? s[i] : $7); en = (e[i] < $8 ? e[i] : $8)
             c = $3; sub(/^chr/, "", c)
             if (st < en && c ~ /^[0-9]+$/) print c, (st - 10 < 0 ? 0 : st - 10), en + 10
           }}' ${refgene} \
      | sort -k1,1 -k2,2n | bedtools merge > ${work}/coding_pad10.bed

    # 1b. gnomAD 2.1.1 sites (ANNOVAR coordinates/alleles) within coding regions
    awk -F'\t' 'NR > 1 && $1 ~ /^[0-9]+$/ {print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' ${db} \
      | bedtools intersect -a - -b ${work}/coding_pad10.bed -u > ${work}/sites.tsv
    echo "coding sites: $(wc -l < ${work}/sites.tsv)"

    # 1c. anchor bases for indels: insertion (Ref "-") = base at Start; deletion (Alt "-") = base at Start - 1
    awk 'BEGIN {OFS="\t"} $4 == "-" {print $1, $2, $2 + 1} $5 == "-" {print $1, $2 - 1, $2}' ${work}/sites.tsv \
      | sort -u -k1,1 -k2,2n > ${work}/anchors.bed
    bedtools getfasta -fi ${ref_fa} -bed ${work}/anchors.bed -tab > ${work}/anchors.tab

    # 1d. PASS site IDs (chr:pos:ref:alt, VCF-style) from the gnomAD 2.1.1 filter file
    zcat ${filter_file} | awk '$2 == "PASS" {print $1}' > ${work}/pass_ids.txt
    echo "PASS sites in filter file: $(wc -l < ${work}/pass_ids.txt)"

    # 1e. one VCF per chromosome, PASS sites only (ID = chr:pos:ref:alt; dummy sample so -vcfinput works)
    rm -f ${work}/vcf/chr*.body
    awk -F'\t' -v w=${work}/vcf 'BEGIN {OFS="\t"}
      FILENAME == ARGV[1] {pass[$1] = 1; next}
      FILENAME == ARGV[2] {a[$1] = toupper($2); next}
      {
        c = $1; s0 = $2 + 0; start = s0 + 1; ref = $4; alt = $5
        ins_key = c ":" s0 "-" (s0 + 1)
        del_key = c ":" (s0 - 1) "-" s0
        if (ref == "-")      {pos = start;     anc = a[ins_key]; vref = anc;     valt = anc alt}
        else if (alt == "-") {pos = start - 1; anc = a[del_key]; vref = anc ref; valt = anc}
        else                 {pos = start;     vref = ref;       valt = alt}
        id = c ":" pos ":" vref ":" valt
        write_ok = (vref != "" && valt != "" && vref != valt && (id in pass))
        if (write_ok) print c, pos, id, vref, valt, ".", "PASS", ".", "GT", "0/1" > (w "/chr" c ".body")
      }' ${work}/pass_ids.txt ${work}/anchors.tab ${work}/sites.tsv

    for c in $(seq 1 22); do
        {
          printf '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDUMMY\n'
          sort -k2,2n -k3,3 -u ${work}/vcf/chr${c}.body
        } > ${work}/vcf/chr${c}.vcf
        rm ${work}/vcf/chr${c}.body
        echo "chr${c}: $(grep -vc '^#' ${work}/vcf/chr${c}.vcf) PASS coding sites"
    done
fi

#-----------------------------------------------------------------------------------------------------------------------
# Stage 2: annotate + format (PBS array, one chromosome per sub-job)
#-----------------------------------------------------------------------------------------------------------------------

if [[ ${MODE} == "annotate" ]]; then

    c=${CHR:-${PBS_ARRAY_INDEX:-}}
    if [[ -z "${c}" ]]; then echo "Set CHR (qsub -v MODE=annotate,CHR=N) or submit as an array (-J)" >&2; exit 1; fi
    infile=${work}/vcf/chr${c}.vcf
    outfile=${work}/annovar/chr${c}

    # identical annotation command to the internal datasets
    ${annovar}/table_annovar.pl "${infile}" ${humandb}/ -vcfinput -buildver hg19 \
        -out "${outfile}" -remove \
        -protocol refGeneWithVer,gnomad211_exome,dbnsfp47a \
        -operation g,f,f -nastring . \
        -arg -exonicsplicing,,

    Rscript ${formatter} ${outfile}.hg19_multianno.txt ${c} ${outdir} ${hcc_bed}
fi
