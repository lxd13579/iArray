#!/bin/bash
# Date: 2021/7/30
# Update: 2021/7/30
# Author: Xiao dong Li
# Version: V1.0

# stop the script when an error occurs
set -ue

#**************************************************
## 1 Make project ----------------------------
#**************************************************

mkdir project

mkdir -p project/reference
mkdir -p project/data
mkdir -p project/scripts
mkdir -p project/results
mkdir -p project/tmp
mkdir -p project/logs
mkdir -p project/test
mkdir -p project/software

tree project

cd project

Species="Arabidopsis thaliana" # NCBI Taxonomy name
Reference=Tair10.chr.fasta   # reference genome sequence, data located in reference directory
RmFASTA=Tair10.repeats.sm.fa # soft-masked genome sequence, data located in reference directory
RmGFF=Tair10.repeats.gff3    # repeats annotation file, data located in reference directory
GffFile=Tair10.gene.gff3     # gene annotation file, data located in reference directory
OutputName="Tair10"          # output file name prefix, which was defined by users
SnpEffDir=project/software/snpEff

# All software can be downloaded from anaconda

#**************************************************
## 2 Data preprocess ----------------------------
#**************************************************

echo "`date '+[iArray] %D %T'` -- Data preprocess start..."

## step1: Mask repeats  ------------------------

mkdir -p reference/RepeatMasker

# Mask a fasta file based on repeats coordinates with bedtools

if [ -f reference/$RmFASTA ]
then
	cp reference/$RmFASTA reference/RepeatMasker/$Reference.masked
elif [ -f reference/$RmGFF ]
then
	bedtools maskfasta \
	-fi reference/$Reference \
	-bed reference/$RmGFF \
	-fo reference/RepeatMasker/$Reference.masked \
	-soft	
else
	RepeatMasker \
	-species $Species \
	-pa 10 \
	-a -inv -xsmall -excln -gff -html \
	-dir reference/RepeatMasker \
	reference/$Reference \
	> logs/$OutputName.RepeatMasker.log	
fi

if [ -s reference/RepeatMasker/$Reference.masked ]
then
	echo "`date '+[iArray] %D %T'` -- The soft-masked genome sequence has been generated successfully."
else
	echo "`date '+[iArray] %D %T'` -- The soft-masked genome sequence does not exist."
	exit
fi

## step2: Identify SNPs ------------------------

mkdir -p results/01.snp_discovery

OutputDir=results/01.snp_discovery

# build index for reference genome
samtools faidx reference/$Reference

samtools dict reference/$Reference -o $(basename -s .fasta reference/$Reference).dict

Mapcaller index reference/$Reference

# SNPs calling
MapCaller -i reference/$Reference -f data/sample1.clean.R1.fq.gz -f2 data/sample1.clean.R2.fq.gz -vcf $OutputDir/sample1.vcf -dup 1 -monomorphic -log logs/sample1.log

bgzip -@ 8 -c $OutputDir/sample1.vcf > $OutputDir/sample1.vcf.gz

bcftools index -f -t $OutputDir/sample1.vcf.gz --threads 32

# Merge multiple VCF files

bcftools merge $OutputDir/sample1.vcf.gz $OutputDir/sample2.vcf.gz -o $OutputDir/Raw.merge.vcf.gz -O z --threads 32
bcftools index -f -t $OutputDir/Raw.merge.vcf.gz --threads 32

## step3: Filter SNPs ------------------------

mkdir -p results/02.snp_filtering

InputDir=results/01.snp_discovery
OutputDir=results/02.snp_filtering
InputFile=$InputDir/Raw.merge.vcf.gz

#samtools faidx reference/$Reference

# some case should be noticed
# 1. The header of VCF file may be not contain chromosomes/contigs information
# 2. The CHROM column of VCF file may be a number rather than character, like 1, 2, 3, .. not chr01, chr02, chr03, .., chrN.
# method: bcftools annotate --rename-chrs chr_map.txt xxx.vcf.gz -o xxx.vcf.gz -O z --threads 10
# 3. The FILTER column of VCF file may be not contain any information, like PASS, q30 field.
# 4. The QUAL column of VCF file may be not contain any information, like 30, 40, and so on.
# 5. A VCF file at least should contain the followed information: CHROM, POS, REF, FIRST_ALT and GT columns.
# 6. If the genotyping data is plink format(.bed, .bim, and .fam), please convert into VCF format using PLINK software.
# 7. If the genotyping data is consisted of multiple files (one chromosome one file), please concat these files using Bcftools software before filtering.
bcftools reheader \
-f reference/$Reference.fai $InputDir/$InputFile | \
bcftools view \
-m2 -M2 -v snps -f PASS - \
-o $OutputDir/$OutputName.snp.rh.bcf \
-O u --threads 10

# MAF and missing filter
bcftools filter \
-i 'F_MISSING<=0.5 & MAF>=0.05 & QUAL>=30' $OutputDir/$OutputName.snp.rh.bcf | \
bcftools annotate \
--set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' - \
-o $OutputDir/$OutputName.snp.maf0.05.mis0.5.vcf.gz -O z \
--threads 10

echo "`date '+[iArray] %D %T'` -- Done."

#***************************************************
## 3 Probes evaluation ----------------------------
#***************************************************

echo "`date '+[iArray] %D %T'` -- Probes evaluation start..."

mkdir -p results/03.snp_evaluate

InputDir=results/02.snp_filtering
OutputDir=results/03.snp_evaluate

## step1: Add 100 bp up and downstream of each SNP
cut -f1,2 reference/$Reference.fai > reference/$OutputName.genome

bcftools query \
-f '%CHROM\t%POS0\t%END\t%ID\n' \
$InputDir/$OutputName.snp.maf0.05.mis0.5.vcf.gz \
> $OutputDir/$OutputName.snp.bed

bedtools slop \
-i $OutputDir/$OutputName.snp.bed \
-g reference/$OutputName.genome -b 100 \
> $OutputDir/$OutputName.snp.100bp.flank.bed

## step2: Extract probes sequence from reference genome and evaluation 
# (GC content: 0.35-0.65, no N base, no low complex/simple repeats)
bedtools getfasta \
-fi reference/RepeatMasker/$Reference.masked \
-bed $OutputDir/$OutputName.snp.100bp.flank.bed \
-fo $OutputDir/$OutputName.probes.fa -nameOnly

seqkit fx2tab \
-B n -B N -B atgc -B ATGC -I \
-g -l -n $OutputDir/$OutputName.probes.fa \
> $OutputDir/$OutputName.probes.stats.txt

awk '{if ($3>=35 && $3 <=65 && $7 == 100) print}' $OutputDir/$OutputName.probes.stats.txt | \
cut -f1 > $OutputDir/$OutputName.probes.gc.txt

seqkit grep \
-f $OutputDir/$OutputName.probes.gc.txt $OutputDir/$OutputName.probes.fa \
> $OutputDir/$OutputName.probes.gc.fa

## step3: Blast against the genome 
# (qcovs >= 50, alignment length >= 100, mismatch number <= 30, match hits <= 2)
mkdir -p $OutputDir/blastdb

makeblastdb \
-in reference/$Reference \
-dbtype nucl \
-parse_seqids \
-out $OutputDir/blastdb/$OutputName.genome

blastn \
-query $OutputDir/$OutputName.probes.gc.fa \
-out $OutputDir/$OutputName.probes.blast.txt \
-db $OutputDir/blastdb/$OutputName.genome \
-outfmt "6 qseqid sseqid evalue bitscore qcovs pident qlen length mismatch gapopen qstart qend sstart send" \
-word_size 10 \
-evalue 1e-5 \
-max_target_seqs 100 \
-num_threads 8

python3 ./scripts/snp_specific.py $OutputDir/$OutputName.probes.blast.txt $OutputDir/$OutputName.probes.specific.txt 50 100 0 30 1

cut -f1 $OutputDir/$OutputName.probes.specific.txt > $OutputDir/$OutputName.specific.id.txt

vcftools \
--gzvcf $InputDir/$OutputName.snp.maf0.05.mis0.5.vcf.gz \
--snps $OutputDir/$OutputName.specific.id.txt \
--out $OutputDir/$OutputName.specific --recode --stdout | bgzip > \
$OutputDir/$OutputName.specific.vcf.gz

## step4: GWAS markers evaluation

bedtools slop \
-i data/$OutputName.gwas.bed \
-g reference/$OutputName.genome -b 100 \
> $OutputDir/$OutputName.gwas.100bp.flank.bed

bedtools getfasta \
-fi reference/RepeatMasker/$Reference.masked \
-bed $OutputDir/$OutputName.gwas.100bp.flank.bed \
-fo $OutputDir/$OutputName.gwas.probes.fa -nameOnly

seqkit fx2tab \
-B n -B N -B atgc -B ATGC -I \
-g -l -n $OutputDir/$OutputName.gwas.probes.fa \
> $OutputDir/$OutputName.gwas.probes.stats.txt

awk '{if ($3>=35 && $3 <=65 && $7 == 100) print}' $OutputDir/$OutputName.gwas.probes.stats.txt | \
cut -f1 > $OutputDir/$OutputName.gwas.probes.gc.txt

seqkit grep \
-f $OutputDir/$OutputName.gwas.probes.gc.txt $OutputDir/$OutputName.gwas.probes.fa \
> $OutputDir/$OutputName.gwas.probes.gc.fa

blastn \
-query $OutputDir/$OutputName.gwas.probes.gc.fa \
-out $OutputDir/$OutputName.gwas.probes.blast.txt \
-db $OutputDir/blastdb/$OutputName.genome \
-outfmt "6 qseqid sseqid evalue bitscore qcovs pident qlen length mismatch gapopen qstart qend sstart send" \
-word_size 10 \
-evalue 1e-5 \
-max_target_seqs 100 \
-num_threads 8

python3 ./scripts/snp_specific.py $OutputDir/$OutputName.gwas.probes.blast.txt $OutputDir/$OutputName.gwas.probes.specific.txt 50 100 0 30 1

cut -f1 $OutputDir/$OutputName.gwas.probes.specific.txt > $OutputDir/$OutputName.gwas.specific.id.txt

echo "`date '+[iArray] %D %T'` -- Done."

#***************************************************
## 4 Effect prediction ----------------------------
#***************************************************

echo "`date '+[iArray] %D %T'` -- Effect prediction start..."

mkdir -p results/04.snp_annotate

# input and output data directory
InputDir=results/03.snp_evaluate
OutputDir=results/04.snp_annotate

## step1: build database and annotation

mkdir -p $SnpEffDir/data/$OutputName

cp reference/$Reference $SnpEffDir/data/$OutputName/sequences.fa

# chromosome code may be transformed into number code to consist with genome, like 1, 2, 3, 4, ...
#paste <(cut -f1 reference/$GffFile | sort | uniq) <(cut -f1 reference/$OutputName.genome | sort) > reference/chr2num.txt
#gffread reference/$GffFile -T -m reference/chr2num.txt -o $SnpEffDir/data/$OutputName/genes.gtf

grep -v "#" reference/$GffFile | gffread - -T -o $SnpEffDir/data/$OutputName/genes.gtf

echo "$OutputName.genome: $OutputName" >> $SnpEffDir/snpEff.config

java -Xmx40G -jar $SnpEffDir/snpEff.jar build -gtf22 -v $OutputName

# html and csv format + ann vcf
java -Xmx40G -jar $SnpEffDir/snpEff.jar $OutputName -v \
-stats $OutputDir/$OutputName.ann.html \
-csvStats $OutputDir/$OutputName.ann.csv \
$InputDir/$OutputName.specific.vcf.gz | bgzip > $OutputDir/$OutputName.ann.vcf.gz

# level1, level2 and level3
zcat $OutputDir/$OutputName.ann.vcf.gz | \
java -jar $SnpEffDir/SnpSift.jar filter "((ANN[0].IMPACT = 'HIGH') | (ANN[0].IMPACT = 'MODERATE'))" | \
java -jar $SnpEffDir/SnpSift.jar extractFields \
-s "," -e "." - CHROM POS REF ALT "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_P" \
> $OutputDir/$OutputName.AnnLevel1.txt

zcat $OutputDir/$OutputName.ann.vcf.gz | \
java -jar $SnpEffDir/SnpSift.jar filter "(ANN[0].IMPACT = 'LOW')" | \
java -jar $SnpEffDir/SnpSift.jar extractFields \
-s "," -e "." - CHROM POS REF ALT "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_P" \
> $OutputDir/$OutputName.AnnLevel2.txt

zcat $OutputDir/$OutputName.ann.vcf.gz | \
java -jar $SnpEffDir/SnpSift.jar filter "(ANN[0].IMPACT = 'MODIFIER')" | \
java -jar $SnpEffDir/SnpSift.jar extractFields \
-s "," -e "." - CHROM POS REF ALT "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_P" \
> $OutputDir/$OutputName.AnnLevel3.txt

echo "`date '+[iArray] %D %T'` -- Done."

#***************************************************
## 5 Marker selection ----------------------------
#***************************************************

echo "`date '+[iArray] %D %T'` -- Marker selection start..."

mkdir -p results/05.snp_select

# input and output data directory
InputDir=results/04.snp_annotate
OutputDir=results/05.snp_select

## step1: generate bed files
sed '1d' $InputDir/$OutputName.AnnLevel1.txt | \
awk '{print $1, $2-1, $2, "Level1"}' OFS="\t" \
> $OutputDir/$OutputName.AnnLevel1.bed

sed '1d' $InputDir/$OutputName.AnnLevel2.txt | \
awk '{print $1, $2-1, $2, "Level2"}' OFS="\t" \
> $OutputDir/$OutputName.AnnLevel2.bed

sed '1d' $InputDir/$OutputName.AnnLevel3.txt | \
awk '{print $1, $2-1, $2, "Level3"}' OFS="\t" \
> $OutputDir/$OutputName.AnnLevel3.bed


## step2: run panel_custom.r script for snp selection
# input: chromo number, the number of SNPs, bin size, level1.bed, level2.bed, level3.bed, ..., leveln.bed
chr_num=5 
window_size=100000
snp_number=50000

head -n $chr_num reference/$OutputName.genome > $OutputDir/$OutputName.chr.genome

for i in {1..3};do
bedtools makewindows -g $OutputDir/$OutputName.chr.genome -w $window_size | \
bedtools window -w 0 -a - -b $OutputDir/$OutputName.AnnLevel$i.bed \
> $OutputDir/$OutputName.AnnLevel$i.win.bed
done

bedtools makewindows -g $OutputDir/$OutputName.chr.genome -w $window_size | \
bedtools intersect \
-nonamecheck \
-a - \
-b $OutputDir/$OutputName.AnnLevel1.bed $OutputDir/$OutputName.AnnLevel2.bed $OutputDir/$OutputName.AnnLevel3.bed \
-C -names Level1 Level2 Level3 > $OutputDir/$OutputName.AnnLevel.count.bed

Rscript ./scripts/panel_custom.r \
--bed $OutputDir/$OutputName.AnnLevel.count.bed \
--indir $OutputDir \
--outdir $OutputDir/50k \
--prefix $OutputName.50k \
--number $snp_number

## step3: generate VCF file and annotation results from all selected markers
grep -h -w -v "Chr" $OutputDir/50k/*SnpSel.tsv | sort -k1,1 -k2,2n > $OutputDir/50k/$OutputName.50k.AllLevel.SnpSel.tsv

vcftools \
--gzvcf $InputDir/$OutputName.ann.vcf.gz \
--positions $OutputDir/50k/$OutputName.50k.AllLevel.SnpSel.tsv \
--out $OutputDir/50k/$OutputName.array50k --recode --recode-INFO-all --stdout | bgzip > \
$OutputDir/50k/$OutputName.array50k.vcf.gz

zcat $OutputDir/50k/$OutputName.array50k.vcf.gz | \
java -jar $SnpEffDir/SnpSift.jar extractFields \
-s "," -e "." - CHROM POS REF ALT "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_P" \
> $OutputDir/50k/$OutputName.array50k.ann.txt

vcftools \
--gzvcf $OutputDir/50k/$OutputName.array50k.vcf.gz \
--freq2 \
--out $OutputDir/50k/$OutputName.array50k.stats

echo "`date '+[iArray] %D %T'` -- Done."

#***************************************************
## 6 Panel summary ----------------------------
#***************************************************

echo "`date '+[iArray] %D %T'` -- Panel summary start..."

mkdir -p results/06.panel_summary

## step1: the coverage rate of bins --------------------
cat <(sed 's/_/\t/' results/03.snp_evaluate/$OutputName.gwas.specific.id.txt) results/05.snp_select/50k/$OutputName.50k.AllLevel.SnpSel.tsv | \
sort -k1,1 -k2,2n | \
uniq | \
awk '{print $1,$2-1,$2,$1"_"$2}' OFS="\t" > results/06.panel_summary/$OutputName.50k.ReqAndGWAS.bed

# final panel size
wc -l results/06.panel_summary/$OutputName.50k.ReqAndGWAS.bed

# coverage rate of bins
bedtools makewindows -g results/05.snp_select/$OutputName.chr.genome -w 100000 | \
bedtools window -w 0 -a - -b results/06.panel_summary/$OutputName.50k.ReqAndGWAS.bed -c \
> results/06.panel_summary/$OutputName.50k.snpCount.bed

awk '
BEGIN{
	a=0;b=0;
	printf "Empty bins\tNonempty bins\tTotal bins\tCoverage rate\n"
}{
	if ($4 == 0) {a++} else {b++}
}
END{
	printf "%s\t%s\t%s\t%s\n",a,b,a+b,b/(a+b)*100
}' results/06.panel_summary/$OutputName.50k.snpCount.bed

# step2: SNPs located in annotated genes --------------------

grep -v "#" reference/$GffFile | \
awk '$3=="gene"' | \
gffread - --gene2exon --bed | \
cut -f1-4 > results/06.panel_summary/$OutputName.genes.bed

bedtools window \
-w 0 \
-a results/06.panel_summary/$OutputName.genes.bed \
-b results/06.panel_summary/$OutputName.50k.ReqAndGWAS.bed \
-c > results/06.panel_summary/$OutputName.50k.geneCount.bed

awk '
BEGIN{
	printf "CHROM\tGenes\tSNPs\n"
}
{
	if ($5 > 0) {A[$1]=A[$1]+$5;B[$1]++}
}
END{
	totalA = 0;
	totalB = 0;
    asorti(A,S);
    # get length
    j = length(S);
    for (i = 1; i <= j; i++) {
		totalA += A[S[i]];
		totalB += B[S[i]];
        printf("%s\t%s\t%s\n", S[i],B[S[i]],A[S[i]])
    };
	printf("Total\t%s\t%s\n", totalB,totalA)
}' results/06.panel_summary/$OutputName.50k.geneCount.bed


# step3: gap and maf summary --------------------

Rscript ./scripts/panel_summary.r $OutputName

echo "`date '+[iArray] %D %T'` -- Done."

exit

##################### END #############################


