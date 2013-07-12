#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
GATK=/opt/GATK/2.4-3-g2a7af43/GenomeAnalysisTK.jar
JAVA=/opt/java/jdk1.6.0_38/bin/java
genome=/rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa
dbSNP=/rhome/cjinfeng/HEG4_cjinfeng/Variations/dbSNP_rice/dbSNP_HEG4/HEG4_dbSNP.vcf
repeat=/rhome/cjinfeng/HEG4_cjinfeng/GenomeAlign/Lastz/input/mask/MSU_r7.fa.RepeatMasker.out.bed
BAM=bam.list
prefix=ALL
mem=1
cpu=6

#echo "ReduceReads for each bam using nway co-reduce, which could keep these reads that not have SNP in some sample in common SNPs site"
#$JAVA -Xmx10g -jar $GATK -T ReduceReads -R $genome -I $BAM -o reduced.bam
#echo "Reduce Done"
echo "Multisample call and filter by hardfilter or VQSR"


###snp call
if [ ! -e $prefix.gatk.snp.raw.vcf ]; then
$JAVA -Xmx10g -jar $GATK \
      -T UnifiedGenotyper \
      --filter_mismatching_base_and_quals \
      -R $genome \
      -I $BAM \
      -o $prefix.gatk.snp.raw.vcf \
      -nct $cpu \
      -stand_call_conf 30 \
      -stand_emit_conf 10 \
      -glm SNP \
      > $prefix.gatk.snp.log 2> $prefix.gatk.snp.log2
fi

###indel call
if [ ! -e $prefix.gatk.indel.raw.vcf ]; then
$JAVA -Xmx10g -jar $GATK \
      -T UnifiedGenotyper \
      --filter_mismatching_base_and_quals \
      -R $genome \
      -I $BAM \
      -o $prefix.gatk.indel.raw.vcf \
      -nct $cpu \
      -stand_call_conf 30 \
      -stand_emit_conf 10 \
      -glm INDEL \
      > $prefix.gatk.indel.log 2> $prefix.gatk.indel.log2
fi

###hardfilter indel
if [ ! -e $prefix.gatk.indel.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.indel.raw.vcf \
      -o $prefix.gatk.indel.hardfilter.vcf \
      --filterExpression "QD < 2.0" \
      --filterName "QDFilter" \
      --filterExpression "ReadPosRankSum < -20.0" \
      --filterName "ReadPosFilter" \
      --filterExpression "FS > 200.0" \
      --filterName "FSFilter" \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
      --filterName "QualFilter"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.indel.hardfilter.vcf -o $prefix.gatk.indel.pass.vcf --excludeFiltered

fi

###hardfilter snp
if [ ! -e $prefix.gatk.snp.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.raw.vcf \
      -o $prefix.gatk.snp.hardfilter.vcf \
      --clusterSize 3 \
      --clusterWindowSize 10 \
      --filterExpression "QD < 2.0" \
      --filterName "QDFilter" \
      --filterExpression "MQ < 40.0" \
      --filterName "MQFilter" \
      --filterExpression "FS > 60.0" \
      --filterName "FSFilter" \
      --filterExpression "HaplotypeScore > 13.0" \
      --filterName "HaplotypeScoreFilter" \
      --filterExpression "MQRankSum < -12.5" \
      --filterName "MQRankSumFilter" \
      --filterExpression "ReadPosRankSum < -8.0" \
      --filterName "ReadPosRankSumFilter" \
      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
      --filterName "StandardFilters" \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --mask $prefix.gatk.indel.pass.vcf \
      --maskName "INDEL"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.snp.hardfilter.vcf -o $prefix.gatk.snp.pass.vcf --excludeFiltered
fi

###VQSR
if [ ! -e $prefix.VQSR.snp.tranches ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantRecalibrator \
      -R $genome \
      -input $prefix.gatk.snp.raw.vcf \
      --maxGaussians 4 \
      --percentBadVariants 0.05 \
      -resource:dbsnp,known=true,training=true,truth=true,prior=10.0 $dbSNP \
      -an QD -an FS -an MQ -an MQRankSum -an HaplotypeScore -an ReadPosRankSum \
      --ts_filter_level 99.0 \
      -mode SNP \
      -recalFile $prefix.VQSR.snp.recal \
      -tranchesFile $prefix.VQSR.snp.tranches \
      -rscriptFile $prefix.VQSR.snp.plots.R
fi

if [ ! -e $prefix.gatk.snp.VQSR.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T ApplyRecalibration \
      -R $genome \
      -input $prefix.gatk.snp.raw.vcf \
      --ts_filter_level 99.0 \
      -tranchesFile $prefix.VQSR.snp.tranches \
      -recalFile $prefix.VQSR.snp.recal \
      -mode SNP \
      -o $prefix.gatk.snp.VQSR.vcf

###mask SNPs cluster and INDEL
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.VQSR.vcf \
      -o $prefix.gatk.snp.VQSR.INDEL.vcf \
      --clusterSize 3 \
      --clusterWindowSize 10 \
      --mask $prefix.gatk.indel.pass.vcf \
      --maskExtension 10 \
      --maskName "INDEL"

###mask repeat
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.VQSR.INDEL.vcf \
      -o $prefix.gatk.snp.VQSR.MASK.vcf \
      --mask $repeat \
      --maskName "REPEAT"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.snp.VQSR.MASK.vcf -o $prefix.gatk.snp.VQSR.pass.vcf --excludeFiltered
fi

###Annotation
REFERENCE="rice7"
DB=/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV/SNPeff/input/
BIN=/rhome/cjinfeng/software/tools/SVcaller/snpEff
INPUT_VCF=$prefix.gatk.snp.VQSR.pass.vcf
VCF_BASE=`basename $INPUT_VCF .vcf`
SNPEFF="$JAVA -Xmx1g -jar $BIN/snpEff.jar"
if [ ! -e $DB/$REFERENCE ]; then
echo "Downloading DB $REFERENCE"
$SNPEFF download -c $BIN/snpEff.config -v $REFERENCE
fi

if [ ! -e $prefix.gatk.snp.VQSR.pass.anno.vcf ]; then

$SNPEFF eff -c $BIN/snpEff.config -v -o gatk $REFERENCE $INPUT_VCF > $VCF_BASE.eff.vcf

$JAVA -Xmx1g -jar $GATK \
      -T VariantAnnotator \
      -R $genome \
      -A SnpEff \
      --variant $prefix.gatk.snp.VQSR.pass.vcf \
      --snpEffFile $prefix.gatk.snp.VQSR.pass.eff.vcf \
      -o $prefix.gatk.snp.VQSR.pass.anno.vcf

fi



###validation

:<<COMMENT
if [ ! -e $prefix.gatk.snp.pass.validation.vcf ]; then

$JAVA -Xmx1g -jar $GATK \
        -T ValidationSiteSelectorWalker \
        -R $genome \
        --variant $prefix.gatk.snp.pass.vcf \
        -o $prefix.gatk.snp.pass.validation.vcf \
        --numValidationSites 100   \
        -sampleMode  POLY_BASED_ON_GT \
        -freqMode KEEP_AF_SPECTRUM \
        -sn HEG4

$JAVA -Xmx1g -jar $GATK \
        -T ValidationSiteSelectorWalker \
        -R $genome \
        --variant $prefix.gatk.indel.pass.vcf \
        -o $prefix.gatk.indel.pass.validation.vcf \
        --numValidationSites 100   \
        -sampleMode  POLY_BASED_ON_GT \
        -freqMode KEEP_AF_SPECTRUM \
        -sn HEG4 

fi
COMMENT

echo "Done!"
