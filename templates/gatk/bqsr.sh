java -jar ${params.GATK} BaseRecalibrator \
  --input ${bam} \
  --reference ${genomeFasta} \
  ${knownSites} \
  --output bqsr_0.table \
  --tmp-dir \$TMPDIR


java -jar ${params.GATK} ApplyBQSR \
  -I ${bam} \
  --bqsr-recal-file bqsr_0.table \
  -O ${bam.baseName}.bqsr1.bam \
  --tmp-dir \$TMPDIR


java -jar ${params.GATK} BaseRecalibrator \
  --input ${bam.baseName}.bqsr1.bam \
  --reference ${genomeFasta} \
  ${knownSites} \
  --output bqsr_1.table \
  --tmp-dir \$TMPDIR

java -jar ${params.GATK} AnalyzeCovariates \
  -before bqsr_0.table \
  -after bqsr_1.table \
  -plots ${outputCovariates} \
  --tmp-dir \$TMPDIR

