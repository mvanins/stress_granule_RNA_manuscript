java -jar ${params.GATK} HaplotypeCaller \
  --reference ${genomeFasta} \
  --input ${bam} \
  --dont-use-soft-clipped-bases true \
  --standard-min-confidence-threshold-for-calling 20.0 \
  --output haplotype_0.vcf \
  --native-pair-hmm-threads 2 \
  --verbosity ERROR \
  -L ${intervalList}/scattered.interval_list \
  --interval-padding ${params.haplotypePadding} \
  --annotation StrandBiasBySample \
  --tmp-dir \$TMPDIR

java -jar ${params.GATK} SelectVariants \
  -V haplotype_0.vcf \
  --intervals ${intervalList}/scattered.interval_list \
  -O ${outputVCF} \
  --tmp-dir \$TMPDIR

rm haplotype_0.vcf
