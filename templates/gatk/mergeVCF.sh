
echo "${vcf}" | tr " " "\\n" > input.list

java -jar ${params.GATK} MergeVcfs \
  --INPUT input.list \
  --OUTPUT ${outputVCF} \
  --TMP_DIR \$TMPDIR

rm input.list
