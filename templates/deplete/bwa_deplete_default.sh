bwa aln -t ${task.cpus} ${depletionIndex}/depletion ${library}_R2.trimmed.fastq.gz > ${library}_aln.sai
bwa samse ${depletionIndex}/depletion ${library}_aln.sai ${library}_R2.trimmed.fastq.gz | samtools view -F 4 | awk '{print \$1}' > ${library}_aln-hits.txt

bwa mem -t ${task.cpus} ${depletionIndex}/depletion ${library}_R2.trimmed.fastq.gz | samtools view -F 4 | awk '{print \$1}' > ${library}_mem-hits.txt

cat ${library}_aln-hits.txt ${library}_mem-hits.txt | sort | uniq | gzip -c > ${library}_depletion-hits.txt.gz

rm ${library}_aln-hits.txt ${library}_mem-hits.txt
rm ${library}_aln.sai

filterFastq.pl --reads ${reads} --hits ${library}_depletion-hits.txt.gz
