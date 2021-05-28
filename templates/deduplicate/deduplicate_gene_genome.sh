umi_tools dedup \
  -I ${bam} \
  -S ${bam.baseName}.dedup.bam \
  -L ${bam.baseName}_dedup.log \
  --spliced-is-unique \
  --per-gene \
  --gene-tag=GX \
  --per-cell \
  --extract-umi-method tag \
  --umi-tag=UR \
  --cell-tag=CB \
  --read-length \
  --no-sort-output

addRG.pl --bam ${bam.baseName}.dedup.bam

rm ${bam.baseName}.dedup.bam
mv ${bam.baseName}.dedup.rg.bam ${bam.baseName}.dedup.bam
