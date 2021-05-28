#!/usr/bin/env nextflow
/*
  Workflow to analyze TRIBE data
  Michael VanInsberghe
*/


def helpMessage() {
  log.info nfcoreHeader()
  log.info"""
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run nf-core/fastqcat --reads '*_R{1,2}.fastq.gz' -profile docker
  Mandatory arguments:
    --reads                       Path to input data (must be surrounded with quotes)
  Options:
    --genomeDir                   Path to previously created genome reference
    --genomeFasta                 Path to genome reference fasta
    --genomeAnnotations           Path to genome annotation gtf
    --genomeVariants              Path to genome variants for filtering
    --depleteFasta                Path to fasta containing sequences to deplete
  """.stripIndent()
}

knownVariants = file(params.knownVariants, checkIfExists:true)
intervalList = file(params.intervalList, checkIfExists:true)
exonsIntervals = file(params.exonsIntervals, checkIfExists:true)
intronsIntervals = file(params.intronsIntervals, checkIfExists:true)
exonsBed = file(params.exonsBed, checkIfExists:true)

if(!params.genomeDir){
  // build genome directory if not supplied

  referenceOut = "${params.outdir}/reference"
  genomeAnnotations = file(params.genomeAnnotations, checkIfExists:true)
  genomeFasta = file(params.genomeFasta, checkIfExists:true)
  if(params.deplete){
    depleteFasta = file(params.depleteFasta, checkIfExists:true)
  }

  process index_genome{
    cpus 1
    time '2h'
    memory '5G'

    publishDir "${referenceOut}/original/", mode:'copy', overwrite: true

    tag "indexing ${params.genomeFasta}"

    input:
    file(genomeFasta)
    file(genomeAnnotations)

    output:
    set file(genomeFasta), file("*.fai") into genome_euk_tRNA, genome_mit_tRNA, genome_prepare, genome_annotations
    file("${genomeAnnotations}")
    file("*.fai") into geome_fai
    file("${genomeFasta}") into genome_fasta

    script:
    """
    samtools faidx ${genomeFasta}
    """
  }

  process eukaryotic_tRNAscan{
    // follows the recoomended paramaters for the updated tRNAscan-SE
    // from https://dx.doi.org/10.1007%2F978-1-4939-9173-0_1
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768409/

    cpus 36
    time '24h'
    memory '50G'

    publishDir "${referenceOut}/tRNAScan/eukaryotic/", mode:'copy', overwrite:true

    tag "Eukaryotic tRNAScan-SE on $params.genomeFasta"

    input:
    set file(genome), file(fai) from genome_euk_tRNA

    output:
    file("${genome.baseName}.Eukaryotic_tRNAs_bp.bed") into eukar_tRNAs_bed
    file("*")

    script:
    """
    tRNAscan-SE \
      -HQ \
      -o# \
      -f# \
      -s# \
      -m# \
      -b# \
      -a# \
      -l# \
      --brief \
      --thread ${task.cpus} \
      -p ${genome.baseName}.Eukaryotic_tRNAs \
      ${genome}

    cat ${genome.baseName}.Eukaryotic_tRNAs.out | tRNAScanToBed.pl > ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
    sed '/^MT\\|^chrM/d' ${genome.baseName}.Eukaryotic_tRNAs_bp.bed > temp.bed
    mv temp.bed ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
    """
  }

  process mitochondrial_tRNAscan{
    // follows the recoomended paramaters for the updated tRNAscan-SE
    // from https://dx.doi.org/10.1007%2F978-1-4939-9173-0_1
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768409/

    cpus 36
    time '24h'
    memory '50G'

    publishDir "${referenceOut}/tRNAScan/mitochondrial/", mode:'copy', overwrite:true

    tag "Eukaryotic tRNAScan-SE on $params.genomeFasta"

    input:
    set file(genome), file(fai) from genome_euk_tRNA

    output:
    file("${genome.baseName}.Mitochondrial_tRNAs_bp.bed") into mito_tRNAs_bed
    file("*")

    script:
    """
    tRNAscan-SE \
      -M vert \
      -Q \
      -o# \
      -f# \
      -m# \
      -b# \
      -a# \
      -l# \
      --brief \
      --thread ${task.cpus} \
      -p ${genome.baseName}.Mitochondrial_tRNAs \
      ${genome}

    cat ${genome.baseName}.Mitochondrial_tRNAs.out | tRNAScanToBed.pl > ${genome.baseName}.Mitochondrial_tRNAs_bp.bed
    """
  }

  process create_genomes{
  
    time '2h'
    memory '50G'
    clusterOptions '--gres=tmpspace:75G' // for sort in mergeGTF
    
    tag "create tRNA-masked genomes for $params.genomeFasta";
    
    publishDir "${referenceOut}/tRNA-masked/", mode:'copy', overwrite: true
    
    input:
    file(eukaryotic) from eukar_tRNAs_bed
    file(mitochondrial) from mito_tRNAs_bed
    set file(genome), file(genome_index) from genome_prepare
    file(annotations) from genomeAnnotations
    
    output:
    file("${genome.baseName}.tRNA_masked_regions.bed")
    // file("${genome.baseName}.tRNA_masked.fa")
    file("${genome.baseName}.pre-tRNAs.bed")
    file("${genome.baseName}.pre-tRNAs.fa")
    file("${genome.baseName}.mature-tRNAs.bed")
    file("${genome.baseName}.mature-tRNAs.fa")
    file("${genome.baseName}.mature-tRNAs.cluster.fa")
    file("${genome.baseName}.tRNACluster.gtf")
    file("${genome.baseName}.tRNA-masked-withmature.fa")
    file("${genome.baseName}.tRNA-masked-withmature.fa.fai")
    file("${genome.baseName}.tRNA-masked-withmature.dict")
    // miR file("${annotations.baseName}.miRNA.tRNA.gtf")
    file("${annotations.baseName}.tRNA.gtf")
    set file("${genome.baseName}.tRNA-masked-withmature.fa"), file("${annotations.baseName}.tRNA.gtf") into star_genome
    set file("${genome.baseName}.tRNA-masked-withmature.fa"), file("${genome.baseName}.tRNA-masked-withmature.fa.fai"), file("${genome.baseName}.tRNA-masked-withmature.dict") into gatk_splitCigar, gatk_bqsr1, gatk_bqsr2, gatk_hap1, gatk_hap2, gatk_annot
    
    script:
    """
    #### tRNAs
    # mask tRNAs in genome
    cat ${eukaryotic} ${mitochondrial} | sort -k 1,1 -k2,2n  > ${genome.baseName}.tRNA_masked_regions.bed
    bedtools maskfasta -fi ${genome} -fo ${genome.baseName}.tRNA_masked.fa -bed ${genome.baseName}.tRNA_masked_regions.bed
    
    ## pre tRNAs
    # remove pseudogenes, add 50 nt 5' and 3' flanking regions
    grep -v "pseudo\\s" ${eukaryotic} | expandBed12.pl --slop 50 --index ${genome_index} > ${eukaryotic.baseName}.pre-tRNAs.bed
    
    # select only MT tRNAs from the MT model, add flaking regions
    sed '/^MT\\|^chrM/!d' ${mitochondrial} > ${mitochondrial.baseName}.only_MT.bed
    grep -v "pseudo\\s" ${mitochondrial.baseName}.only_MT.bed | expandBed12.pl --slop 50 --index ${genome_index} > ${mitochondrial.baseName}.pre-tRNAs.bed
    
    # combine to pre-tRNAs.bed
    cat ${eukaryotic.baseName}.pre-tRNAs.bed ${mitochondrial.baseName}.pre-tRNAs.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.pre-tRNAs.bed
    
    # extract pre-tRNA sequences
    bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.pre-tRNAs.bed -fo ${genome.baseName}.pre-tRNAs.fa
    
    ## mature tRNAs
    cat ${eukaryotic} ${mitochondrial.baseName}.only_MT.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.mature-tRNAs.bed
    bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.mature-tRNAs.bed | appendFasta.pl --append cca > ${genome.baseName}.mature-tRNAs.fa
  
    ## cluster mature tRNAs
    collapseSequences.pl ${genome.baseName}.mature-tRNAs.fa > ${genome.baseName}.mature-tRNAs.cluster.fa
    # produces cluster_info.gtf
    mv cluster_info.gtf ${genome.baseName}.tRNACluster.gtf
    samtools faidx ${genome.baseName}.mature-tRNAs.cluster.fa
    
    ## assemble final genome
    cat ${genome.baseName}.tRNA_masked.fa ${genome.baseName}.mature-tRNAs.cluster.fa > ${genome.baseName}.tRNA-masked-withmature.fa
    samtools faidx ${genome.baseName}.tRNA-masked-withmature.fa

    java -jar ${params.GATK} CreateSequenceDictionary \
      --REFERENCE ${genome.baseName}.tRNA-masked-withmature.fa 
    
    #### Annotations
    mergeGTF.pl ${annotations} ${genome.baseName}.tRNACluster.gtf > ${annotations.baseName}.tRNA.gtf
    
    # geneinfo file for later analysis
    #extractGeneInfo.pl ${annotations.baseName}.tRNA.gtf > ${annotations.baseName}.tRNA.geneinfo
    """
  }

  process star_reference{
  
    cpus 24
    time '12h'
    memory '100G'
    
    tag "STAR reference"
    
    publishDir "${referenceOut}/", mode:'copy', overwrite: true
    
    input:
    set file(masked_genome), file(annotations) from star_genome
    
    output:
    file("star_tRNAmasked_${params.starOverhang}")
    file("star_tRNAmasked_${params.starOverhang}") into starIndex
    
    script:
    """
    mkdir ./star_tRNAmasked_${params.starOverhang}
    
    STAR \
    	--runMode genomeGenerate \
    	--runThreadN ${task.cpus} \
    	--genomeDir ./star_tRNAmasked_${params.starOverhang} \
    	--genomeFastaFiles ${masked_genome} \
    	--limitGenomeGenerateRAM 75161927680 \
    	--sjdbGTFfile ${annotations} \
    	--sjdbOverhang ${params.starOverhang}
    """
  }

  process extract_annotations{

    cpus 1
    time '5h'
    memory '75G'

    tag "Extract annotations for ${params.genomeFasta}"


    beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
    afterScript "deactivate"
    
    publishDir "${referenceOut}/annotations/", mode:'copy', overwrite: true
    
    input:
    set file(genome), file(genome_index) from genome_annotations
    file(annotations) from genomeAnnotations
    
    output:
    file(annotations)
    file("*.{csv.gz,feather}")
    file("*.pickle")
    file("${annotations.baseName}.annotations.{csv.gz,feather}") into annotationFeather
    file("*.annotations.pickle") into annotationPickle
    
    
    script:
    """
    extractAnnotations.py -g ${annotations}
    extractCodons.py -f ${genome} -t ${annotations.baseName}.annotations.pickle
    """
  }


  process depletion_reference{

    cpus 12
    time '12h'
    memory '50G'

    tag 'Reference for depletion'

    publishDir "${referenceOut}/depletion/", mode:'copy', overwrite: true

    input:
    file(depleteFasta)

    output:
    file("bwa_depletion")
    file("bwa_depletion") into depletionIndex


    script:
    """
    mkdir ./bwa_depletion
    bwa index -p bwa_depletion/depletion $depleteFasta
    """
  }
}else{
  // genome has already been prepared, set variables
  // annotationFeather = file("${params.genomeDir}/annotations/*.annotations.{csv.gz,feather}", checkIfExists: true).first()
  // annotationPickle = file("${params.genomeDir}/annotations/*.annotations.pickle", checkIfExists: true).first()
  // genomeFasta = file("${params.genomeDir}/original/*.{fa,fasta}", checkIfExists: true).first()
  // genomeFai =  file("${params.genomeDir}/original/*.fai", checkIfExists: true).first()
  // starIndex = file("${params.genomeDir}/star_tRNAmasked_${params.starOverhang}", checkIfExists: true)

  Channel
    .fromPath("${params.genomeDir}/annotations/*.annotations.{csv.gz,feather}")
    .first()
    .ifEmpty{ exit 1, "Cannot find annotation csv/feather in: ${params.genomeDir}/annotations\nGenome directory not properly formed" }
    .set{annotationFeather}

  Channel
    .fromPath("${params.genomeDir}/annotations/*.annotations.pickle")
    .first()
    .ifEmpty{ exit 1, "Cannot find annotation pickle in: ${params.genomeDir}/annotations\nGenome directory not properly formed" }
    .set{annotationPickle}

  Channel
    .fromPath("${params.genomeDir}/original/*.{fa,fasta}")
    .first()
    .ifEmpty{ exit 1, "Cannot find genome fasta in ${params.genomeDir}/original\nGenome directory not properly formed" }
    .set{genomeFasta}

  Channel
    .fromPath("${params.genomeDir}/original/*.fai")
    .first()
    .ifEmpty{ exit 1, "Cannot find genome index in ${params.genomeDir}/original\nGenome directory not properly formed" }
    .set{genomeFai}

  Channel
    .fromPath("${params.genomeDir}/star_tRNAmasked_${params.starOverhang}")
    .ifEmpty{ exit 1, "Cannot find STAR index with overhang ${params.starOverhang} in: ${params.genomeDir}\nGenome directory not properly formed" }
    .set{starIndex}

    Channel
      .fromPath("${params.genomeDir}/rRNA_depletion")
      .ifEmpty{ exit 1, "Cannot find rRNA depletion index in ${params.genomeDir}\nGenome directory not properly formed" }
      .set(rRNAIndex)
    
    Channel
      .fromPath("${params.genomeDir}/rRNA_depletion")
      .ifEmpty{ exit 1, "Cannot find rRNA depletion index in ${params.genomeDir}\nGenome directory not properly formed" }
      .set(tRNAIndex)




  if(params.deplete){
    depletionIndex = file("${params.genomeDir}/bwa_depletion", checkIfExists: true)
  }
}

if(!params.bulk){
  whitelist = file(params.whitelist, checkIfExists: true)
}


Channel
  .fromFilePairs (params.reads, size: params.bulk ? 1:2)
  .ifEmpty{ exit 1, "Cannot find any reads matching: ${params.reads}\nPath needs to be enclosed in quotes!\nPath requires at least one * wildcard!" }
  .into { raw_reads_fastqc; raw_reads_trim }

// FastQCs
// raw, trimmed, and depleted if deplete
process raw_fastqc {
  cpus 6
  time '5h'
  memory '5G'
  
  tag "FastQC on raw $library"
  publishDir "${params.outdir}/QC/raw_fastqc", mode:'copy', overwrite: true
  
  input:
  set name, file(reads) from raw_reads_fastqc
  
  output:
  file("*_fastqc.{zip,html}") into multiqc_fastqc
  
  script:
  library = name.toString().tokenize('_').get(0)
  """
  fastqc -t ${task.cpus} -f fastq -q ${reads} --nogroup
  """
}

if(params.trimMethod){
  process trim {
    time '5h'
    cpus 8
    memory '10G'
    
    tag "trimming $library"
    
    beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv36/bin/activate"
    afterScript "deactivate"
    
    input:
    set name, file(reads) from raw_reads_trim
    // set library, file(read1), file(read2) from extracted_fastq
    
    output:
    set library, file("*.trimmed.fastq.gz") into trimmed_fastq, trimmed_reads_fastqc
    //set library, file("${library}_R1_trimmed.fastq.gz"), file("${library}_R2_trimmed.fastq.gz") into trimmed_fastq
    file("${library}.trim_report.txt") into multiqc_cutadapt
    
    script:
    library = name.toString().tokenize('_').get(0)
    template "${params.trimMethod}"
  
  }
}else{
  raw_reads_trim.set{ trimmed_reads }
}

if(params.trimMethod){
  process trimmed_fastqc {
    cpus 6
    time '5h'
    memory '5G'
    
    tag "FastQC on trimmed $library"
    publishDir "${params.outdir}/QC/trimmed_fastqc", mode:'copy', overwrite: true
    
    input:
    set library, file(reads) from trimmed_reads_fastqc
    
    output:
    file("*_fastqc.{zip,html}") into multiqc_trimmed_fastqc
    
    script:
    """
    fastqc -t ${task.cpus} -f fastq -q ${reads} --nogroup
    """
  }
}else{
  multiqc_trimmed_fastqc = Channel.empty()
}


if (!params.deplete){
  trimmed_reads .set { depleted_reads }
  // multiqc_bwa_depletion = Channel.empty()
  depleted_reads_fastqc = Channel.empty()
}else{
  process bwa_depletion {
    cpus 12
    time '4h'
    memory '50G'
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"
    
    tag "depleting contaminant sequences on $library"
    
    input:
    file(depletionIndex)
    set library, file(reads) from trimmed_fastq
    
    output:
    set library, file("*.depleted.fastq.gz") into depleted_reads_first, depleted_reads_second, depleted_reads_fastqc
    // LOGFILES set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_bwa_depletion
    
    script:
    template "${params.depleteMethod}"
  }
}

if(params.deplete){
  process depleted_fastqc {
    cpus 6
    time '5h'
    memory '5G'
    
    tag "FastQC on depleted $library"
    publishDir "${params.outdir}/QC/depleted_fastqc", mode:'copy', overwrite: true
    
    input:
    set library, file(reads) from depleted_reads_fastqc
    
    output:
    file("*_fastqc.{zip,html}") into multiqc_depleted_fastqc
    
    script:
    """
    fastqc -t ${task.cpus} -f fastq -q ${reads} --nogroup
    """
  }
}else{
  multiqc_depleted_fastqc = Channel.empty()
}

// Alignment 
if(params.bulk){
  process STAR{
    cpus 24
    time '4h'
    memory '75G'
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"
    
    tag "STAR Alignment ${library}"
    
    publishDir "${params.outdir}/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy', overwrite: true
    publishDir "${params.outdir}/alignment", pattern: "*.ba*", mode:'copy', overwrite: true
    
    input:
    file starIndex
    set library, file(reads) from depleted_reads_first 
    
    output:
    set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star
    set library, file("${library}_Aligned.sortedByCoord.out.bam"), file("${library}_Aligned.sortedByCoord.out.bam.bai") into aligned_genome
    set library, file("${library}_Aligned.toTranscriptome.out.bam"), file("${library}_Aligned.toTranscriptome.out.bam.bai") into aligned_transcriptome
    	
    script:
    """
    STAR \
    	${params.STARmapParams} \
    	--runThreadN ${task.cpus} \
    	--genomeDir ${starIndex} \
    	--readFilesIn ${reads} \
    	--outFileNamePrefix ${library}_ \
    	--outSAMattributes ${params.STARoutSamAttributes}
    
    samtools index ${library}_Aligned.sortedByCoord.out.bam
    
    sambamba sort \
    	-m 2G \
    	-t ${task.cpus} \
    	${library}_Aligned.toTranscriptome.out.bam
    
    rm ${library}_Aligned.toTranscriptome.out.bam
    
    mv ${library}_Aligned.toTranscriptome.out.sorted.bam ${library}_Aligned.toTranscriptome.out.bam
    mv ${library}_Aligned.toTranscriptome.out.sorted.bam.bai ${library}_Aligned.toTranscriptome.out.bam.bai
    
    """
  }
}else{
  process STARsolo_first {
    cpus 24
    time '4h'
    memory '85G'
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"
    
    tag "STARSolo first pass $library"
    
    beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
    afterScript "deactivate"
    
    input:
    file whitelist
    file starIndex
    set library, file(reads) from depleted_reads_first
    
    output:
    file("*_SJ.out.tab") into sjtab_first

    script:
    if(params.cDNARead == "R2"){
      readOrder = "${reads[1]} ${reads[0]}"
    }else if(params.cDNARead == "R1"){
      readOrder = "${reads[0]} ${reads[1]}"
    }
    """
    STAR \
      ${params.STARmapParams} \
      --runThreadN ${task.cpus} \
      --genomeDir ${starIndex} \
      --readFilesIn ${readOrder} \
      --outFileNamePrefix ${library}_ \
      --outSAMattributes ${params.STARoutSamAttributes} ${params.STARsoloSamAttributes} \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist ${whitelist} \
      --soloCBstart ${params.CBstart} \
      --soloCBlen ${params.CBlen} \
      --soloUMIstart ${params.UMIstart} \
      --soloUMIlen ${params.UMIlen} \
      --soloBarcodeReadLength 0 \
      --soloStrand ${params.readStrand} \
      --soloFeatures Gene GeneFull SJ Velocyto \
      --soloUMIdedup 1MM_Directional \
      --soloCellFilter None
  
    rm -rf *.bam 
    rm -rf *.bai
    rm -rf *.out
    rm -rf *STARtmp
    
    """
  }

  process STARsolo_second {
    cpus 24
    time '4h'
    memory '85G'
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"
    
    tag "STARSolo second pass $library"
    
    beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
    afterScript "deactivate"
    
    publishDir "${params.outdir}/alignment", pattern: "*.ba*", mode:'copy', overwrite: true
    publishDir "${params.outdir}/quantification", pattern: "*_Solo.out", mode:'copy', overwrite: true
    publishDir "${params.outdir}/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy' , overwrite: true
    
    input:
    file whitelist
    file starIndex
    set library, file(reads) from depleted_reads_second
    file sjtabs from sjtab_first.collect()
    
    output:
    file("*_Solo.out")
    set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star
    set library, file("${library}_Aligned.sortedByCoord.out.bam"), file("${library}_Aligned.sortedByCoord.out.bam.bai") into aligned_genome
    set library, file("${library}_Aligned.toTranscriptome.out.bam"), file("${library}_Aligned.toTranscriptome.out.bam.bai") into aligned_transcriptome
    	
    script:
    if(params.cDNARead == "R2"){
      readOrder = "${reads[1]} ${reads[0]}"
    }else if(params.cDNARead == "R1"){
      readOrder = "${reads[0]} ${reads[1]}"
    }
    """
    # filter from https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I?pli=1
    cat ${sjtabs} | awk '(\$5 > 0 && \$7 > 2 && \$6==0)' | cut -f1-6 | sort | uniq > merged_SJ.tab

    STAR \
      ${params.STARmapParams} \
      --sjdbFileChrStartEnd merged_SJ.tab \
      --runThreadN ${task.cpus} \
      --genomeDir ${starIndex} \
      --readFilesIn ${readOrder} \
      --outFileNamePrefix ${library}_ \
      --outSAMattributes ${params.STARoutSamAttributes} ${params.STARsoloSamAttributes} \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist ${whitelist} \
      --soloCBstart ${params.CBstart} \
      --soloCBlen ${params.CBlen} \
      --soloUMIstart ${params.UMIstart} \
      --soloUMIlen ${params.UMIlen} \
      --soloBarcodeReadLength 0 \
      --soloStrand ${params.readStrand} \
      --soloFeatures Gene GeneFull SJ Velocyto \
      --soloUMIdedup 1MM_Directional \
      --soloCellFilter None

    #addRG.pl --bam ${library}_Aligned.sortedByCoord.out.bam
    #rm ${library}_Aligned.sortedByCoord.out.bam
    #mv ${library}_Aligned.sortedByCoord.out.rg.bam ${library}_Aligned.sortedByCoord.out.bam
    
    sambamba index -t ${task.cpus} ${library}_Aligned.sortedByCoord.out.bam
    
    addCB.py \
      -b ${library}_Aligned.toTranscriptome.out.bam \
      -t ${task.cpus} \
      -c ${whitelist} 
    
    sambamba sort \
      -m 2G \
      -t ${task.cpus} \
      ${library}_Aligned.toTranscriptome.out_CB.bam
    
    rm ${library}_Aligned.toTranscriptome.out.bam
    rm ${library}_Aligned.toTranscriptome.out_CB.bam
    
    mv ${library}_Aligned.toTranscriptome.out_CB.sorted.bam ${library}_Aligned.toTranscriptome.out.bam
    mv ${library}_Aligned.toTranscriptome.out_CB.sorted.bam.bai ${library}_Aligned.toTranscriptome.out.bam.bai
    
    find . -type f \\( -name '*.tsv' -o -name '*.mtx' \\) -exec gzip "{}" \\;
    """
  }
}

process dedup_genome {
  memory '200G'
  time '50h'
  // clusterOptions "--gres=tmpspace:50G"
  
  tag "umi_tools dedup on $library"
  
  beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
  afterScript "deactivate"
  
  input:
  set library, file(bam), file(bai) from aligned_genome
  
  output:
  set library, file("${bam.baseName}.dedup.bam") into deduplicated_genome
  file("*_dedup.log") into multiqc_deduplicate_genome
  
  script:
  template "${params.deduplicateGenomeMethod}"
}

process sort_dedup_genome {
  memory '50G'
  cpus 12
  time '2h'
  clusterOptions "--gres=tmpspace:50G"
  
  tag "sorting deduped $library"

  publishDir "${params.outdir}/deduplicate", mode:'copy', overwrite: true

  input:
  set library, file(bam) from deduplicated_genome

  output:
  set file("${bam.baseName}.sorted.bam"), file("${bam.baseName}.sorted.bam.bai") into sorted_deduplicated_genome
  file("${bam.baseName}.sorted.bam")
  file("${bam.baseName}.sorted.bam.bai")
  
  script:
  """
  sambamba sort \
    -m 4G \
    -t ${task.cpus} \
    ${bam}

  """
}

process merge_genome{
  memory '50G'
  cpus 12
  time '2h'
  clusterOptions "--gres=tmpspace:50G"
  
  tag "merging genome alignments"

  input:
  file("*") from sorted_deduplicated_genome.collect()

  output:
  set file("merged.bam"), file("merged.bam.bai") into merged_genome
  
  script:
  """
  sambamba merge \
    -t ${task.cpus} \
    merged.bam \
    *.bam
  """
}



process scatter_splitCigar {
  memory '10G'
  time '1h'

  tag "scatter splitCigar"

  input:
  set file(bam), file(bai) from merged_genome 

  output:
  file("split_scatter/*") into splitCigar_scatter

  script:
  """
  mkdir split_scatter

  java -jar ${params.GATK} SplitSamByNumberOfReads \
    --INPUT ${bam} \
    --OUTPUT split_scatter/ \
    --SPLIT_TO_N_FILES ${params.scatterSplitCigarFiles}

  """
}

process splitCigar {
  memory '10G'
  time '2h'
  clusterOptions "--gres=tmpspace:10G"

  tag "splitCigar on ${scatter}"

  input:
  file(scatter) from splitCigar_scatter.flatten()
  set file(genomeFasta), file(genomeFai), file(genomeDict) from gatk_splitCigar

  output:
  file("*.split") into splitCigar_gather

  script:
  """
  java -jar ${params.GATK} SplitNCigarReads \
    --reference ${genomeFasta} \
    --input ${scatter} \
    --output ${scatter}.split \
    --tmp-dir \$TMPDIR

  """
}

process gather_splitCigar{
  memory '50G'
  cpus 12
  time '2h'
  clusterOptions "--gres=tmpspace:50G"
  
  tag "merging genome alignments"
  publishDir "${params.outdir}/gatk/splitCigar", mode:'copy', overwrite: true

  input:
  file("*") from splitCigar_gather.collect()

  output:
  set file("merged_split.bam"), file("merged_split.bam.bai") into gathered_splitCigar1, gathered_splitCigar2
  file("merged_split.bam")
  file("merged_split.bam.bai")
  
  script:
  """
  sambamba merge \
    -t ${task.cpus} \
    merged_split.bam \
    *.bam.split
  """
}

process scatter_intervals {
  cpus 1
  memory '5G'
  time '1h'

  tag "scatter intervals ${intervalList}"

  input:
  file(intervalList) from intervalList

  output:
  file("split_intervals/*") into intervals_hap1, intervals_hap2

  script:
  """
  mkdir split_intervals

  java -jar ${params.GATK} IntervalListTools \
    --SCATTER_COUNT ${params.scatterHaplotype} \
    --SUBDIVISION_MODE INTERVAL_SUBDIVISION \
    --UNIQUE true \
    --SORT true \
    --INPUT ${intervalList} \
    --OUTPUT split_intervals/
  """
}

process index_knownVariants{
  cpus 1
  memory '10G'
  time '1h'

  tag "indexing knownVariants"

  input:
  file(knownVariants) from knownVariants

  output:
  set file("${knownVariants}"), file("${knownVariants}.tbi") into knownVar_bqsr1, knownVar_bqsr2, knownVar_annot

  script:
  """
  tabix -p vcf ${knownVariants}
  """
}

process bqsr_first {
  cpus 1
  memory '50G'
  time '5h'
  clusterOptions "--gres=tmpspace:50G"
  
  tag "BQSR round 1"
  publishDir "${params.outdir}/gatk/bqsr", pattern: "*.pdf", mode:'copy', overwrite: true

  input:
  set file(bam), file(bai) from gathered_splitCigar1
  set file(genomeFasta), file(genomeFai), file(genomeDict) from gatk_bqsr1
  set file(knownVariants), file(knownVarIndex) from knownVar_bqsr1

  output:
  set file("*.bqsr1.bam"), file("*.bqsr1.bai") into recalibrated_first
  file("bqsr_first_AnalyzeCovariates.pdf")

  script:
  knownSites = "--known-sites ${knownVariants}" 
  outputCovariates = "bqsr_first_AnalyzeCovariates.pdf"
  template 'gatk/bqsr.sh'
}

process haplotype_first {
  cpus 1
  memory '50G'
  time '5h'
  clusterOptions "--gres=tmpspace:5G"

  tag "HaplotypeCaller round 1 ${intervalList}"

  input:
  set file(bam), file(bai) from recalibrated_first
  set file(genomeFasta), file(genomeFai), file(genomeDict) from gatk_hap1
  file(intervalList) from intervals_hap1.flatten()

  output:
  file("${intervalList.baseName}_first.vcf") into hap_first

  script:
  outputVCF = "${intervalList.baseName}_first.vcf"
  template 'gatk/haplotypeCaller.sh'

}

process merge_first {
  cpus 1
  memory '50G'
  time '1h'
  clusterOptions "--gres=tmpspace:50G"

  tag "Merging first round VCF"
  publishDir "${params.outdir}/gatk/haplotype", mode:'copy', overwrite: true

  input:
  file(vcf) from hap_first.collect()

  output:
  file("haplotype_first_merged.vcf.gz") 
  set file("haplotype_first_merged.vcf.gz"), file("haplotype_first_merged.vcf.gz.tbi") into hap_first_merged

  script:
  outputVCF = "haplotype_first_merged.vcf.gz"
  template 'gatk/mergeVCF.sh'

}

process bqsr_second {
  cpus 1
  memory '50G'
  time '5h'
  clusterOptions "--gres=tmpspace:50G"
  
  tag "BQSR round 2"
  publishDir "${params.outdir}/gatk/bqsr", pattern: "*.pdf", mode:'copy', overwrite: true

  input:
  set file(bam), file(bai) from gathered_splitCigar2
  set file(genomeFasta), file(genomeFai), file(genomeDict) from gatk_bqsr2
  set file(firstVariants), file(firstVarIndex) from hap_first_merged
  set file(knownVariants), file(knownVarIndex) from knownVar_bqsr2


  output:
  set file("*.bqsr1.bam"), file("*.bqsr1.bai") into recalibrated_second
  file("bqsr_second_AnalyzeCovariates.pdf")

  script:
  knownSites = "--known-sites ${knownVariants} --known-sites ${firstVariants}"
  outputCovariates = "bqsr_second_AnalyzeCovariates.pdf"
  template 'gatk/bqsr.sh'
}

process haplotype_second {
  cpus 1
  memory '50G'
  time '5h'
  clusterOptions "--gres=tmpspace:5G"

  tag "HaplotypeCaller round 2 ${intervalList}"

  input:
  set file(bam), file(bai) from recalibrated_second
  set file(genomeFasta), file(genomeFai), file(genomeDict) from gatk_hap2
  file(intervalList) from intervals_hap2.flatten()

  output:
  file("${intervalList.baseName}_second.vcf") into hap_second

  script:
  outputVCF = "${intervalList.baseName}_second.vcf"
  template 'gatk/haplotypeCaller.sh'

}

process merge_second {
  cpus 1
  memory '50G'
  time '1h'
  clusterOptions "--gres=tmpspace:50G"

  tag "Merging second round VCF"
  publishDir "${params.outdir}/gatk/haplotype", mode:'copy', overwrite: true

  input:
  file(vcf) from hap_second.collect()

  output:
  file("haplotype_second_merged.vcf.gz") 
  set file("haplotype_second_merged.vcf.gz"), file("haplotype_second_merged.vcf.gz.tbi") into hap_second_merged

  script:
  outputVCF = "haplotype_second_merged.vcf.gz"
  template 'gatk/mergeVCF.sh'

}


process select_annotate_variants{
  cpus 1
  memory '20G'
  time '10h'
  clusterOptions "--gres=tmpspace:5G"

  tag "Selecting variants"

  publishDir "${params.outdir}/gatk/haplotype", mode:'copy', overwrite: true

  input:
  set file(variants), file(variantsIndex) from hap_second_merged
  file(exonsIntervals) from exonsIntervals
  file(intronsIntervals) from intronsIntervals
  set file(genomeFasta), file(genomeFai), file(genomeDict) from gatk_annot
  file(exonsBed) from exonsBed
  set file(knownVariants), file(knownVarIndex) from knownVar_annot

  output:
  file("haplotype_merged_known_filtered_introns.vcf.gz")
  file("haplotype_merged_filtered_exons_transcriptome_annot.vcf.gz")
  file("haplotype_merged_filtered_exons_transcriptome_annot_SB.txt.gz")


  script:
  """
  ## select exonic and intronic variants
  java -jar ${params.GATK} SelectVariants \
    -V ${variants} \
    --intervals ${exonsIntervals} \
    -O temp_filtered_exons.vcf.gz 

  java -jar ${params.GATK} SelectVariants \
    -V ${variants} \
    --intervals ${intronsIntervals} \
    -O temp_filtered_introns.vcf.gz

  java -jar ${params.GATK} SelectVariants \
    -V temp_filtered_exons.vcf.gz \
    --discordance ${knownVariants} \
    -O temp_known_filtered_exons.vcf.gz
  
  java -jar ${params.GATK} SelectVariants \
    -V temp_filtered_introns.vcf.gz \
    --discordance ${knownVariants} \
    -O haplotype_merged_known_filtered_introns.vcf.gz

  ## annotate
  echo '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand from bed/gtf file">' > annot.hdr
  echo '##INFO=<ID=GENE,Number=1,Type=String,Description="ENSG__Gene__Biotype from bed/gtf file">' >> annot.hdr

  tabix -p bed ${exonsBed}
  bcftools annotate \
    -a  ${exonsBed}\
    -h annot.hdr \
    -c CHROM,FROM,TO,STRAND,GENE temp_known_filtered_exons.vcf.gz | \
  gzip -c > temp_transannotated.vcf.gz

  annotateSequenceContext.pl \
    -r ${genomeFasta} \
    -v temp_transannotated.vcf.gz \
    -o haplotype_merged_filtered_exons_transcriptome_annot.vcf.gz

  java -jar ${params.GATK} VariantsToTable \
    -V haplotype_merged_filtered_exons_transcriptome_annot.vcf.gz \
    -F CHROM -F POS -F REF -F ALT -F TYPE -F FILTER -F QUAL -F STRAND -F GENE -F CONTEXT -GF SB \
    -O haplotype_merged_filtered_exons_transcriptome_annot_SB.txt \
    --moltenize true

  sed -i '/\\tNA\$/d' haplotype_merged_filtered_exons_transcriptome_annot_SB.txt

  gzip haplotype_merged_filtered_exons_transcriptome_annot_SB.txt

  rm temp_*

  """
}

process multiQC {
  time '1h'
  memory '5G'

  publishDir "${params.outdir}/QC", mode:'copy', overwrite: true

  beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/mqcdev/bin/activate"
  afterScript "deactivate"
  
  input:
  file('*') from multiqc_fastqc.collect().ifEmpty([])
  file('*') from multiqc_trimmed_fastqc.collect().ifEmpty([])
        file('*') from multiqc_depleted_fastqc.collect().ifEmpty([]) 
  file('*') from multiqc_cutadapt.collect().ifEmpty([])
  file('*') from multiqc_star.collect().ifEmpty([])
  file('*') from multiqc_deduplicate_genome.collect().ifEmpty([])

  output:
  file "multiqc_report.html"
  file "multiqc_data"

  script:
  """
  multiqc .
  """
}

