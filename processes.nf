process PETRIM {
	cpus = params.threads_trim
	input:
	tuple val(SEQ), val(seq), path(fastqs), val(lane)
	
	output:
	tuple val(SEQ), path("*_t_*P.fastq.gz"), path("*_t_1U.fastq.gz"), path("*_t_2U.fastq.gz"), val(lane) emit: trims
	tuple val(SEQ), path("*.tstat"), val(lane) emit: tstat
	
	script:
	"""	
	trimmomatic PE \
	-threads ${task.cpus} \
	-summary ${seq}_L${lane}.tstat \
	${fastqs} \
	-baseout ${seq}_L${lane}_t.fastq.gz \
	HEADCROP:1 ILLUMINACLIP:/nfs/home/lmajeres/tools/adapters.fa:2:30:10 \
    LEADING:20 SLIDINGWINDOW:15:20 MINLEN:75
	"""
}

process BWA {
    cpus = params.threads_bwa
    memory = params.mem_bwa
    input:
    tuple val(SEQ), val(seq), val(samp), path(trim_pair), path(trim_fwd), path(trim_rev), val(lane)

    output:
    tuple val(SEQ), path("*.bam"), val(lane) emit: bams
    tuple val(SEQ), val(seq), path("*P.bam"), val(lane) emit: fs_in
    
    script:
    """
    RG=\$(echo -e "@RG\tID:${seq}:${lane}\tPL:illumina\tLB:${samp}\tSM:${samp}")

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_pair} |
    samtools collate -@ ${task.cpus} -O -u - |
    samtools fixmate -@ ${task.cpus} -m -u - - |
    samtools sort -@ ${task.cpus} -u - ${seq}_L${lane}_P.bam


    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_fwd} |
    samtools sort -@ ${task.cpus} -u - ${seq}_L${lane}_1U.bam

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_rev} |
    samtools sort -@ ${task.cpus} -u - ${seq}_L${lane}_2U.bam
    """
}

process FLAGSTAT {
    cpus = params.threads_fs
    input:
    tuple val(SEQ), val(seq), path(pair_bam), val(lane)

    output:
    tuple val(SEQ), path("*.flagstat"), val(lane)

    script:
    """
    samtools flagstat -@ ${task.cpus} ${pair_bam} > ${seq}_L${lane}.flagstat
    """
}