process BWA {
    cpus params.threads_bwa
    memory params.mem_bwa
    maxForks params.maxfork_bwa
    input:
    tuple val(SEQ), val(seq), val(samp), path(trim_p_fwd), path(trim_p_rev), path(trim_u_fwd), path(trim_u_rev), val(lane)

    output:
    tuple val(SEQ), path("${seq}_L${lane}_P.bam"), path("${seq}_L${lane}_1U.bam"), path("${seq}_L${lane}_2U.bam"), emit: bams
    tuple val(SEQ), val(seq), path("${seq}_L${lane}_P.bam"), val(lane), emit: fs_in
    
    script:
    """
    RG=\$(echo -e "@RG\tID:${seq}:${lane}\tPL:illumina\tLB:${samp}\tSM:${samp}")
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_p_fwd} ${trim_p_rev} |
    samtools collate -@ ${task.cpus} -O -u - |
    samtools fixmate -@ ${task.cpus} -m -u - - |
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_P.bam

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_u_fwd} |
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_1U.bam

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_u_rev} |
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_2U.bam
    """

    stub:
    """
    RG=\$(echo -e "@RG\tID:${seq}:${lane}\tPL:illumina\tLB:${samp}\tSM:${samp}")
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    bwa-mem2 mem \
    -R "\${RG}" \
    -t 1 \
    ${params.bwaIndex} \
    ${trim_p_fwd} ${trim_p_rev} |
    samtools collate -@ ${task.cpus} -O -u - |
    samtools fixmate -@ ${task.cpus} -m -u - - |
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_P.bam
    touch ${seq}_L${lane}_1U.bam
    touch ${seq}_L${lane}_2U.bam
    """
}