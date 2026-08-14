process MERGEMD {
    cpus params.threads_merge
    input:
    tuple val(samp), val(dv_flag), path(bam_p), path(bam_1u), path(bam_2u)

    output:
    tuple val(samp), val(dv_flag), path("${samp}_MD.bam"), emit: mdbam
    tuple val(samp), path("${samp}_MD.stat"), emit: mdstat

    script:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    samtools merge -@ ${task.cpus} -u -o - ${bam_p} ${bam_1u} ${bam_2u} |
    samtools markdup -@ ${task.cpus} -S -d 2500 \
    -s -f ${samp}_MD.stat - ${samp}_MD.bam
    """

    stub:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    samtools markdup -S -d 2500 -s -f ${samp}_MD.stat ${bam_p} ${samp}_MD.bam
    """
}