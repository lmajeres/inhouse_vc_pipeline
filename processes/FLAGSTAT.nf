process FLAGSTAT {
    cpus params.threads_fs
    input:
    tuple val(SEQ), val(seq), path(pair_bam), val(lane)

    output:
    tuple val(SEQ), path("${seq}_L${lane}.flagstat"), val(lane), emit: fs

    script:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    samtools flagstat -@ ${task.cpus} ${pair_bam} > ${seq}_L${lane}.flagstat
    """
}