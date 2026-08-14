process XYRAT {
    cpus params.threads_xy
    input:
    tuple val(samp), val(dv_flag), path(mdbam)
    
    output:
    tuple val(samp), val(dv_flag), path(mdbam), path("${samp}_MD.bam.bai"), env(SEX_CALL), emit: sexed_bam
    tuple val(samp), env(SEX_CALL), env(RATIO), emit: sex_qc
    
    script:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    samtools index -b ${mdbam}
    samtools idxstats ${mdbam} | grep 'NC' > tmp.idxstat
    x=\$(grep '${params.x_chr}' tmp.idxstat | cut -f 3 || echo 0)
    y=\$(grep '${params.y_chr}' tmp.idxstat | cut -f 3 || echo 0)
    
    xd=\$(echo "scale=3; \${x}*150/139009144" | bc -l)
    yd=\$(echo "scale=3; \${y}*150/59476289" | bc -l)
    
    if (( \$(echo "\$yd > 0" | bc -l) )); then
        RATIO=\$(echo "scale=3; \${xd}/\${yd}" | bc -l)
    else
        if (( \$(echo "\$xd > 0" | bc -l) )); then
            RATIO="999"  # X coverage present, Y absent - likely female, assigning arbitrarily high ratio
        else
            RATIO="0"    # Both zero - potential QC issue
        fi
    fi

    if (( \$(echo "\$RATIO > ${params.x_thresh}" | bc -l) )); then
        SEX_CALL="FEMALE"
    elif (( \$(echo "\$RATIO < ${params.y_thresh} && \$RATIO > 0.5" | bc -l) )); then
        SEX_CALL="MALE"
    else
        SEX_CALL="UNDETERMINED" # includes between X and Y (ambiguous) and < 0.5 (likely qual issue)
    fi
    """

    stub:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}   

    samtools index -b ${mdbam}
    SEX_CALL="UNDETERMINED"
    RATIO="0"
    """
}