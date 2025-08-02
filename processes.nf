process PETRIM {
    cpus = params.threads_trim
    maxForks = params.maxfork_trim
    input:
    tuple val(SEQ), val(seq), path(fastq_r1), path(fastq_r2), val(lane)
        
    output:
    tuple val(SEQ),path("${seq}_L${lane}_t_1P.fastq.gz"), path("${seq}_L${lane}_t_2P.fastq.gz"),
        path("${seq}_L${lane}_t_1U.fastq.gz"), path("${seq}_L${lane}_t_2U.fastq.gz"), val(lane), emit: trims
    tuple val(SEQ), path("${seq}_L${lane}.tstat"), val(lane), emit: tstat
        
    script:
    """
    trimmomatic PE \
    -threads ${task.cpus} \
    -summary ${seq}_L${lane}.tstat \
    ${fastq_r1} ${fastq_r2} \
    -baseout ${seq}_L${lane}_t.fastq.gz \
    HEADCROP:1 ILLUMINACLIP:${params.adapters}:2:30:10 \
    LEADING:20 SLIDINGWINDOW:15:20 MINLEN:75
    """

    stub:
    """
    zcat ${fastq_r1} | head -n 40 > fwd_head.fastq
    zcat ${fastq_r2} | head -n 40 > rev_head.fastq
    trimmomatic PE \
    -threads 1 \
    -summary ${seq}_L${lane}.tstat \
    fwd_head.fastq rev_head.fastq \
    -baseout ${seq}_L${lane}_t.fastq \
    HEADCROP:1 ILLUMINACLIP:${params.adapters}:2:30:10 \
    LEADING:20 SLIDINGWINDOW:15:20 MINLEN:75
    touch ${seq}_L${lane}_t_1U.fastq
    touch ${seq}_L${lane}_t_2U.fastq
    """
}

process BWA {
    cpus = params.threads_bwa
    memory = params.mem_bwa
    maxForks = params.maxfork_bwa
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

process FLAGSTAT {
    cpus = params.threads_fs
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

process MERGEMD {
    cpus = params.threads_merge
    input:
    tuple val(samp), path(bam_p), path(bam_1u), path(bam_2u)

    output:
    tuple val(samp), path("${samp}_MD.bam"), emit: mdbam
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

process XYRAT {
    cpus = params.threads_xy
    input:
    tuple val(samp), path(mdbam)
    
    output:
    tuple val(samp), path(mdbam), path("${samp}_MD.bam.bai"), env(SEX_CALL), emit: sexed_bam
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

process DEEPVARIANT {
    cpus = params.threads_dv
    memory = params.mem_dv
    maxForks = params.maxfork_dv
    input:
    tuple val(samp), path(mdbam), path(index), val(sex)

    output:
    tuple path("${samp}.g.vcf.gz"), path("${samp}.g.vcf.gz.tbi"), emit: gvcf
    tuple val(samp), path("${samp}.visual_report.html"), emit: vcf_stat

    script:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    if [ "${sex}" == "MALE" ]; then
        REGIONS="${params.autosomes} ${params.x_chr} ${params.y_chr}"
        HAPLOID_FLAG="--haploid_contigs=${params.x_chr},${params.y_chr}"
    elif [ "${sex}" == "FEMALE" ]; then
        REGIONS="${params.autosomes} ${params.x_chr}"
        HAPLOID_FLAG=""
    else
        REGIONS="${params.autosomes}"
        HAPLOID_FLAG=""
    fi

    singularity exec --cleanenv ~/tools/deepvariant_latest.sif \
    run_deepvariant \
    --model_type=WGS \
    --vcf_stats_report=true \
    --postprocess_cpus=0 \
    --make_examples_extra_args='small_model_call_multiallelics=false' \
    --ref=${params.dvRef} \
    --regions="\${REGIONS}" \
    --reads=${mdbam} \
    --output_vcf=${samp}.vcf.gz \
    --output_gvcf=${samp}.g.vcf.gz \
    --num_shards=${task.cpus} \
    \$HAPLOID_FLAG
    """

    stub:
    """
    singularity exec --cleanenv ~/tools/deepvariant_latest.sif \
    run_deepvariant \
    --model_type=WGS \
    --vcf_stats_report=true \
    --postprocess_cpus=0 \
    --make_examples_extra_args='small_model_call_multiallelics=false' \
    --ref=${params.dvRef} \
    --regions="${params.autosomes}" \
    --reads=${mdbam} \
    --output_vcf=${samp}.vcf.gz \
    --output_gvcf=${samp}.g.vcf.gz \
    --num_shards=${task.cpus} \
    --dry_run=true
    touch ${samp}.g.vcf.gz
    touch ${samp}.g.vcf.gz.tbi
    touch ${samp}.visual_report.html
    """
}

process GLNEXUS {
    cpus = params.threads_glnex
    memory = params.mem_glnex
    clusterOptions = '--nodelist=clrv1203'
    input:
    path(manifest)

    output:
    path("${params.outName}.bcf"), emit: bcf

    script:
    """
    export TMPDIR=\${PWD}/tmp_dir
    mkdir -p \${TMPDIR}

    singularity exec --cleanenv ~/tools/glnexus_v1.4.1.sif \
    glnexus_cli \
    -t ${task.cpus} \
    -m ${params.budg_glnex} \
    -a -c DeepVariant \
    --list ${manifest} > ${params.outName}.bcf
    """

    stub:
    """
    touch ${params.outName}.bcf
    """
}

process PARSE_DVQC {
    cpus = params.threads_dvqc
    input:
    tuple val(samp), path(html)

    output:
    tuple val(samp), path("${samp}.json"), emit: json

    script:
    """
    python ${params.dvqcParser} --input ${html} --output ${samp}.json
    """

    stub:
    """
    echo "{
    \"total_variants\": 0,
    \"refcall_variants\": 0,
    \"nonrefcall_variants\": 0,
    \"pct_refcall\": 0,
    \"titv_ratio\": 0,
    \"mean_gq\": 0,
    \"median_gq\": 0,
    \"pct_gq_gt10\": 0,
    \"pct_gq_gt20\": 0,
    \"pct_gq_gt30\": 0
}" > ${samp}.json
    """
}
