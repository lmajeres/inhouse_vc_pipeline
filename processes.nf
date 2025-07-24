process PETRIM {
    cpus = params.threads_trim
    input:
    tuple val(SEQ), val(seq), path(fastqs), val(lane)
        
    output:
    tuple val(SEQ), path("*_t_*P.fastq.gz"), path("*_t_1U.fastq.gz"), path("*_t_2U.fastq.gz"), val(lane), emit: trims
    tuple val(SEQ), path("*.tstat"), val(lane), emit: tstat
        
    script:
    """
    trimmomatic PE \
    -threads ${task.cpus} \
    -summary ${seq}_L${lane}.tstat \
    ${fastqs} \
    -baseout ${seq}_L${lane}_t.fastq.gz \
    HEADCROP:1 ILLUMINACLIP:${params.adapters}:2:30:10 \
    LEADING:20 SLIDINGWINDOW:15:20 MINLEN:75
    """
}

process BWA {
    cpus = params.threads_bwa
    memory = params.mem_bwa
    input:
    tuple val(SEQ), val(seq), val(samp), path(trim_pair), path(trim_fwd), path(trim_rev), val(lane)

    output:
    tuple val(SEQ), path("*P.bam"), path("*1U.bam"), path("*2U.bam"), emit: bams
    tuple val(SEQ), val(seq), path("*P.bam"), val(lane), emit: fs_in
    
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
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_P.bam

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_fwd} |
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_1U.bam

    bwa-mem2 mem \
    -R "\${RG}" \
    -t ${task.cpus} \
    ${params.bwaIndex} \
    ${trim_rev} |
    samtools sort -@ ${task.cpus} -o ${seq}_L${lane}_2U.bam
    """
}

process FLAGSTAT {
    cpus = params.threads_fs
    input:
    tuple val(SEQ), val(seq), path(pair_bam), val(lane)

    output:
    tuple val(SEQ), path("*.flagstat"), val(lane), emit: fs

    script:
    """
    samtools flagstat -@ ${task.cpus} ${pair_bam} > ${seq}_L${lane}.flagstat
    """
}

process MERGEMD {
    cpus = params.threads_merge
    input:
    tuple val(SEQ), val(samp), path(bam_p), path(bam_1u), path(bam_2u)

    output:
    tuple val(samp), path("*MD.bam"), emit: mdbam
    tuple val(samp), path("*MD.stat"), emit: mdstat

    script:
    """
    samtools merge -@ ${task.cpus} -u -o - ${bam_p} ${bam_1u} ${bam_2u} |
    samtools markdup -@ ${task.cpus} -S -d 2500 \
    -s -f ${samp}_MD.stat - ${samp}_MD.bam
    """
}

process XYRAT {
    cpus = params.threads_xy
    input:
    tuple val(samp), path(mdbam)
    
    output:
    tuple val(samp), path(mdbam), path("*.bai"), env(SEX_CALL), emit: sexed_bam
    tuple val(samp), env(SEX_CALL), env(RATIO), emit: sex_qc
    
    script:
    """
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
}

process DEEPVARIANT {
    cpus = params.threads_dv
    memory = params.mem_dv
    input:
    tuple val(samp), path(mdbam), path(index), val(sex)

    output:
    path("*.g.vcf.gz"), emit: gvcf
    tuple val(samp), path("*.html"), emit: vcf_stat

    script:
    """
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
    --regions=\${REGIONS} \
    --reads=${mdbam} \
    --output_vcf=${samp}.vcf.gz \
    --output_gvcf=${samp}.g.vcf.gz \
    --num_shards=${task.cpus} \
    \$HAPLOID_FLAG
    """
}

process GLNEXUS {
    cpus = params.threads_glnex
    mem = params.mem_glnex
    input:
    path(manifest)

    output:
    path("*.bcf"), emit: bcf

    script:
    """
    singularity exec --cleanenv ~/tools/glnexus_v1.4.1.sif \
    glnexus_cli \
    -t ${task.cpus} \
    -m ${params.budg_glnex} \
    -l ${manifest} \
    -a -c DeepVariant > ${params.outName}.bcf
    """
}

process PARSE_DVQC {
    cpus = params.threads_dvqc
    input:
    tuple val(samp), path(html)

    output:
    tuple val(samp), path("*.json"), emit: json

    script:
    """
    python ${params.dvqcParser} --input ${html} --output ${samp}.json
    """
}
