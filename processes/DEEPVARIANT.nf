process DEEPVARIANT {
    cpus params.threads_dv
    memory { dv_flag ? params.mem_dv_flag : params.mem_dv }
    maxForks params.maxfork_dv
    input:
    tuple val(samp), val(dv_flag), path(mdbam), path(index), val(sex)

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
    --intermediate_results_dir \${TMPDIR} \
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