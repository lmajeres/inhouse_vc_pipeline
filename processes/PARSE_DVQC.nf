process PARSE_DVQC {
    cpus params.threads_dvqc
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