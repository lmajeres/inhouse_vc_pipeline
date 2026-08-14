process PETRIM {
    cpus params.threads_trim
    maxForks params.maxfork_trim
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