process GLNEXUS {
    cpus params.threads_glnex
    memory params.mem_glnex
//  clusterOptions = '--nodelist=clrv1203' // This doesn't work
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