/*///////////////////////////
//                         //
//   Cow VC Pipeline       //
//                         //
///////////////////////////*/

 /*=============================================
/   Global Variables, Libraries, and Scripts  /
============================================*/
// Libraries
import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import nextflow.io.ValueObject
import nextflow.util.KryoHelper

nextflow.preview.output = true

include { PETRIM ; BWA ; FLAGSTAT ; MERGEMD ; XYRAT ; DEEPVARIANT ; GLNEXUS ; PARSE_DVQC } from './processes.nf'

 /*=============================
/   Classes and Functions     /
============================*/

@ValueObject
class META {
    String seq_file        // original associated file(s)
    String sample_id       // given sample id
    String group_dir       // directory containing files
    Boolean dv_flag        // flag marking sample as needing more memory for DeepVariant
    Integer warn           // Count of warnings accrued; starts at 0 and counts up if WARN is triggered. Currently unused.

//    void flag() { this.warn += 1 }

//    String describe() {
//        return "META(seq_file=$seq_file, sample_id=$sample_id, group_dir=$group_dir, warn=$warn)"
//    }
}
KryoHelper.register(META)

def PARSE(csv) {
    def out = []
    new File(csv).eachLine { line, index ->
        if (index == 1) return // skip header
        def info = line.split(",").collect { it.trim() }
        out << new META(info[0], info[1], info[2], info[3].toLowerCase()=="true", 0)
    }
    return out
}

 /*================
/   Workflow     /
===============*/

workflow{
    main:
    if (!params.rawCsv && !params.bamCsv && !params.vcfCsv){
        error "Need at least one input CSV (rawCsv, bamCsv, or vcfCsv)."
    }

    rawInChannel = params.rawCsv ? Channel.fromList(PARSE(params.rawCsv)) : Channel.empty()
    bamInChannel = params.bamCsv ? Channel.fromList(PARSE(params.bamCsv)) : Channel.empty()
    vcfInChannel = params.vcfCsv ? Channel.fromList(PARSE(params.vcfCsv)) : Channel.empty()

    // Fetch lanes
    rawPairs = rawInChannel.map { META ->
            def r1s = file("${params.raws}/${META.group_dir}/${META.seq_file}_*_R1_001.fastq.gz")
            def lanes = r1s.collect { r1 ->
                def match = (r1.name =~ /_L00(\d+)_R1_/)
                def lane = match ? match[0][1] : "0" }
            return [META, r1s, lanes]
        }
        .transpose() // Flatten from run level to file level
        .map { META, r1, lane ->
            def r2 = file(r1.toString().replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'))
            def seq = META.seq_file
            return [META, seq, r1, r2, lane]
        }

    // Begin processing file pairs
    PETRIM(raw_pairs)
    align_in = PETRIM.out.trims.map { META, trim_p_fwd, trim_p_rev, trim_u_fwd, trim_u_rev, lane ->
        def seq = META.seq_file
        def samp = META.sample_id
        return [META, seq, samp, trim_p_fwd, trim_p_rev, trim_u_fwd, trim_u_rev, lane]
    }
    BWA(align_in)

    // QC collection for file pairs
    FLAGSTAT(BWA.out.fs_in)
    trim_qc = PETRIM.out.tstat.map { META, tstat_file, lane ->
        def stats = tstat_file.text
        def both_surv = (stats =~ /Both Surviving Reads: (\d+)/)[0][1] as Integer
        def both_surv_pct = (stats =~ /Both Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        def fwd_pct = (stats =~ /Forward Only Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        def rev_pct = (stats =~ /Reverse Only Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        // construct hashmap
        def key = "${META.sample_id}_${META.seq_file}_${lane}"
        def qc_data = [
            sample: META.sample_id,
            group: META.group_dir,
            run: META.seq_file,
            lane: lane,
            paired_survival_pct: both_surv_pct,
            r1_only_survival_pct: fwd_pct,
            r2_only_survival_pct: rev_pct,
            num_trim_pairs: both_surv
        ]
        return [key, qc_data]
    }
    flagstat_qc = FLAGSTAT.out.fs.map { META, flagstat_file, lane ->
        def num_aligned_reads = (flagstat_file.text =~ /(\d+) \+ \d+ primary mapped/)[0][1] as Integer
        def key = "${META.sample_id}_${META.seq_file}_${lane}"
        def qc_data = [num_aligned_reads: num_aligned_reads]
        return [key, qc_data]
    }
    file_qc = trim_qc.join(flagstat_qc).map { key, trim_data, flagstat_data ->
        def merged_data = trim_data + flagstat_data // merge hashmaps
        merged_data.aligned_pct = (merged_data.num_aligned_reads / 2) / merged_data.num_trim_pairs // calculate alignment rate
        return [key, merged_data] // new hashmap
    }

    // Now we bring things to sample-level
    bamIn2Merge = bamInChannel.map { META ->
        def samp = META.sample_id
        def bam_p = file("${params.bams}/${META.group_dir}/${META.seq_file}_*_P.bam")
        def bam_1u = file("${params.bams}/${META.group_dir}/${META.seq_file}_*_1U.bam")
        def bam_2u = file("${params.bams}/${META.group_dir}/${META.seq_file}_*_2U.bam")
        return [META, samp, bam_p, bam_1u, bam_2u]
    }
    merge_in = BWA.out.bams.map { META, bam_p, bam_1u, bam_2u -> 
            def samp = META.sample_id
            return [META, samp, bam_p, bam_1u, bam_2u]
        }
        .mix(bamIn2Merge)
        .groupTuple(by:1) // This is a bottleneck in the pipeline; it has to wait until all bams are done to be sure it got them all
        .map { META, samp, bam_p, bam_1u, bam_2u ->
            def dv_flag = META.any { it.dv_flag }
            def sort_p = bam_p.flatten()
                .sort { a,b -> a.name.compareTo(b.name) }
            def sort_1u = bam_1u.flatten()
                .sort { a,b -> a.name.compareTo(b.name) }
            def sort_2u = bam_2u.flatten
                .sort { a,b -> a.name.compareTo(b.name) }
            return [samp, dv_flag, sort_p, sort_1u, sort_2u]      
        }
    MERGEMD(merge_in)
    XYRAT(MERGEMD.out.mdbam)
    DEEPVARIANT(XYRAT.out.sexed_bam)
    gvcf2JointCall = vcfInChannel.map { META ->
            def gvcf_path = file("${params.vcfs}/${META.group_dir}/${META.seq_file}.g.vcf.gz")
            def gvcf_index_path = file("${params.vcfs}/${META.group_dir}/${META.seq_file}.g.vcf.gz.tbi")
            return [gvcf_path, gvcf_index_path]
        }
        .mix(DEEPVARIANT.out.gvcf)
    gvcf_manifest = gvcf2JointCall.map { gvcf_path, gvcf_index_path -> gvcf_path.toString() }
        .collectFile(name: 'gvcf_manifest.txt', newLine: true)
    GLNEXUS(gvcf_manifest)

    // QC collection for sample-level
    md_qc = MERGEMD.out.mdstat.map { samp, mdstat_file ->
        def stats = mdstat_file.text
        def exam = (stats =~ /EXAMINED: (\d+)/)[0][1] as Integer
        def cover = (exam*150.0)/params.genome_size
        def dups = (stats =~ /DUPLICATE TOTAL: (\d+)/)[0][1] as Integer
        def dup_rate = dups/exam
        def lib_size = (stats =~ /ESTIMATED_LIBRARY_SIZE: (\d+)/)[0][1] as Long
        // construct hashmap
        def qc_data = [
            sample: samp,
            coverage: cover,
            dup_pct: dup_rate,
            est_lib_size: lib_size
        ]
        return [samp, qc_data]
    }
    xy_qc = XYRAT.out.sex_qc.map { samp, sex_call, ratio ->
        def qc_data = [
            sex_called: sex_call,
            xy_ratio: ratio
        ]
        return [samp, qc_data]
    }
    samp1_qc = md_qc.join(xy_qc).map { samp, md_data, xy_data ->
        def merge1_data = md_data + xy_data // merge hashmaps
        return [samp, merge1_data] // new hashmap
    }
    PARSE_DVQC(DEEPVARIANT.out.vcf_stat)
    dv_qc = PARSE_DVQC.out.json.map { samp, json_file ->
        def stats_text = json_file.text
        def stats = new groovy.json.JsonSlurper().parseText(stats_text)
        // construct hashmap
        def qc_data = [
            total_variants: stats.total_variants,
            refcall_variants: stats.refcall_variants,
            nonrefcall_variants: stats.nonrefcall_variants,
            pct_refcall: stats.pct_refcall,
            titv_ratio: stats.titv_ratio,
            mean_gq: stats.mean_gq,
            median_gq: stats.median_gq,
            pct_gq_gt10: stats.pct_gq_gt10,
            pct_gq_gt20: stats.pct_gq_gt20,
            pct_gq_gt30: stats.pct_gq_gt30
        ]
        return [samp, qc_data]
    }
    samp2_qc = samp1_qc.join(dv_qc).map { samp, merge1_data, dv_data ->
        def merge2_data = merge1_data + dv_data
        return [samp, merge2_data]
    }

    publish:
    unmerged_bams = BWA.out.bams // [metadata, [path_P, path_1U, path_2U]]
    final_bams = XYRAT.out.sexed_bam // [sample_id, path(mdbam), path(bai), sex call]
    gvcfs = DEEPVARIANT.out.gvcf // [path(gvcf)]
    bcf = GLNEXUS.out.bcf // [path(bcf)] (singular)
    run_qc = file_qc.map { key, qc -> return qc} // hashmap indexed by run_lane_sample
    samp_qc = samp2_qc.map { samp, qc -> return qc} // hashmap indexed by sample
}

output{
    unmerged_bams {
        path "unmerged_bams"
	    mode 'symlink'
    }
    final_bams {
        path "md_bams"
	    mode 'symlink'
    }
    gvcfs {
        path "gvcfs"
	    mode 'symlink'
    }
    bcf {
        path "."
	    mode 'move'
    }
    run_qc {
	    path "."
	    mode 'move'
	    index {
		    path "${params.outName}_runqc.csv"
		    header true
	    }
    }
    samp_qc {
        path "."
	    mode 'move'
        index {
		    path "${params.outName}_sampqc.csv"
		    header true
            }
    }
}
