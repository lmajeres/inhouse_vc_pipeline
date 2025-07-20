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

include { PETRIM ; BWA ; FLAGSTAT ; MERGEMD ; XYRAT ; DEEPVARIANT ; GLNEXUS ; PARSE_DVQC } from './processes.nf'

 /*=============================
/   Classes and Functions     /
============================*/

class SEQ {
    String raw_seq_id      // original file name
    String sample_id       // given sample id
    String group           // directory containing files
    Integer warn           // Count of warnings accrued; starts at 0 and counts up if WARN is triggered. Currently unused.

    SEQ(String raw_seq_id, String sample_id, String group) {
        this.raw_seq_id = raw_seq_id
        this.sample_id  = sample_id
        this.group = group
        this.warn = 0
    }

    void flag() { this.warn += 1 }

    String describe() {
        return "SEQ(raw_seq_id=$raw_seq_id, sample_id=$sample_id, group=$group, warn=$warn)"
    }
}

def PARSE(csv) {
    def out = []
    new File(csv).eachLine { line, index ->
        if (index == 1) return // skip header
        def info = line.split(",").collect { it.trim() }
        out << new SEQ(info[0], info[1], info[2])
    }
    return out
}

 /*================
/   Workflow     /
===============*/

workflow{
    main:
    // Starts at run-level (multiple possible sequencing runs/libraries per sample, but doesn't consider technical reps (ie lanes))
    SEQlist = PARSE("${params.inCsv}") // Needs a csv of files we want to parse! format as file_basename, sample_id, dir_of_group

    // Fetch lanes
    raw_pairs = Channel.fromList(SEQlist)
        .map { SEQ ->
            def r1s = file("${params.raws}/${SEQ.group}/${SEQ.raw_seq_id}_*_R1_001.fastq.gz")
            def lanes = r1s.collect { r1 ->
                def match = (r1.name =~ /_L00(\d+)_R1_/)
                def lane = match ? match[0][1] : "no_lane" }
            return [SEQ, r1s, lanes]
        }
        .transpose() // Flatten from run level to file level
        .map { SEQ, r1, lane ->
            def r2 = file(r1.toString().replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'))
            def seq = SEQ.raw_seq_id
            return [SEQ, seq, [r1, r2], lane]
        }

    // Begin processing file pairs
    PETRIM(raw_pairs)
    align_in = PETRIM.out.trims.map { SEQ, trim_pair, trim_fwd, trim_rev, lane ->
        def seq = SEQ.raw_seq_id
        def samp = SEQ.sample_id
        return [SEQ, seq, samp, trim_pair, trim_fwd, trim_rev, lane]
    }
    BWA(align_in)

    // QC collection for file pairs ** TODO: Test if this works
    FLAGSTAT(BWA.out.fs_in)
    trim_qc = PETRIM.out.tstat.map { SEQ, tstat_file, lane ->
        def stats = tstat_file.text // get file
        // extract info
        def both_surv = (stats =~ /Both Surviving Reads: (\d+)/)[0][1] as Integer
        def both_surv_pct = (stats =~ /Both Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        def fwd_pct = (stats =~ /Forward Only Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        def rev_pct = (stats =~ /Reverse Only Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        // construct hashmap
        def key = "${SEQ.sample_id}_${SEQ.raw_seq_id}_${lane}"
        def qc_data = [
            sample: SEQ.sample_id,
            group: SEQ.group,
            run: SEQ.raw_seq_id,
            lane: lane,
            paired_survival_pct: both_surv_pct,
            r1_only_survival_pct: fwd_pct,
            r2_only_survival_pct: rev_pct,
            num_trim_pairs: both_surv
        ]
        return [key, qc_data]
    }
    flagstat_qc = FLAGSTAT.out.fs.map { SEQ, flagstat_file, lane ->
        def num_aligned_reads = (flagstat_file.text =~ /(\d+) \+ \d+ primary mapped/)[0][1] as Integer
        def key = "${SEQ.sample_id}_${SEQ.raw_seq_id}_${lane}"
        def qc_data = [num_aligned_reads: num_aligned_reads]
        return [key, qc_data]
    }
    file_qc = trim_qc.join(flagstat_qc).map { key, trim_data, flagstat_data ->
        def merged_data = trim_data + flagstat_data // merge hashmaps
        merged_data.aligned_pct = (merged_data.num_aligned_reads / 2) / merged_data.num_trim_pairs // calculate alignment rate
        return [key, merged_data] // new hashmap
    }

    // Now we bring things to sample-level
    merge_in = BWA.out.bams.map { SEQ, bam_paths -> 
            def samp = SEQ.sample_id
            return [SEQ, samp, bam_paths]
        }
        .groupTuple(by:1) // This will be a bottleneck in the pipeline, since it will have to wait until all bams are done to be sure it got them all
    MERGEMD(merge_in)
    XYRAT(MERGEMD.out.mdbam)
    DEEPVARIANT(XYRAT.out.sexed_bam)
    gvcf_manifest = DEEPVARIANT.out.gvcf.collectFile(name: 'gvcf_manifest.txt', newLine: true)
    GLNEXUS(gvcf_manifest)

    // QC collection for sample-level
    md_qc = MERGEMD.out.mdstat.map { samp, mdstat_file ->
        def stats = mdstat_file.text // get file
        // extract info
        def exam = (stats =~ /EXAMINED: (\d+)/)[0][1] as Integer
        def cover = (exam*150)/2770669782
        def dups = (stats =~ /DUPLICATE TOTAL: (\d+)/)[0][1] as Integer
        def dup_rate = dups/exam
        def lib_size = (stats =~ /ESTIMATED_LIBRARY_SIZE: (\d+)/)[0][1] as Integer
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
    // Outputs ** PUBLISH: Run-level bams (in-case there is an issue, since we aren't doing auto-qc), sample-level MD bams, gvcfs, final bcf, qc csvs
    // run bams, sample bams, and gvcfs should perhaps live in scratch? Maybe ask Trevor about if I can keep the gvcfs somewhere more permanent
    // final bcf and csvs can be returned in main directory.
    // instead of doing collectfile, make an index file from the hashmap?
    publish:
    unmerged_bams = BWA.out.bams // [metadata, [path_P, path_1U, path_2U]]
    final_bams = MERGEMD.out.mdbam // [sample_id, path(mdbam)]
    gvcfs = DEEPVARIANT.out.gvcf // [path(gvcf)]
    bcf = GLNEXUS.out.bcf // [path(bcf)] (singular)
    run_qc = file_qc.map { key, qc -> return qc} // hashmap indexed by run, lane, sample
    samp_qc = samp2_qc.map { samp, qc -> return qc} // hashmap indexed by sample
}

output{
    unmerged_bams {
        path "${params.outputScratch}/unmerged_bams"
    }
    final_bams {
        path "${params.outputScratch}/md_bams"
    }
    gvcfs {
        path "${params.outputScratch}/gvcfs"
    }
    bcf {
        path "${outputDir}"
    }
    run_qc {
       path "${outputDir}"
       index { 
        path "${params.outName}_runqc.csv"
        header true
       }
    }
    samp_qc {
        path "${outputDir}"
        index { 
            path "${params.outName}_sampqc.csv"
            header true
       }
    }
}