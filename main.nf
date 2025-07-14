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

include { PETRIM; BWA; FLAGSTAT } from "processes.nf"
params.raws="/lustre/isaac24/proj/UTK0204/lmajeres/ihvc_cattle/raw"

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
    SEQlist = PARSE(${params.input_csv}) // Needs a csv of files we want to parse! format as file_basename, sample_id, dir_of_group

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
        def seq = SEQ.raw_sed_id
        def samp = SEQ.sample_id
        return [SEQ, seq, samp, trim_pair, trim_fwd, trim_rev, lane]
    }
    BWA(align_in)

    // QC collection for file pairs
    FLAGSTAT(BWA.out.fs_in)
    trim_qc = PETRIM.out.tstat.map { SEQ, tstat_file, lane ->
        def stats = tstat_file.text // get file
        // extract info
        def both_surv = (stats =~ /Both Surviving Reads: (\d+)/)[0][1] as Integer
        def both_surv_pct = (stats =~ /Both Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        def fwd_pct = (stats =~ /Forward Only Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        def rev_pct = (stats =~ /Reverse Only Surviving Read Percent: ([\d.]+)/)[0][1] as Double
        // construct hashmap
        def key = "${SEQ.sample_id}_${SEQ.seq}_${lane}"
        def qc_data = [
            sample: SEQ.sample_id,
            group: SEQ.group,
            run: SEQ.seq,
            lane: lane,
            paired_survival_pct: both_surv_pct,
            r1_only_survival_pct: fwd_pct,
            r2_only_survival_pct: rev_pct,
            num_trim_pairs: both_surv
        ]
        return [key, qc_data]
    }
    flagstat_qc = FLAGSTAT.out.stats.map { SEQ, flagstat_file, lane ->
        def stats = flagstat_file.text
        def num_aligned_reads = (stats =~ /(\d+) \+ \d+ primary mapped/)[0][1] as Integer
        def key = "${SEQ.sample_id}_${SEQ.seq}_${lane}"
        def qc_data = [num_aligned_reads: num_aligned_reads]
        return [key, qc_data]
    }
    combined_qc = trim_qc.join(flagstat_qc).map { key, trim_data, flagstat_data ->
        def merged_data = trim_data + flagstat_data // merge hashmaps
        merged_data.aligned_pct = (merged_data.num_aligned_reads / 2) / merged_data.num_trim_pairs // calculate alignment rate
        return [key, merged_data] // new hashmap
    }
    combined_qc.collectFile(
        name: 'trim_align_qc.csv',
        newLine: true,
        storeDir: params.outdir,
        sort: true
    ) { key, data ->
        "${data.sample},${data.run},${data.lane},${data.paired_survival_pct},${data.r1_only_survival_pct},${data.r2_only_survival_pct},${data.num_trim_pairs},${data.aligned_pct},${data.num_aligned_reads}"
    }

    // Now we bring things to sample-level

    // Outputs
    publish:

}

output{

}