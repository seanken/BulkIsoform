version 1.0

import "BulkPipe.wdl" as BulkPipe
import "MakeRef.wdl" as MakeRef

workflow Pipeline {
    input {
        # MakeRef inputs
        File? fa
        File? gtf_for_ref
        String? ref_outdir
        
        # BulkPipe inputs
        Array[File]? fq1
        Array[File]? fq2
        String? pipeline_outdir
        String? ref_comb
        String? isoMethod
        Int? s
        String? gene_name
        Int? numGibbsSamples
    }
    
    # Determine which workflow to run based on inputs
    Boolean run_makeref = defined(fa) && defined(gtf_for_ref)
    Boolean run_bulkpipe = defined(fq1) && defined(fq2)
    Boolean both_workflows = run_makeref && run_bulkpipe && !defined(ref_comb)
    
    if (run_makeref) {
        call MakeRef.MakeRef {
            input:
                fa = select_first([fa]),
                gtf = select_first([gtf_for_ref]),
                outdir = select_first([ref_outdir, "ref"])
        }
    }
    
    # Determine ref_comb: use MakeRef output if both workflows run, otherwise use provided ref_comb
    String final_ref_comb = if both_workflows then select_first([MakeRef.output_directory]) else select_first([ref_comb, ""])
    
    if (run_bulkpipe) {
        call BulkPipe.BulkPipeline {
            input:
                fq1 = select_first([fq1]),
                fq2 = select_first([fq2]),
                outdir = select_first([pipeline_outdir, "results"]),
                ref_comb = final_ref_comb,
                isoMethod = select_first([isoMethod, "Salmon"]),
                s = select_first([s, 2]),
                gene_name = select_first([gene_name, "gene_name"]),
                numGibbsSamples = select_first([numGibbsSamples, 0])
        }
    }
    
    output {
        # MakeRef outputs
        File? ref_directory = MakeRef.output_directory
        
        # BulkPipe outputs
        File? mapped_bam_out = BulkPipeline.mapped_bam_out
        File? mapper_qc = BulkPipeline.mapper_qc
        File? exon_counts = BulkPipeline.exon_counts
        File? exon_counts_sum = BulkPipeline.exon_counts_sum
        File? with_intron_counts = BulkPipeline.with_intron_counts
        File? with_intron_counts_sum = BulkPipeline.with_intron_counts_sum
        File? picard_qc = BulkPipeline.picard_qc
        File? picard_dup = BulkPipeline.picard_dup
        File? rsem_out = BulkPipeline.rsem_out
        File? salmon_out = BulkPipeline.salmon_out
    }
}
