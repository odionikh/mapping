#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { BWA_INDEX } from './modules/bwa'
include { BWA_MEM } from './modules/bwa'
include { SAMTOOLS_SORT } from './modules/samtools.nf'
include { SAMTOOLS_INDEX } from './modules/samtools.nf'
include { SAMTOOLS_FAIDX } from './modules/samtools.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools.nf'
include { SAMTOOLS_IDXSTATS } from './modules/samtools.nf'
include { SAMTOOLS_DEPTH } from './modules/samtools.nf'
include { BCFTOOLS_MPILEUP } from './modules/bcftools.nf'
include { BCFTOOLS_VIEW } from './modules/bcftools.nf'


// Define channels from input CSV and workflows

workflow {
    input_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> 
            def sample_id = row.sample_id
            def reference = file(row.reference)
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = file(row.fastq_2)
            [sample_id, reference, fastq_1, fastq_2]
        }

    // Create reference_ch with unique references
    reference_ch = input_ch
        .map { it[1] }
        .unique()

    // Index reference genomes
    BWA_INDEX(reference_ch)

    // Combine samples with their corresponding indexed reference
    aligned_input_ch = input_ch
        .map { sample_id, reference, fastq_1, fastq_2 -> 
            [reference.simpleName, sample_id, fastq_1, fastq_2]
        }
        .combine(BWA_INDEX.out.index, by: 0)
        .map { ref_name, sample_id, fastq_1, fastq_2, reference, index_files ->
            [sample_id, fastq_1, fastq_2, ref_name, reference, index_files]
        }

    // Perform alignment

    bwa_mem_ch = BWA_MEM(aligned_input_ch)

    //Run SAMTOOLSORT
    samtools_sort_ch = SAMTOOLS_SORT(bwa_mem_ch.sam)

    //Run SAMTOOLS_INDEX
    samtools_index_ch = SAMTOOLS_INDEX(samtools_sort_ch.bam)

    //Run SAMTOOLS_FAIDX
    SAMTOOLS_FAIDX(reference_ch)

    //Run SAMTOOLS_FLAGSTAT
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.bam)

    //Run SAMTOOLS_IDXSTATS
    SAMTOOLS_IDXSTATS(SAMTOOLS_INDEX.out.bam_bai)

    //Run SAMTOOLS_DEPTH
    SAMTOOLS_DEPTH(SAMTOOLS_SORT.out.bam)

    //Input for BCFTOOLS
    mpileup_input_ch = samtools_sort_ch.bam
        .join(input_ch.map {sample_id, reference, fastq_1, fastq_2 -> [sample_id, reference] })
        .map { sample_id, bam, reference  ->
            [ sample_id, bam, reference ]
        }

    // Run BCFTOOLS_MPILEUP
    bcftools_mpileup_ch = BCFTOOLS_MPILEUP(mpileup_input_ch)

    // Run BCFTOOLS_VIEW
    BCFTOOLS_VIEW(bcftools_mpileup_ch.bcf)

}
