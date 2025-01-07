process BWA_INDEX {
    tag "$reference.simpleName"

    publishDir "${params.outdir}/bwa_index/", mode: 'copy'
    
    input:
    path reference

    output:
    tuple val(reference.simpleName), path(reference), path("${reference}.*"), emit: index

    script:
    """
    echo "Starting BWA_INDEX for $reference"
    bwa index ${reference}
    """
}

process BWA_MEM {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_1), path(fastq_2), val(ref_name), path(reference), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}_${reference.simpleName}.sam"), emit: sam

    script:
    """
    echo "Starting BWA_MEM for sample: $sample_id"
    echo "Reference: $reference"
    echo "FASTQ1: $fastq_1"
    echo "FASTQ2: $fastq_2"
    bwa mem -t ${task.cpus} ${reference} ${fastq_1} ${fastq_2} > ${sample_id}_${reference.simpleName}.sam
    """
}

