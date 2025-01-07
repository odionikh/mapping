process SAMTOOLS_SORT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${sam_file}
    """
}

process SAMTOOLS_INDEX {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam_bai

    script:
    """
    samtools index ${bam}
    """
}

process SAMTOOLS_FAIDX {
    tag "$reference.simpleName"

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.fai"), emit: faidx

    script:
    """
    samtools faidx ${reference}
    """
}

process SAMTOOLS_FLAGSTAT {
    tag "$sample_id"

    publishDir "${params.outdir}/$sample_id/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.flagstat"), emit: flagstat

    script:
    """
    samtools flagstat $bam > ${sample_id}.flagstat
    """
}

process SAMTOOLS_IDXSTATS {
    tag "$sample_id"

    publishDir "${params.outdir}/$sample_id/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bam_bai)

    output:
    tuple val(sample_id), path("${sample_id}.idxstats"), emit: idxstats

    script:
    """
    samtools idxstats $bam > ${sample_id}.idxstats
    """
}

process SAMTOOLS_DEPTH {
    tag "$sample_id"

    publishDir "${params.outdir}/$sample_id/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.depth"), emit: depth

    script:
    """
    samtools depth -aa $bam > ${sample_id}.depth
    """
}

