process BCFTOOLS_MPILEUP {
    tag "$sample_id"

    publishDir "${params.outdir}/$sample_id/", mode: 'copy'
    input:
    tuple val(sample_id),  path(bam), path(reference)  // tuple val(sample_id), path(bam), path(bam_bai), path(reference), path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}_${reference.simpleName}.mpileup.bcf"), emit: bcf

    script:
    """
    bcftools mpileup -Ou -f ${reference} ${bam} | \
    bcftools call -mv -Ob -o ${sample_id}_${reference.simpleName}.mpileup.bcf
    """
}


process BCFTOOLS_VIEW {
    tag "$sample_id"

    publishDir "${params.outdir}/$sample_id/", mode: 'copy'

    input:
    tuple val(sample_id), path(bcf_file)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf"), emit: filtered_vcf

    script:
    """
    bcftools view -i 'QUAL>=10' ${bcf_file} -Ov -o ${sample_id}.filtered.vcf
    """
}

