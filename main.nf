nextflow.enable.dsl = 2

params.sampleID = ""
params.fastqR1 = ""
params.fastqR2 = ""
params.ref = ""
params.outDir = ""

sampleID = params.sampleID
fastqR1 = Channel.fromPath(params.fastqR1)
fastqR2 = Channel.fromPath(params.fastqR2)
ref = Channel.fromPath(params.ref)

workflow {
    main:
        mapping( sampleID, fastqR1, fastqR2, ref )
        markdup( sampleID, mapping.out.bam )
        bamIndex( mapping.out.bam )
    emit:
        mapping.out.bam
        bamIndex.out.bai
}

process mapping {
    machineType "mem1_ssd1_v2_x36"
    container "quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path fastqR1
        path fastqR2
        path ref
    output:
        path "${sampleID}.bam", emit: bam
        path "${sampleID}.bam.bai", emit: bai
    script:
        """
        mkdir -p genome
        zcat "${ref}" | tar xvf - -C genome
        genome_file=`ls ./genome/*.bwt`
        genome_file="\${genome_file%.bwt}"
        bwa mem -K 100000000 -t 36 -M \
            \${genome_file} \
            ${fastqR1} ${fastqR2} \
            -R "@RG\\tID:None\\tPL:None\\tPU:None\\tLB:None\\tSM:${sampleID}" | \
            samtools sort --threads 36 - > ${sampleID}.bam
        samtools index ${sampleID}.bam
        """
}

process markdup {
    machineType "mem1_ssd1_v2_x36"
    container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    publishDir "${params.outDir}"
    input:
        val sampleID
        path bam
    output:
        path "./${sampleID}.md.bam", emit: bam_md
        path "./${sampleID}.md.bam.bai", emit: bai_md
    script:
        """
        gatk --java-options "-Xmx64g" \
            MarkDuplicates \
                --INPUT ${bam} \
                --METRICS_FILE ${sampleID}.md.bam.metrics \
                --TMP_DIR . \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_INDEX true \
                --OUTPUT ${sampleID}.md.bam
        mv ${sampleID}.md.bai ${sampleID}.md.bam.bai
        """
}

process bamIndex {
    machineType "mem1_ssd1_v2_x2"
    container "quay.io/biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:c17ce694dd57ab0ac1a2b86bb214e65fedef760e-0"
    publishDir "${params.outDir}"
    input:
        path bam
    output:
        path "bam/*.bai", emit: bai
    script:
        """
        mkdir -p bam
        cp ${bam} ./bam/
        bam=`ls ./bam/*.bam`
        samtools index \${bam}
        """
}

