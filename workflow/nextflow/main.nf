// Nextflow version of the CBIO 7501 capstone: steps 2-6, the same as workflow/Snakefile.
// Reads the test data that step 1 fetched into data/raw/test.
// Run from the project folder: pixi run run-nextflow

params.data   = "data/raw/test"
params.region = "chr22:34127223-34838003"
params.tumour = "HG008-T"
params.normal = "HG008-N-D"
params.outdir = "results/capstone_nextflow"

// Step 2: prepare the reference piece (renamed: GATK can misread a colon)
process PREPARE_REFERENCE {
    input:
    path fasta

    output:
    tuple path("ref.fa"), path("ref.fa.*"), path("ref.dict")

    script:
    def name = params.region.replace(":", "_").replace("-", "_")
    """
    sed 's/^>.*/>${name}/' ${fasta} > ref.fa
    samtools faidx ref.fa
    gatk CreateSequenceDictionary -R ref.fa -O ref.dict
    bwa index ref.fa
    """
}

// Step 3: align one sample's reads (bwa mem, with the sample name as read group)
process ALIGN {
    tag "${sample}"
    cpus 4

    input:
    tuple val(sample), path(r1), path(r2)
    tuple path(ref), path(ref_index), path(dict)

    output:
    tuple val(sample), path("${sample}.unsorted.bam")

    script:
    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' ${ref} ${r1} ${r2} \\
      | samtools view -b -o ${sample}.unsorted.bam -
    """
}

// Step 4: mark duplicates, sort by position, index
process MARKDUP {
    tag "${sample}"
    cpus 4

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script:
    """
    samtools fixmate -m ${bam} - | samtools sort -@ ${task.cpus} - | samtools markdup - ${sample}.bam
    samtools index ${sample}.bam
    """
}

// Step 5a: Mutect2, tumour against normal
process CALL {
    input:
    tuple path(ref), path(ref_index), path(dict)
    tuple val(t), path(tumour_bam), path(tumour_bai)
    tuple val(n), path(normal_bam), path(normal_bai)

    output:
    tuple path("unfiltered.vcf.gz"), path("unfiltered.vcf.gz.tbi"), path("unfiltered.vcf.gz.stats")

    script:
    """
    gatk Mutect2 -R ${ref} -I ${tumour_bam} -I ${normal_bam} -normal ${params.normal} -O unfiltered.vcf.gz
    """
}

// Step 5b: filter, keep PASS
process FILTER {
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(ref), path(ref_index), path(dict)
    tuple path(vcf), path(tbi), path(stats)

    output:
    path "pass.vcf.gz"

    script:
    """
    gatk FilterMutectCalls -R ${ref} -V ${vcf} -O filtered.vcf.gz
    gatk SelectVariants -R ${ref} -V filtered.vcf.gz --exclude-filtered -O pass.vcf.gz
    """
}

// Step 6: compare with NIST's known mutations (the same script as Snakemake)
process COMPARE {
    publishDir params.outdir, mode: 'copy'

    input:
    path calls
    path known
    path trusted
    path compare_script

    output:
    path "comparison"

    script:
    """
    python ${compare_script} ${calls} ${known} ${trusted} ${params.region} comparison
    """
}

workflow {
    reference = PREPARE_REFERENCE(file("${params.data}/reference.fa"))

    reads = channel.of(params.tumour, params.normal)
        .map { s -> tuple(s, file("${params.data}/${s}_R1.fq.gz"), file("${params.data}/${s}_R2.fq.gz")) }

    bams = MARKDUP(ALIGN(reads, reference))
    tumour = bams.filter { b -> b[0] == params.tumour }
    normal = bams.filter { b -> b[0] == params.normal }

    passed = FILTER(reference, CALL(reference, tumour, normal))
    COMPARE(passed, file("${params.data}/known_mutations.vcf"),
            file("${params.data}/trusted_regions.bed"), file("scripts/compare_calls.py"))
}
