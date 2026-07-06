
params.dir="/home/zl436/palmer_scratch/oyster_genome"
params.reads="${params.dir}/HiC/*_{1,2}.fastq.gz"

process get_trimming_parameters{

        container 'chesterli214/gsexp_amd64_v0.6'

        input:
        tuple val(replicateID),path(QC_results_1),path(QC_results_2)

        output:
        tuple val(replicateID),env(trimming_params)

        script:
        """
        unzip ${QC_results_1}
        unzip ${QC_results_2}
        qual1=\$(Get_trimming_parameter_qual ${QC_results_1.baseName}/fastqc_data.txt)
        qual2=\$(Get_trimming_parameter_qual ${QC_results_1.baseName}/fastqc_data.txt)
        trimming_params="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:10:30 MINLEN:50"
        """
}

process quality_trimming{

        container 'quay.io/biocontainers/trimmomatic:0.36--3'

        cpus 5

        input:
        tuple val(replicateID),path(reads)

        output:
        tuple val(replicateID),path("${replicateID}_*_trimmed.fastq.gz")

        script:
        """
        trimmomatic PE -threads ${task.cpus} ${reads} ${replicateID}_1_trimmed.fastq.gz ${replicateID}_1_unpaired.fastq.gz ${replicateID}_2_trimmed.fastq.gz ${replicateID}_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:15 TRAILING:5 SLIDINGWINDOW:10:30 MINLEN:50
        """
}

process remove_polyG{
        container 'quay.io/biocontainers/fastp:0.24.0--heae3180_1'

        cpus 4

        input:
        tuple val(replicateID),path(reads)

        output:
        tuple val(replicateID),path("${replicateID}_*_trimmed_grmv.fastq.gz")


        script:
        """
        fastp -i ${replicateID}_1_trimmed.fastq.gz -I ${replicateID}_2_trimmed.fastq.gz -o ${replicateID}_1_trimmed_grmv.fastq.gz -O ${replicateID}_2_trimmed_grmv.fastq.gz --trim_poly_g --trim_poly_x --poly_x_min_len 10 -l 50 --thread ${task.cpus}
        """
}

process reads_quality_control_cp1{

        container 'quay.io/biocontainers/fastqc:0.11.9--0'

        cpus 4

        input:
        tuple val(replicateID),path(reads)

        output:
        tuple val(replicateID),path("*_1_trimmed_grmv_fastqc.zip"),path("*_2_trimmed_grmv_fastqc.zip")

        script:
        """
        fastqc -t ${task.cpus} ${reads}
        """
}

workflow{
        reads_ch=Channel.fromFilePairs(params.reads,checkIfExists: true)
        reads_ch.view()

        quality_trimming(reads_ch)
        quality_trimming.out.view()

        remove_polyG(quality_trimming.out)
        trimmed_reads_ch=remove_polyG.out
        trimmed_reads_ch.view()

        trimmed_qc_ch=trimmed_reads_ch|reads_quality_control_cp1
        trimmed_qc_ch.view()

}





