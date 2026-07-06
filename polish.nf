params.dir="/home/zl436/palmer_scratch/oyster_genome"
params.reads_hifi = "${params.dir}/PacBio/hifi_reads/11S_pb.fastq.gz"
params.contigs="${params.dir}/Contigs/*.p_ctg.gfa"
params.reads_illumina = "${params.dir}/WGS/*{1,2}_trimmed_grmv.fastq.gz"
params.reads_hic = "${params.dir}/HiC/11S_HiC_{1,2}_trimmed_grmv.fastq.gz"


process gfa_to_fa{
	container 'quay.io/biocontainers/gfatools:0.5--h577a1d6_5'

	cpus 4

	input:
	tuple val(assembly_name),path(gfa_file)

	output:
	tuple val(assembly_name),path("${assembly_name}.fa")	

	script:
	"""
	gfatools gfa2fa ${gfa_file} > ${assembly_name}.fa
	"""
}




process evaluate_contigs{
	container 'quay.io/biocontainers/quast:5.3.0--py311pl5321hc84137b_1'

	cpus 8

	input:
	tuple val(assembly_name),path(fa_files)

	output:
	tuple val(assembly_name),path("${assembly_name}_quast_report")

	script:
	"""
	quast -t ${task.cpus} -o ${assembly_name}_quast_report ${fa_files}
	"""

}


process evaluate_genes{
	container 'quay.io/biocontainers/busco:5.8.2--pyhdfd78af_0'

	cpus 8

	input:
	tuple val(assembly_name),path(fa_files)

	output:
	tuple val(assembly_name),path("${assembly_name}_busco_report")

	script:
	"""
	busco -c ${task.cpus} -m genome -i ${fa_files} -o ${assembly_name}_busco_report -l mollusca_odb10
	"""

}

process get_kmer_content{
	container 'quay.io/biocontainers/meryl:1.4.1--h9948957_2'

	cpus 8

	input:
	tuple val(sample_id),path(reads1),path(reads2)

	output:
	path("${sample_id}_kmer_db.meryl")
	
	script:
	"""
	meryl k=21 count output ${sample_id}_kmer_db.meryl ${reads1} ${reads2}
	"""

}

process evaluate_kmers{
	
	container 'danylmb/merqury:1.4.1-build2'

	cpus 8

	input:
	tuple val(assembly_name),path(fa_files),path(kmer_db)

	output:
	tuple val(assembly_name),path("${assembly_name}_merqury_output*")

	script:
	"""
	merqury.sh ${kmer_db} ${fa_files} ${assembly_name}_merqury_output
	"""
}

process prepare_genome_bwa{

        container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

        input:
        path genome

        output:
        path "${genome.baseName}.*"

        script:
        """
        bwa index ${genome} -p ${genome.baseName}
        """
}

process prepare_genome_samtools{

        container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

        input:
        path genome

        output:
        path "${genome}.fai"

        script:
        """
        samtools faidx ${genome}
        """
}

process reads_mapping{
        container 'quay.io/biocontainers/bwa:0.7.18--he4a0461_1'

        cpus 4

        input:
        path genome_index
        tuple val(replicateID),path(trimmed_reads)

        output:
        tuple val(replicateID),path("*.sam")

        script:
        """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${replicateID}\\tSM:${replicateID}.SM\\tLB:${replicateID}.LB' -o ${replicateID}.sam ${genome_index[0].baseName} ${trimmed_reads}
        rm ${trimmed_reads}
        """
}

process post_alignment_processing{
        container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

        cpus 4

        input:
        tuple val(replicateID),path(samfile)

        output:
        tuple val(replicateID),path("${replicateID}_st_dr_filtered_gs.bam"),path("${replicateID}_st_dr_filtered_gs.bam.bai")

        script:
        """
        samtools fixmate -@ ${task.cpus} -O BAM ${replicateID}.sam ${replicateID}.bam
        rm ${replicateID}.sam
        samtools sort -@ ${task.cpus} -O BAM -o ${replicateID}_sorted.bam ${replicateID}.bam
        rm ${replicateID}.bam
        samtools rmdup ${replicateID}_sorted.bam ${replicateID}_st_dr.bam
        rm ${replicateID}_sorted.bam
        samtools view -@ ${task.cpus} -O BAM -q 30 -o ${replicateID}_st_dr_filtered_gs.bam ${replicateID}_st_dr.bam
	rm ${replicateID}_st_dr.bam
        samtools index ${replicateID}_st_dr_filtered_gs.bam
        """
}

process polish_assembly{
    
	container 'pvstodghill/nextpolish:latest'

	cpus 16
	memory '32 GB'

	input:
	path(contigs)
	path(hifi_reads)
	tuple val(sample_id),path(read1),path(read2)

	output:
	path("*01_rundir/genome.nextpolish.fasta")

	script:
	"""
	ls ${read1} ${read2}>sgs.fofn
	ls ${hifi_reads} ${contigs}>hifi.fofn


	cat << 'EOF' > run.cfg
	[General]
	job_type = local
	job_prefix = nextPolish
	task = best
	rewrite = yes
	rerun = 3
	parallel_jobs = 6
	multithread_jobs = 5
	genome = ./${contigs}
	genome_size = auto
	workdir = ./01_rundir
	polish_options = -p {multithread_jobs}

	[sgs_option]
	sgs_fofn = ./sgs.fofn
	sgs_options = -max_depth 100 -bwa

	[hifi_option]
	hifi_fofn = ./hifi.fofn
	hifi_options = -min_read_len 1k -max_depth 100
	hifi_minimap2_options = -x map-pb
	EOF

	nextPolish run.cfg
	"""
}

workflow{

	Channel.from([["hap1", file("${params.dir}/Contigs/contigs.bp.hap1.p_ctg.gfa")],["hap2", file("${params.dir}/Contigs/contigs.bp.hap2.p_ctg.gfa")],["primary", file("${params.dir}/Contigs/contigs.bp.p_ctg.gfa")]]).set{gfa_files}
	Channel.fromFilePairs(params.reads_illumina,flat:true).set{illumina_ch}
	
	gfa_files.view()
	illumina_ch.view()

	fa_ch=gfa_files|gfa_to_fa
	fa_ch.view()

	fa_ch|evaluate_contigs
	evaluate_contigs.out.view()

	fa_ch|evaluate_genes
	evaluate_genes.out.view()

	kmer_ch=illumina_ch|get_kmer_content
	kmer_ch.view()	
	
	comb_kmer_ch=fa_ch.combine(kmer_ch)
	comb_kmer_ch.view()	

	evaluate_kmers(comb_kmer_ch)
	evaluate_kmers.out.view()

	contig_ch=Channel.from(file("${params.dir}/Contigs/primary.fa"))	
	contig_ch.view()	
	hifi_ch=Channel.fromPath("${params.reads_hifi}")
	hifi_ch.view()

	polished_contig_ch=polish_assembly(contig_ch,hifi_ch,illumina_ch)

}



