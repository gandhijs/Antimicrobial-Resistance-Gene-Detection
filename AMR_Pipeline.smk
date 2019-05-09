SAMPLE, NUM = glob_wildcards("/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample, [^_]+}_{num}.fastq")

rule all:
    input:
        # get the raw paired end reads
        expand("/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_{num}.fastq", zip, sample=SAMPLE, num=[1, 2]),
        # mash concatenated fastq file
        expand("/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}.fastq", sample=SAMPLE),
        # mash output
        expand("/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}.fastq.msh", sample=SAMPLE),
        # top mash hit
        expand("/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}_mash_distances.tab", sample=SAMPLE),
        # mash information
        expand("/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}_mash_info.txt", sample=SAMPLE),  
        # readMetrics from CG_Pipeline results for PE read 1
        expand("/home/gandhijs/isolates_flagged_mcr-1/CG_pipeline_results/{sample}_1_read_metrics.txt", sample=SAMPLE),
        # readMetrics from CG_Pipeline results for PE read 2 
        expand("/home/gandhijs/isolates_flagged_mcr-1/CG_pipeline_results/{sample}_2_read_metrics.txt", sample=SAMPLE),
	# trimmomatic with PE read 1 and PE read 2 to generate 4 output files
        expand("/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_1P.fastq", sample=SAMPLE),
	expand("/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_1U.fastq", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_2P.fastq", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_2U.fastq", sample=SAMPLE),
        # shovill results 
        expand("/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/contigs.fa", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/shovill.corrections", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/spades.fasta", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/contigs.gfa", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/shovill.log", sample=SAMPLE),
        # quast results
        expand("/home/gandhijs/isolates_flagged_mcr-1/quast/{sample}_quast_output/report.txt", sample=SAMPLE),
        # resfinder results
        #expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/Hit_in_genome_seq.fsa", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/results_table.txt", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/results.txt", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/Resistance_gene_seq.fsa", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/results_tab.txt", sample=SAMPLE),
        expand("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/tmp", sample=SAMPLE),


"""Mash will take as input two paired-end raw fastq files and refseq database. It will first concatenate the paired end reads. The next step will be to sketch the concatenated reads. This is followed by performing the mash info to verify the contents and finally estimate pairwise distances. The final step of this rule will be to sort the top hits by the p-value and write the top hit to a csv file. 
"""

rule mash:
    input:
        read1="/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_1.fastq",
        read2="/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_2.fastq",
        refseq="/home/gandhijs/bin/mash-Linux64-v2.1/refseq.genomes.k21s1000.msh"
    output:
        combined="/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}.fastq",
        msh="/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}.fastq.msh",
        tab="/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}_mash_distances.tab",
        mash_info="/home/gandhijs/isolates_flagged_mcr-1/mash_output/{sample}_mash_info.txt"
    shell:
        """
        set +o pipefail
        cat {input.read1} {input.read2} > {output.combined}
        mash sketch -m 2 {output.combined}
        mash info {output.msh} > {output.mash_info}
        mash dist {input.refseq} {output.msh} > {output.tab}
        sh consolidate_mash.sh
        """
"""run_assembly_readMetrics.pl will take as input two paired-end raw fastq files and output read data to text file. Finally the read metrics will be consolidated to a csv file.
"""

rule readMetrics:
    input:
        read1="/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_1.fastq",
	read2="/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_2.fastq"
    output:
        read1="/home/gandhijs/isolates_flagged_mcr-1/CG_pipeline_results/{sample}_1_read_metrics.txt",
	read2="/home/gandhijs/isolates_flagged_mcr-1/CG_pipeline_results/{sample}_2_read_metrics.txt"
    threads:
        1
    shell:
        """
        run_assembly_readMetrics.pl --fast {input.read1} -e \"5000000" > {output.read1}
        run_assembly_readMetrics.pl --fast {input.read2} -e \"5000000" > {output.read2}
        sh consolidate_read_metrics.sh
        """

"""Trimmomatic will take as input two paired-end raw fastq files and output four output files, two for the 'paired' output in which both of the reads survived the processing, and two for the 'unpaired' output in which a read survived but the partner read did not. 
"""
rule trimmomatic:
    input:
        read1="/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_1.fastq",
        read2="/home/gandhijs/isolates_flagged_mcr-1/raw_isolate_data/{sample}_2.fastq"
    output:
        read1_P="/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_1P.fastq",
        read1_U="/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_1U.fastq",
        read2_P="/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_2P.fastq",
        read2_U="/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_2U.fastq"
    threads:
        1    
    shell:
        "sh trimmomatic PE {input.read1} {input.read2} {output.read1_P} {output.read1_U} {output.read2_P}"
        " {output.read2_U} SLIDINGWINDOW:4:30"

"""Shovill will take as input two trimmed 'paired' end read files to conduct genome assembly. 
"""
rule shovill:
    input:
        read1="/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_1P.fastq",
        read2="/home/gandhijs/isolates_flagged_mcr-1/trimmomatic_results/{sample}_2P.fastq"    
    threads:
        48
    output:
        #directory("/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output")
        "/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/contigs.fa",
        "/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/shovill.corrections",
        "/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/spades.fasta",
        "/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/contigs.gfa",
        "/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/shovill.log"
    shell:
        """
        set +o pipefail
        shovill --cpu {threads} --outdir /home/gandhijs/isolates_flagged_mcr-1/assembly_output/{wildcards.sample}_shovill_output --R1 {input.read1} --R2 {input.read2}
        """

"""Quast will evaluate the quality of the genome assembly.
"""
rule quast:
    input:
        contigs="/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/contigs.fa"
    output:
        report = "/home/gandhijs/isolates_flagged_mcr-1/quast/{sample}_quast_output/report.txt",
    threads:
        1
    shell:
        """
        quast.py {input.contigs} --output-dir /home/gandhijs/isolates_flagged_mcr-1/quast/{wildcards.sample}_quast_output
        """

"""ResFinder will identify acquired antimicrobial resistant genes from the assembled contigs.
"""
rule resfinder:
    input:
        contig="/home/gandhijs/isolates_flagged_mcr-1/assembly_output/{sample}_shovill_output/contigs.fa"
    output:
        "/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/Hit_in_genome_seq.fsa",
        "/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/results_table.txt",
        "/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/results.txt",
        "/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/Resistance_gene_seq.fsa",
        "/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/results_tab.txt",
        directory("/home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{sample}_resfinder_output/tmp")
    shell:
        """
        set +o pipefail
        python3 ~/bin/resfinder/resfinder.py -i {input.contig} -o /home/gandhijs/isolates_flagged_mcr-1/resfinder_output/{wildcards.sample}_resfinder_output -p ~/bin/resfinder/resfinder_db -b blastn -t 0.90 -l 0.60
        """
