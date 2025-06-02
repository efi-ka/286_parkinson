SAMPLES, = glob_wildcards("parksamples/{sample}.vcf") # Global wildcards for matching the samples' names based on the wildcard pattern


rule all:
	input:
		expand("parksamples/{sample}_mdec_vep_GQ_rare_pdgenes_spliceai_CADD.xlsx", sample = SAMPLES) #specifies the final output targets. Creating a list of filenames



rule decompose_multiallelic_vt:                 #rule for decomposing multiallelic regions
    input:
        "parksamples/{sample}.vcf"
    output:
        "parksamples/{sample}_mdec.vcf"
    shell:
        "vt decompose -s {input} -o {output}"

rule vep_annotate_file:                        #rule for annotating the vcf file with vep running in singularity
	input:
		"parksamples/{sample}_mdec.vcf",
	output:
		"parksamples/{sample}_mdec_vep.vcf"
	shell:
		"export SINGULARITY_BIND+=/home/efthymia/singularity_VEP ; singularity exec /home/efthymia/singularity_VEP/image/vep.sif vep\
		--dir_cache /home/efthymia/singularity_VEP/cache\
		--vcf\
		--format vcf\
		--force_overwrite\
		--offline\
		--hgvs\
		--everything\
		--pick\
		--fork 4\
		--failed 1\
		--input_file {input}\
		--output_file {output}\
		--fasta /home/efthymia/singularity_VEP/cache/homo_sapiens/110_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz\
	  	--plugin dbNSFP,/home/efthymia/singularity_VEP/plugins/dbNSFP4.4a_grch38.gz,/home/efthymia/singularity_VEP/dbNSFP_replacement_logic,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_score,MutationAssessor_pred,\
		FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,Reliability_index,M-CAP_score,M-CAP_pred,REVEL_score,MutPred_score,CADD_raw,CADD_phred,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,DANN_score,VEST3_score,GERP++_NR,GERP++_RS,clinvar_clnsig,clinvar_trait\
		--plugin CADD,snv=/home/efthymia/singularity_VEP/plugins/whole_genome_SNVs.tsv.gz,indels=/home/efthymia/singularity_VEP/plugins/gnomad.genomes.r4.0.indel.tsv.gz\
		--fields 'Allele,Consequence,IMPACT,Existing_variation,Gene,SYMBOL,Codons,Amino_acids,Feature_type,Feature,HGVSc,cDNA_position,CDS_position,EXON,INTRON,HGVSp,Protein_position,BIOTYPE,CANONICAL,STRAND,DISTANCE,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,ENSP,SWISSPROT,TREMBL,UNIPARC,GENE_PHENO,DOMAINS,CCDS,CLIN_SIG,PHENO,PUBMED,SIFT,PolyPhen,AF,EUR_AF,SAS_AF,EAS_AF,gnomADe_AF,gnomADe_NFE_AF,gnomADe_FIN_AF,gnomADg_AF,gnomADg_NFE_AF,gnomADg_FIN_AF,\
		SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,M-CAP_score,M-CAP_pred,MutPred_score,REVEL_score,CADD_PHRED,DANN_score,GERP++_RS,clinvar_clnsig,clinvar_trait,GERP++_NR'"


rule filt_GQ_greaterthan20:                    #rule for selecting variants with genotype quality more than 20
	input:
		"parksamples/{sample}_mdec_vep.vcf"
	output:
		"parksamples/{sample}_mdec_vep_GQ.vcf"
	shell:
		"bcftools filter -i 'FMT/GQ>19' {input} -o {output}"

rule filt_rare_variants:                     #rule for selecting variants with MAF <= 0.02 in gnomAD non-Finnish Europeans population (genomes and exomes)
	input:
		"parksamples/{sample}_mdec_vep_GQ.vcf"
	output:
		"parksamples/{sample}_mdec_vep_GQ_rare.vcf"
	shell:
		"filter_vep -i {input} -o {output} -filter '(gnomADg_NFE_AF <= 0.02 or not gnomADg_NFE_AF) and (gnomADe_NFE_AF <= 0.02 or not gnomADe_NFE_AF)'"

rule filt_parkgenes:                        #rule for selecting variants located in genes that have been associated with PD
	input:
		"parksamples/{sample}_mdec_vep_GQ_rare.vcf"
	output:
		"parksamples/{sample}_mdec_vep_GQ_rare_pdgenes.vcf"
	shell:
		"python3 /home/efthymia/snakemake_park/scriptforgenesfun.py {input} > {output}"


rule filt_spliceai:                         #rule for annotating the vcf file with spliceai score
	input:
		file= "parksamples/{sample}_mdec_vep_GQ_rare_pdgenes.vcf",
		fasta="/home/efthymia/Downloads/hg38.fa"
	output:
		"parksamples/{sample}_mdec_vep_GQ_rare_pdgenes_spliceai.vcf"
	shell:
		"spliceai -I {input.file} -O {output} -R {input.fasta} -A grch38"

rule filt_CADD:                            #rule for filtering based on CADD_phred score
	input:
		"parksamples/{sample}_mdec_vep_GQ_rare_pdgenes_spliceai.vcf"
	output:
		"parksamples/{sample}_mdec_vep_GQ_rare_pdgenes_spliceai_CADD.vcf"
	shell:
		"filter_vep -i {input} -o {output} -filter 'CADD_PHRED >= 20'"


rule makeexcel:                            #rule for writting the results to an excel file
	input:
		"parksamples/{sample}_mdec_vep_GQ_rare_pdgenes_spliceai_CADD.vcf"
	output:
		"parksamples/{sample}_mdec_vep_GQ_rare_pdgenes_spliceai_CADD.xlsx"
	run:
		with open(input[0], 'r') as f:      #Check if there are any variant lines
			has_variants = any(not line.startswith('#') and line.rstrip() != '' for line in f)

		if not has_variants:
			print("no variants")
			open(output[0], 'w').close()
		else:
			shell("python3 /home/efthymia/snakemake_park/scriptforexcelnew.py {input} {output}")
