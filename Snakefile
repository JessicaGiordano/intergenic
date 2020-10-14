rule basic_annotation_gff3:
	output:
	 gff3="../task/gencode/{organism}/{version}/basic.annotation.gff3.gz"  
	shell:"wget -O {output.gff3} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{wildcards.organism}/release_{wildcards.version}/gencode.v{wildcards.version}.basic.annotation.gff3.gz"

rule basic_annotation_gtf:
	output:
	 gtf="../task/gencode/{organism}/{version}/basic.annotation.gtf.gz"  
	shell:"wget -O {output.gtf} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{wildcards.organism}/release_{wildcards.version}/gencode.v{wildcards.version}.basic.annotation.gtf.gz"

rule genePred:
	input:"../task/gencode/{organism}/{version}/basic.annotation.gtf" 
	output:"../task/gencode/{organism}/{version}/basic.annotation.genePred" 
	shell:"gtfToGenePred {input} {output}"

rule primary_assembly_annotation:
	output:
	 gz="../task/gencode/{organism}/{version}/primary_assembly.annotation.gff3.gz"  
	shell: "wget -O {output.gz} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{wildcards.organism}/release_{wildcards.version}/gencode.v{wildcards.version}.primary_assembly.annotation.gff3.gz && zcat {output.gz}"

rule zcat:
	input:"../task/gencode/{organism}/{version}/{filename}.annotation.gff3.gz"  
	output:"../task/gencode/{organism}/{version}/{filename}.annotation.gff3"
	shell: "zcat {input} > {output}"  

rule gtf2bed:
	input:"{filename}.annotation.gff3"
	output:"{filename}.annotation.bed"
	shell:"""gff2bed < {input} | sort -k1,1 -k2,2n  > {output}"""

rule gene_content:	
	input:"../task/gencode/{organism}/{version}/{filename}.annotation.bed"
	output:"../task/gencode/{organism}/{version}/{filename}.annotation.gene.bed"
	shell:"""awk 'BEGIN {{OFS="\\t"}} {{ if ($8 == "gene") {{ print}} }}' {input} > {output} """

rule flat_bed:
	input: "../FLAT/gencode.v33.basic.annotation.pred_for_picard"
	output: "../FLAT/gencode.v33.basic.annotation.pred_for_picard.bed"
	shell: """awk 'BEGIN {{OFS="\\t"}} {{print $3,$5,$6,$1}}' {input} > {output}"""

rule chrom_size:
	output: "../task/genome/{ref_genome}/{ref_genome}.chrom.sizes"
	shell: "wget -c https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.ref_genome}/bigZips/{wildcards.ref_genome}.chrom.sizes -O {output}"

rule chrom_coords:
	input: "../task/genome/{ref_genome}/{ref_genome}.chrom.sizes"
	output:"../task/genome/{ref_genome}/{ref_genome}.chrom.coords"
	shell: """awk 'BEGIN {{OFS="\\t" }} {{print $1,1, $2}}' {input} > {output}"""

rule intergenic_region:
	input:
	 whole_genome="../task/genome/{ref_genome}/{ref_genome}.chrom.coords",
	 gene_coords="../task/gencode/{organism}/{version}/basic.annotation.genePred"
	output:
	 intergenic_bed="../task/gencode/{organism}/{version}/{ref_genome}/intergenic.bed"
	shell:
	 "bedtools subtract -a {input.whole_genome} -b {input.gene_coords} > {output.intergenic_bed}"

rule integenic_content:
	input:
         bam="../{sample}/{sample}.Aligned.sortedByCoord.out.bam",
         intergenic_bed="../task/gencode/{organism}/{version}/{ref_genome}/intergenic.bed"
	output:"../results/{sample}/{organism}/{ref_genome}/{version}/{sample}.Aligned.sortedByCoord.out.intergenic.bed"
	shell:"bedtools intersect -a {input.bam} -b {input.intergenic_bed} -wo -bed > {output}"

rule full_intergenic_reads:
	input:"../results/{sample}/{organism}/{ref_genome}/{version}/{sample}.Aligned.sortedByCoord.out.intergenic.bed"
	output:"../results/{sample}/{organism}/{ref_genome}/{version}/{sample}.Aligned.sortedByCoord.out.intergenic.full.bed"
	shell:"""awk 'BEGIN {{OFS="\\t"}} {{ if ($16 == $3-$2) {{ print}} }}' {input} > {output} """

