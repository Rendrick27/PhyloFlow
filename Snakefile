import Bio
import glob
import multiprocessing
# Dataset
sample_names = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob("ASCN/*.txt")]

# Define rule to generate all desired outputs
rule all:
    input:
        expand("fasta/{sample}_align.fasta", sample=sample_names)   

#Alignment rules
rule ASCN_Download:
    input:
        ASCN="ASCN/{sample}.txt"
    output:
        fasta="fasta/{sample}.fasta"
    params:
        library="nuccore"
    shell:
        """
        python3 ASCN_Download.py {params.library} {input.ASCN} > {output.fasta}
        """
rule Alignment:
    input:
        fasta="fasta/{sample}.fasta"
    output:
        align="fasta/{sample}_align.fasta"
    params:
        threads=multiprocessing.cpu_count()    
    shell:
        """
        mafft --auto --thread -{params.threads} {input.fasta} > {output.align}
        """