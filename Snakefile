import Bio

# Dataset
sample_names = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob("dataset/*.fasta")]

# Define rule to generate all desired outputs
rule all:
    input:
        expand("results/mafft/{sample}/tree/{sample}Tree.svg", sample=sample_names),    
        expand("results/clustal_omega/{sample}/tree/{sample}Tree.svg", sample=sample_names),
        expand("results/Muscle/{sample}/tree/{sample}Tree.svg", sample=sample_names),
        expand("iqtree/{sample}/mafft/mafft.log", sample=sample_names),
        expand("iqtree/{sample}/clustal_omega/clustal_omega.log", sample=sample_names),
        expand("iqtree/{sample}/Muscle/Muscle.log", sample=sample_names)

#Alignment rules
rule ASCN_Download:
    input:
        ASCN="ASCN/{sample}.txt"
    output:
        fasta="results/mafft/{sample}/{sample}fileAligned.fasta"
    params:
        library=nuccore
    shell:
        """
        python3 ASCN_Download.py nuccore ASCN > {sample}.fasta
        """
