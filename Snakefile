import os
import glob
import multiprocessing
from python import down_ascn 
from python import fasta_cleaner
from python import concatenate_fastas
from python import newick_tree_visualizer


# Dataset
sample_names = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob("ascn/*.txt")]

rule all:
    input:
        "Tree.svg"

rule Download_ASCN:
    input:
        ascn="ascn/{sample}.txt"
    output:
        fasta="fasta/{sample}.fasta"  
    run:
        down_ascn.download_sequences_from_asn_file(input.ascn,output.fasta)

rule Clean:
    input:
        fasta="fasta/{sample}.fasta"  
    output:
        clean="fasta_clean/{sample}.fasta"  
    run:
        fasta_cleaner.clean_fasta(input.fasta,output.clean)

rule Alignment:
    input:
        clean="fasta_clean/{sample}.fasta"  
    output:
        align="fasta_align/{sample}.fasta"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        mafft --auto --thread {params.threads} {input.clean} > {output.align}
        """

rule Trim:
    input:
        align="fasta_align/{sample}.fasta"
    output:
        trim="trim/{sample}.fasta"
    params:
        threshold=80    
    shell:
        """
        python3 python/trim_fasta_edges.py {input.align} {output.trim} {params.threshold}
        """

rule Best_AICc_model:
    input:
        trim="trim/{sample}.fasta"
    output:
        "trim/{sample}.txt",
        "trim/{sample}.fasta.ckp",
        "trim/{sample}.fasta.log",
        "trim/{sample}.fasta.out",
        "trim/{sample}.fasta.topos",
        "trim/{sample}.fasta.tree"
    conda:
        "envs/yamlfile.yaml"
    shell:
        "modeltest-ng -i {input.trim} -t ml > {output}"  

rule obtain_best_AICc_model:
    input:
        data="trim/{sample}.txt"
    output:
        "trim/{sample}Model.txt"
    shell:
        "grep 'Best model according to AICc' -A 2 {input.data} | tail -n 1 | sed 's/^.* //' > {output}"

rule add_seq_length_to_model:
    input:
        trim=expand("trim/{sample}Model.txt", sample=sample_names)
    output:
        prim="tree/prim.part"
    shell:
        "python3 python/sequence_model_processor.py ./trim {output.prim}"

rule Concatenated:
    input:
        prim="tree/prim.part"
    output:
        align="tree/concatenated.fasta"
    run:
        concatenate_fastas.concatenate_fastas("./trim",output.align)

rule maximum_likelihood_tree_step_1:
    input:
        model="tree/prim.part",
        msa="tree/concatenated.fasta"
    output:
        "tree/concatenated.fasta.raxml.bestTree",
        "tree/concatenated.fasta.raxml.bestModel",
        "tree/concatenated.fasta.raxml.mlTrees",
        "tree/concatenated.fasta.raxml.rba",
        "tree/concatenated.fasta.raxml.startTree"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count(),
        outgroup="Macrobiotus_rybaki" 
    shell:
        """
        raxml-ng --msa {input.msa} --model {input.model} --threads {params.threads} --seed 333 --tree pars{{100}},rand{{100}} --force perf_threads --outgroup {params.outgroup}
        """
rule maximum_likelihood_tree_step_2:
    input:
        model="tree/prim.part",
        msa="tree/concatenated.fasta"
    output:
        "tree/concatenated.fasta.raxml.bootstraps",
        "tree/concatenated.fasta.raxml.log",
        "tree/concatenated.fasta.raxml.rba"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count(),
        bootstrap_trees=10000000,
        outgroup="Macrobiotus_rybaki" 
        #outgroup="Macrobiotus_rybaki,Sisubiotus_spectabilis,Mesobiotus_datanlanicus" 
    shell: 
        """
        raxml-ng --bootstrap --msa {input.msa} --model {input.model} --threads {params.threads} --seed 1 --bs-trees autoMRE{{{params.bootstrap_trees}}} --force perf_threads --outgroup {params.outgroup}
        """
rule maximum_likelihood_tree_step_3:
    input:
        tree="tree/concatenated.fasta.raxml.bestTree",
        bootstraps="tree/concatenated.fasta.raxml.bootstraps"
    output:
        "tree/concatenated.fasta.raxml.bestTree.raxml.support",
        "tree/concatenated.fasta.raxml.bestTree.raxml.log"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count(),
    shell:
        "raxml-ng --support --tree {input.tree} --bs-trees {input.bootstraps} --threads {params.threads} --force perf_threads"  

rule build_tree_mafft:
    input:
        tree="tree/concatenated.fasta.raxml.bestTree.raxml.support"
    output:
        img="Tree.svg"
    run:
        newick_tree_visualizer.tree_generator(input.tree, output.img)        