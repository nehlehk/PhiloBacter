# PhiloBacter
PhiloBacter is a software package that detects recombination events and reconstructs a phylogenetic tree of bacterial genomes.

In practice, PhiloBacter consists of two tools, the first can be used independently as a recombination detection tool. To simplify the usage of these tools, validate their result, and compare them to the other state-of-the-art methods, a pipeline has been built using Nextflow software.

# Pipeline introduction

The application of this pipeline can be classified into two modes. The first mode is for all users who have sequences of a group of bacteria and are interested to learn about recombination events or want to have an authentic and reliable phylogenetic tree of those bacteria that were not impacted by recombination. These users can choose the analysis mode to use the pipeline. 


### Dependencies

Nextflow

conda


### Instllation
```
 # Download repo
git clone https://github.com/nehlehk/PhiloBacter.git
```

If Nextflow doesn't appear to create the conda environment properly. Create manually.
```
conda env create -f PhiloBacter.yml
conda activate PhiloBacter
```


### 1) Analysis mode
In this case, the pipeline consists of two steps. The first step builds an initial tree from the sequence using RAxML. In the second step, suppose the user is only looking for recombination events and their boundaries in the genome or, in other words, looks for the mosaic pattern of the genome. In that case, the pipeline can respond to this request with good accuracy. In addition, the other option is also available if the user is interested in the phylogeny tree. However, the user does not need to set any parameters; if one runs the pipeline using the default settings, all the steps will be done automatically. 

The bonus capability is that this pipeline is not only specific to PhiloBacter. This is possible, if the user also wants to use two other well-known tools in this field, such as [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML) and [Gubbins](https://github.com/nickjcroucher/gubbins), and collect the final trees of all three methods. Here is the command to use the pipeline for real datasets.

```
./nextflow main.nf --mode Analysis --seq genome.fasta  --method pb,cfml,gub
```
"--seq" provides an option to introduce the desired aligned sequence (fasta format) to the pipeline, and "-- method" is used to specify the method for the analysis. pb for PhiloBacter, cfml for ClonalFrameML and gub for Gubbins.







### 2) Simulation mode
The second mode of the pipeline is specified for experts and developers of recombination detection tools in bacterial genomes who want to compare different approaches. In this case, the pipeline consists of five main steps.
##### Step 1- Simulation: The first step is to simulate the clonal and local trees.
##### Step 2- Generating Sequences: To generate the alignment, we use the local and clonal trees created in the last step as input to [Seq-Gen](https://github.com/rambaut/Seq-Gen) software. The default evolution model is GTR.
##### Step 3- Constructing an initial tree: We have used [RAxML](https://github.com/stamatak/standard-RAxML) to build the initial tree based on the sequences that are the output of Seq-Gen.

##### Step 4- Recombination detection and tree inference: In this step, in addition to PhiloBacter, Gubbins and CFML tools can be run on the same simulated dataset. The output of this step is recombination events and tree inference.

##### Step 5- Analysis: The last step of the pipeline consists in evaluating recombination event and phylogenetic tree estimates. Specifically, inferred tree are compared to the clonal (true) tree.

There are various parameters for data simulation that can be adjusted according to the user's requirement, such as genome and recombination length or recombination rate, and time to the most recent ancestor (tMRCA). 

```
./nextflow main.nf --mode sim --genome 10 --genomelen 100000 --recomlen 500 --tMRCA 0.01 --recomrate 0.01 --nu_sim 0.05
