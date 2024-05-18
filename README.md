# Phylo Flow 

## Description
This tool, developed within the scope of ASB.
![Pipeline](./extras/pictures/pipeline.png)


## Requirements
* at least one txt files with ascn number;
* <a href= "https://www.python.org/"> Python </a> & <a href= "hhttps://biopython.org/"> BioPython </a> 
* <a href= "https://mafft.cbrc.jp/alignment/server/index.html"> Mafft </a>;
* <a href= "https://github.com/ddarriba/modeltest"> Modeltest-ng </a>;
* <a href= "https://github.com/amkozlov/raxml-ng"> Raxml-ng </a>;
* <a href= "https://toytree.readthedocs.io/en/latest/index.html"> Toytree </a>;
* <a href= "https://snakemake.readthedocs.io/en/stable/#"> Snakemake </a>;
  

### txt file formart
```bash
#TxT formart
Sequence_name;ascn_number

#Example
Paramacrobiotus_gadabouti_sp._nov._MD50.1;OP394210
```
## Installation
```bash
# Download the project
wget https://github.com/Rendrick27/PhyloFlow/archive/refs/heads/main.zip

# Unzip the folder
unzip main.zip
```
Then copy your .txt files into ascn folder.

## Usage
```bash
# Navigate to the Snakemake directory
cd main.zip
```
Then, run the following command:
```bash
snakemake --use-conda all --cores 1
```
After that it will show a Tree.svg

## Docker
### Build it
```bash
docker build -t {image_name} .

docker run -it --name {container_name} {image_name} /bin/bash
```
### Pull docker image
```bash
docker pull rendrick27/phylo_flow:latest
```

## Settings
You may adjust settings in the Snakemake file, such as threads and bootstraps in params, but remember that using more threads may cause more issues.

## Credits
<p> <a href= "https://github.com/Rendrick27"> Rendrick Carreira - 201901365 </a> </p>

## License
GPLv3
