# mu_calculator
Repo for gene weight (mu) calculation 

## Installation
You can either use BigPipeline environment or follow these steps:
1. git clone https://github.com/LichtargeLab/mu_calculator.git
2. conda env create -f ./mu_calculator/environment.yml
3. conda activate pyMu


## Usage
Required arguments:
| Argument                | Descripion |
| ---------------------- |--------------------- |
| --VCF                | Path to annotated VCF file |
| --samples            |Path to one-column tab seperated (.txt, .tsv) file with sample IDs in the first column. There should be no header row|
| --maxaf  | sets maximum allele frequency threshold, so mu will be calculated for variants with frequency less than this threshold |
| --savepath           | Path for output file |
| --cores              | number of cpus to use |

Optional arguments:
| Argument                 | Descripion |
| ---------------------- |--------------------- |
| --Ann      | Variant annotation pipeline used (options: ANNOVAR, VEP / default: VEP) |
| --ref      | genome reference (options: hg19, hg38 / default: hg38) |
| --GeneLength      | Path to gene length file. ( default: ./refs/gene_length.csv) |



## Command line example
```bash
#set your working directory to mu_calculator
cd ./mu_calculator
#run mu.py
python mu.py --VCF Path/to/vcf_file.vcf.gz --samples Path/to/samples_file.tsv --savepath save/directory/ --cores 20 --maxaf 0.01
```

## Output
Output is a dataframe with gene names in the first column, gene length in the second column, sum EA in the third column, and mu in the fourth column. 


