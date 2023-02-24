# mu_calculator
Repo for gene weight (mu) calculation 

## Installation
You can either use BigPipeline environment or follow these steps:
1. Clone repo `git clone https://github.com/LichtargeLab/mu_calculator.git`
2. Create conda environment `conda env create -f ./mu_calculator/environment.yml`
3. Activate conda environment `conda activate pyMu`
4. Install package `pip install -e mu_calculator`


## Usage
Required arguments:
| Argument                | Descripion |
| ---------------------- |--------------------- |
| --VCF                | Path to annotated VCF file |
| --samples            |Path to two-column comma separated file (.csv) with sample IDs in the first column and patient labels (cases=1, controls=0) in the second column. There should be no header row|
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
python -m MuDiff --VCF Path/to/vcf_file.vcf.gz --samples Path/to/samples_file.csv --savepath save/directory/ --cores 20 --maxaf 0.01
```

## Output
1. Output "mu.tsv" is a dataframe with gene names in the first column, mu of the cases in the second column and mu of the controls in the third column.
2. Output "distance_matrix.tsv" is a dataframe with gene names in the first column and distance (mu_control - mu_case) in the second column. The difference in mu for controls and cases for 1000 random shuffled labels are presented in column names "1" to "1000" and z-scores are presented in the last column.
