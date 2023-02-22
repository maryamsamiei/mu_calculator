import argparse
from joblib import Parallel, delayed
import numpy as np 
import os
import pandas as pd 
from pathlib import Path
import random
from tqdm import tqdm
from .vcf import *


def path(relative_path: str) -> Path:
    """
    Return path relative to module parent
    """
    return Path(__file__).parent / relative_path

def parse_args():
    """
    Parses the arguments.
    """
    parser = argparse.ArgumentParser(description="Mu arguments")
    parser.add_argument("--VCF", nargs="?", default="./",
        help="Location of Cohort VCF")
    parser.add_argument("--samples", nargs="?", default=path("samples.txt"),
        help="samples file path")
    parser.add_argument("--GeneLength", nargs="?", 
        default=path("refs/gene_length.csv"), 
        help="gene length file path")
    parser.add_argument("--ref", nargs="?", default="hg38", 
        choices=("hg19", "hg38"), help="genome reference file")
    parser.add_argument("--noX", default=False, action="store_true",
        help="Do not run analysis on chromosome X")
    parser.add_argument("--savepath", nargs="?", default="./",
        help="save path for output")
    parser.add_argument("--maxaf", type=float,
        help="maximum allele frequency cutoff")
    parser.add_argument("--Ann",default="VEP", choices=("VEP", "ANNOVAR"),
        help="EA annotation method")
    parser.add_argument("--degenerate",default=False, action="store_true",
        help="Use number of variant sites (degenerate distribution)")
    parser.add_argument("--cores", type=int, default=1,
        help="number of CPUs to use for multiprocessing")
    return parser.parse_args()

def _SumEA_degenerate(
        design_matrix: pd.DataFrame, 
        samples: list
        ) -> pd.Series:
    """
    Return degenerate SumEA based on variants in samples for every gene
    :design_matrix: variant matrix
        index: variant column (variant id from pysam.VariantRecord)
        columns: [ gene, variant, EA, *samples ]
    """
    # Only keep variants that appear at least once in samples
    variants = design_matrix[samples].any(axis=1, bool_only=True)
    sample_variants = variants[variants == True]
    dmatrix_sample = design_matrix[["gene", "EA"]]\
        .reindex(sample_variants.index)
    return dmatrix_sample.groupby("gene").EA.sum()

def compute_mu_diff(
        cases: list, 
        controls: list,
        gene_length: pd.DataFrame,
        design_matrix: pd.DataFrame,
        degenerate: bool
        ) -> tuple:
    """
    Compute Mu-diff for a given set of cases and controls
    :cases: list of case ids
    :controls: list of control ids
    :gene_length: Pandas DataFrame containing gene length
        :index: gene name (only those used for analysis)
        :gene_length: int
    :design_matrix: Pandas DataFrame containing SumEA per gene per sample
        - rows are samples
        - columns are genes

    :return: tuple containing pd.Series for mu-cases and mu-controls
    """
    # Subset sumEA matrix to genes and samples of study (cases and controls)
    genes = gene_length.index.to_list()
    if degenerate:
        SumEA_genes_case = _SumEA_degenerate(design_matrix, cases)
        SumEA_genes_case = SumEA_genes_case.reindex(genes)
        SumEA_genes_control = _SumEA_degenerate(design_matrix, controls)\
            .reindex(genes)
    else:
        design_matrix_case = design_matrix.loc[cases, genes]
        design_matrix_control = design_matrix.loc[controls, genes]
        # Compute SumEA associated to each gene for cases and controls
        SumEA_genes_case = np.sum(design_matrix_case, axis=0)
        SumEA_genes_control = np.sum(design_matrix_control, axis=0)

    # Compute Mu in cases
    expected_energy_case = np.sum(SumEA_genes_case) / \
                           np.sum(gene_length.gene_length)
    observed_energy_case = SumEA_genes_case / gene_length.gene_length
    mu_case = expected_energy_case / observed_energy_case
    
    # Compute Mu in controls
    expected_energy_control = np.sum(SumEA_genes_control) / \
                              np.sum(gene_length.gene_length)
    observed_energy_control = SumEA_genes_control / gene_length.gene_length
    mu_control = expected_energy_control / observed_energy_control

    return mu_case, mu_control

def main(args):
    # Create output directory if non existant
    os.makedirs(args.savepath, exist_ok=True)

    if args.Ann=="ANNOVAR":
        if args.ref=="hg19":
            ref = pd.read_csv(path("refs/refGene-lite_hg19.May2013.txt"), 
                              delimiter="\t", header=0, index_col="gene")
        elif args.ref=="hg38":
            ref = pd.read_csv(path("refs/refGene-lite_hg38.June2017.txt"), 
                              delimiter="\t", header=0, index_col="gene")
    elif args.Ann=="VEP":
        if args.ref=="hg19":
            ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh37.v75.txt"), 
                              delimiter="\t", header=0, index_col="gene")
        elif args.ref=="hg38":
            if args.noX:
                ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh38.v94.noX.txt"), 
                                  delimiter="\t", header=0, index_col="gene")
            else:
                ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh38.v94.txt"), 
                                  delimiter="\t", header=0, index_col="gene")

    samples = pd.read_csv(args.samples, header=None,index_col=0)
    controls = samples[samples.iloc[:,0]==0].index.astype(str).tolist()
    cases = samples[samples.iloc[:,0]==1].index.astype(str).tolist()
    total_samples = samples.index.astype(str).tolist()

    # Build SumEA matrix (sample in rows, genes in columns)
    if args.Ann=="ANNOVAR":
        matrix = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)\
            (args.VCF, gene, ref.loc[gene], total_samples, min_af=0,
             max_af=args.maxaf, af_field="AF", EA_parser="canonical")\
                for gene in tqdm(ref.index.unique()))

    if args.Ann=="VEP":
        if args.degenerate:
            matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP_degenerate)\
                (args.VCF, gene, ref.loc[gene], total_samples, 
                 min_af=0, max_af=args.maxaf) \
                    for gene in tqdm(ref.index.unique()))
        else:
            matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP)\
                (args.VCF, gene, ref.loc[gene], total_samples, 
                  min_af=0, max_af=args.maxaf) \
                    for gene in tqdm(ref.index.unique()))
             
    if args.degenerate:
        # row = variants; col = samples
        design_matrix = pd.concat(matrix, axis=0)
        matrix_genes = design_matrix.gene.unique().tolist()
    else:
        # row = samples; col = genes
        design_matrix = pd.concat(matrix, axis=1)
        matrix_genes = design_matrix.columns.tolist()

    ## reading gene length file
    gene_length = pd.read_csv(args.GeneLength, index_col=0)
    genes = set(matrix_genes)\
        .intersection(set(gene_length.index.tolist()))
    gene_length = gene_length.loc[genes]
    
    mu_case, mu_control = compute_mu_diff(cases, controls, gene_length, 
                                          design_matrix, args.degenerate)

    mu_matrix = pd.DataFrame(np.zeros((len(genes), 2)), index=genes, 
                             columns=["mu_case", "mu_control"])
    mu_matrix["mu_case"] = mu_case.copy()
    mu_matrix["mu_control"] = mu_control.copy()
    mu_matrix.to_csv(os.path.join(args.savepath, "mu.tsv"), sep="\t", 
                     header=True, index=True) 
    
    distance_matrix = pd.DataFrame(np.zeros((len(genes), 1)), index=genes, 
                                   columns=["distance"])
    distance_matrix["distance"] = mu_control - mu_case
    
    print("Performing randomization")
    for i in tqdm(range(1000)):
        cases1 = random.sample(total_samples, len(cases))
        controls1 = list(set(total_samples) - set(cases1))
        mu_case, mu_control = compute_mu_diff(cases1, controls1, gene_length,
                                              design_matrix, args.degenerate)
        distance_matrix[str(i)] = mu_control - mu_case
  
    # distance_matrix.to_csv(os.path.join(args.savepath, "distance_matrix.tsv"), 
    #                        sep="\t", header=True, index=True) 
    distance_matrix1 = distance_matrix.drop(columns="distance")
    mean = np.mean(distance_matrix1, axis=1)
    std = np.std(distance_matrix1, axis=1)
    distance_matrix["zscore"] = (distance_matrix["distance"] - mean) / std
    distance_matrix.to_csv(os.path.join(args.savepath, "distance_matrix.tsv"), 
                           sep="\t", header=True, index=True)
    
if __name__ == "__main__":
    args = parse_args()
    main(args)        