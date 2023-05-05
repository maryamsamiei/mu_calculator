import argparse
from joblib import Parallel, delayed
import numpy as np 
import os
import pandas as pd 
from pathlib import Path
import random
import scipy.sparse as sp
from tqdm import tqdm
from typing import Union, Tuple
from .vcf import *


def path(relative_path: str) -> Path:
    """
    Return path relative to module parent
    """
    return Path(__file__).parent / relative_path

def parse_args() -> argparse.Namespace:
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
    parser.add_argument("--GeneExons", nargs="?", 
        default=path("refs/gene_exons.tsv"), 
        help="gene number of exons file path")
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
    parser.add_argument("--method", type=str, default="EA",
        choices=("EA", "spliceAI", "EA-spliceAI", "synonymous"),
        help="number of CPUs to use for multiprocessing")
    parser.add_argument("--spliceAI_thres", type=float, default=0.5,
        help="SpliceAI threshold")
    parser.add_argument("--degenerate",default=False, action="store_true",
        help="Use number of variant sites (degenerate distribution)")
    parser.add_argument("--cores", type=int, default=1,
        help="number of CPUs to use for multiprocessing")
    return parser.parse_args()

def _SumEA_degenerate(
        samples: list,
        total_samples: list,
        ea_matrix: pd.DataFrame, 
        gt_matrix: sp.csc_matrix, 
        ) -> pd.Series:
    """
    Return degenerate SumEA based on variants in samples for every gene
    ### Parameters
    - samples: list of sample ids onto which perform analysis
    - total_samples: list of all sample ids
    - ea_matrix: Pandas DataFrame containing EA per gene variant
    - gt_matrix: Scipy.sparse csc_matrix containing genotype per sample
        - rows are variants
        - columns are samples
    ---
    ### Return
    Degenerate sumEA for each gene
    """
    # Only keep variants that appear at least once in samples
    samples = set(samples)
    sample_vector = np.array([1 if s in samples else 0 for s in total_samples])
    variant_matrix = gt_matrix.multiply(sample_vector)
    sample_variants = variant_matrix.sum(axis=1, dtype=np.int8).A1 > 0
    dmatrix_sample = ea_matrix[["gene", "EA"]].loc[sample_variants]
    return dmatrix_sample.groupby("gene").EA.sum()

def compute_mu_diff(
        cases: list, 
        controls: list,
        samples: list,
        gene_length: pd.DataFrame,
        design_matrix: pd.DataFrame,
        ea_matrix: pd.DataFrame,
        gt_matrix: sp.csc_matrix,
        degenerate: bool,
        ) -> Tuple[pd.Series, pd.Series]:
    """
    Compute Mu-diff for a given set of cases and controls
    ### Parameters
    - cases : list of case ids
    - controls: list of control ids
    - samples: list of sample ids
    - gene_length: Pandas DataFrame containing gene length
        - index: gene name (only those used for analysis)
        - gene_length: int
    - design_matrix: Pandas DataFrame containing pEA per gene per sample
        - rows are samples
        - columns are genes
    - ea_matrix: Pandas DataFrame containing EA per gene variant
    - gt_matrix: Scipy.sparse csc_matrix containing genotype per sample
        - rows are variants
        - columns are samples
    ----
    ### Return
    Tuple containing pd.Series for mu-cases and mu-controls
    """
    # Subset sumEA matrix to genes and samples of study (cases and controls)
    genes = gene_length.index.to_list()
    if degenerate:
        SumEA_genes_case = _SumEA_degenerate(cases, samples, 
                                             ea_matrix, gt_matrix)\
                                            .reindex(genes)
        SumEA_genes_control = _SumEA_degenerate(controls, samples, 
                                                ea_matrix, gt_matrix)\
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

def compute_dmatrix(
        ref: pd.DataFrame,
        samples: list, 
        args: argparse.Namespace
        ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, sp.csc_matrix]]:
    """
    Returns matrices needed for Mu computation
        - For non-degenerate, it returns matrix containing pEA scores 
          for each gene in each sample
            rows: samples
            columns: genes
        - For degenerate, it returns two matrices:
            - ea_matrix: variant (rows) EA scores and associated genes
                rows: variants
                columns: gene, variant (VariantRecord.id), EA
            - gt_matrix: genotype for each sample for each studied variant
                rows: variants
                columns: samples
    """
    # Build SumEA matrix (sample in rows, genes in columns)
    tqdm_desc = "EA matrix"
    if args.Ann=="ANNOVAR":
        matrix = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)\
            (args.VCF, gene, ref.loc[gene], samples, min_af=0,
             max_af=args.maxaf, af_field="AF", EA_parser="canonical")\
                for gene in tqdm(ref.index.unique(), desc=tqdm_desc))

    if args.Ann=="VEP":
        if args.degenerate:
            matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP_degenerate)\
                (args.VCF, gene, ref.loc[gene], samples, 
                 min_af=0, max_af=args.maxaf) \
                    for gene in tqdm(ref.index.unique(), desc=tqdm_desc))
            ea_matrix, gt_matrix = list(zip(*matrix))
            ea_matrix = [ m for m in ea_matrix if m is not None ]
            gt_matrix = [ m for m in gt_matrix if m is not None ]
            assert len(ea_matrix) == len(gt_matrix), "EA and GT matrix mismatch"
            ea_matrix = pd.concat(ea_matrix, axis=0, ignore_index=True)
            gt_matrix = sp.vstack(gt_matrix, format="csc", dtype=np.int8)
            ea_n = ea_matrix.shape[0]
            gt_n = gt_matrix.get_shape()[0]
            assert ea_n == gt_n,\
                "EA and GT matrix mismatch number of variants {ea_n} != {gt_n}"
            return ea_matrix, gt_matrix
        else:
            tqdm_desc = f"{args.method} matrix"
            if args.method in ("EA", "EA-spliceAI"):
                spliceAI = args.method == "EA-spliceAI"
                print(args.method, spliceAI, args.spliceAI_thres)
                matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP)\
                    (args.VCF, gene, ref.loc[gene], samples, 
                    min_af=0, max_af=args.maxaf, include_spliceAI=spliceAI,
                    spliceAI_thres=args.spliceAI_thres) \
                        for gene in tqdm(ref.index.unique(), desc=tqdm_desc))
            elif args.method == "spliceAI":
                matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP_spliceAI)\
                    (args.VCF, gene, ref.loc[gene], samples, 
                    min_af=0, max_af=args.maxaf,
                    spliceAI_thres=args.spliceAI_thres) \
                        for gene in tqdm(ref.index.unique(), desc=tqdm_desc))
            elif args.method == "synonymous":
                matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP_synonymous)\
                    (args.VCF, gene, ref.loc[gene], samples, 
                    min_af=0, max_af=args.maxaf) \
                        for gene in tqdm(ref.index.unique(), desc=tqdm_desc))
    return pd.concat(matrix, axis=1)

def main(args: argparse.Namespace) -> None:
    """
    Run mu-diff analysis given CLI arguments
    """

    # Create output directory if non existant
    os.makedirs(args.savepath, exist_ok=True)

    if args.Ann == "ANNOVAR":
        if args.ref == "hg19":
            ref = pd.read_csv(path("refs/refGene-lite_hg19.May2013.txt"), 
                              delimiter="\t", header=0, index_col="gene")
        elif args.ref == "hg38":
            ref = pd.read_csv(path("refs/refGene-lite_hg38.June2017.txt"), 
                              delimiter="\t", header=0, index_col="gene")
    elif args.Ann == "VEP":
        if args.ref == "hg19":
            ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh37.v75.txt"), 
                              delimiter="\t", header=0, index_col="gene")
        elif args.ref == "hg38":
            if args.noX:
                ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh38.v94.noX.txt"), 
                                  delimiter="\t", header=0, index_col="gene")
            else:
                ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh38.v94.txt"), 
                                  delimiter="\t", header=0, index_col="gene")


    # ref = pd.read_csv("anln.txt", delimiter="\t", header=0, index_col="gene")

    samples = pd.read_csv(args.samples, header=None, index_col=0, sep="\t")
    controls = samples[samples.iloc[:,0]==0].index.astype(str).tolist()
    cases = samples[samples.iloc[:,0]==1].index.astype(str).tolist()
    total_samples = samples.index.astype(str).tolist()

    design_matrix, ea_matrix, gt_matrix = None, None, None
    design_matrix_path = os.path.join(args.savepath, "design_matrix.tsv")
    try:
        design_matrix = pd.read_csv(design_matrix_path, sep="\t", index_col=0)
        design_matrix.index = design_matrix.index.astype(str)
        matrix_genes = design_matrix.columns.tolist()
    except:
        if args.degenerate:
            ea_matrix, gt_matrix = compute_dmatrix(ref, total_samples, args)
            matrix_genes = ea_matrix.gene.unique().tolist()
        else:
            design_matrix = compute_dmatrix(ref, total_samples, args)
            matrix_genes = design_matrix.columns.tolist()

        design_matrix.to_csv(design_matrix_path, sep="\t")

    ## reading gene length file
    if args.method == "spliceAI":
        gene_length = pd.read_csv(args.GeneExons, index_col=0, sep="\t")
    else:
        gene_length = pd.read_csv(args.GeneLength, index_col=0)

    genes = set(matrix_genes)\
        .intersection(set(gene_length.index.tolist()))
    gene_length = gene_length.loc[genes]
    
    mu_case, mu_control = compute_mu_diff(cases, 
                                          controls, 
                                          total_samples,
                                          gene_length,
                                          design_matrix, 
                                          ea_matrix,
                                          gt_matrix,
                                          args.degenerate)

    mu_matrix = pd.DataFrame(np.zeros((len(genes), 2)), index=genes, 
                             columns=["mu_case", "mu_control"])
    mu_matrix["mu_case"] = mu_case.copy()
    mu_matrix["mu_control"] = mu_control.copy()
    mu_matrix.to_csv(os.path.join(args.savepath, "mu.tsv"), sep="\t", 
                     header=True, index=True) 
    
    distance_matrix = pd.DataFrame(np.zeros((len(genes), 1)), index=genes, 
                                   columns=["distance"])
    distance_matrix["distance"] = mu_control - mu_case
    
    for i in tqdm(range(1000), desc="Randomization"):
        cases1 = random.sample(total_samples, len(cases))
        controls1 = list(set(total_samples) - set(cases1))

        mu_case, mu_control = compute_mu_diff(cases1, 
                                              controls1, 
                                              total_samples,
                                              gene_length,
                                              design_matrix, 
                                              ea_matrix,
                                              gt_matrix,
                                              args.degenerate)

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
