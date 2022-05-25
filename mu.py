#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 14:52:38 2022

@author: maryamsamieinasab
"""

import pandas as pd 
import numpy as np 
from pysam import VariantFile
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
from vcf import *


def parse_args():
    """
    Parses the arguments.
    """
    parser = argparse.ArgumentParser(description="Big Pipeline arguments")
    parser.add_argument('--VCF', nargs='?', default='./',help='Location of Cohort VCF')
    parser.add_argument('--samples', nargs='?', default='./samples.txt',help='samples file path')
    parser.add_argument('--GeneLength', nargs='?', default='./refs/gene_length.csv',help='gene length file path')
    parser.add_argument('--ref', nargs='?', default='hg38', choices=('hg19', 'hg38'), help='genome reference file')
    parser.add_argument('--savepath', nargs='?', default='./',help='save path for output')
    parser.add_argument('--maxaf', type=float, help='maximum allele frequency cutoff')
    parser.add_argument('--Ann', default='VEP', choices=('VEP', 'ANNOVAR'),help='EA annotation method')
    parser.add_argument('--cores', type=int, default=1, help='number of CPUs to use for multiprocessing')


    return parser.parse_args()

def main(args):
    if args.Ann=='ANNOVAR':
        if args.ref=='hg19':
            ref = pd.read_csv('./refs/refGene-lite_hg19.May2013.txt', delimiter='\t', header=0, index_col='gene')
        elif args.ref=='hg38':
            ref = pd.read_csv('./refs/refGene-lite_hg38.June2017.txt', delimiter='\t', header=0, index_col='gene')
    elif args.Ann=='VEP':
        if args.ref=='hg19':
            ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh37.v75.txt', delimiter='\t', header=0, index_col='gene')
        elif args.ref=='hg38':
            ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.txt', delimiter='\t', header=0, index_col='gene')



    samples =pd.read_csv(args.samples, sep='\t', header=None,index_col=0).index.astype(str).tolist()
    if args.Ann=='ANNOVAR':
        matrix = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)(args.VCF, gene, ref.loc[gene], samples, min_af=0, max_af=args.maxaf,af_field='AF', 
                                                    EA_parser='canonical') for gene in tqdm(ref.index.unique()))
    if args.Ann=='VEP':
        matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP)(args.VCF, gene, ref.loc[gene], samples, max_af=args.maxaf, min_af=0) for gene in tqdm(ref.index.unique()))
    design_matrix = pd.concat(matrix, axis=1)
    SumEA_genes = np.sum(design_matrix,axis=0)
    ## reading gene length file
    gene_length = pd.read_csv(args.GeneLength, index_col=0)

    genes = set(SumEA_genes.index.tolist()).intersection(set(gene_length.index.tolist()))
    SumEA_genes = SumEA_genes.loc[genes]
    gene_length = gene_length.loc[genes]
    
    ## mu calculation
    expected_energy = np.sum(SumEA_genes)/np.sum(gene_length['gene_length'])
    mu_matrix = pd.DataFrame(np.zeros((len(genes), 3)), index=genes, columns=['sumEA', 'gene_length', 'mu'])
    mu_matrix['sumEA']= SumEA_genes
    mu_matrix['gene_length'] = gene_length
    mu_matrix['mu'] = expected_energy/(mu_matrix['sumEA']/mu_matrix['gene_length'])
    
    mu_matrix.to_csv(args.savepath+'mu.tsv', sep='\t', header=True, index=True) 
    
    

if __name__ == "__main__":
    args = parse_args()
    main(args)        
