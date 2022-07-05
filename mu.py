import pandas as pd 
import numpy as np 
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
from vcf import *
import random


def parse_args():
    """
    Parses the arguments.
    """
    parser = argparse.ArgumentParser(description="Mu arguments")
    parser.add_argument('--VCF', nargs='?', default='./',help='Location of Cohort VCF')
    parser.add_argument('--samples', nargs='?', default='./samples.txt',help='samples file path')
    parser.add_argument('--GeneLength', nargs='?', default='./refs/gene_length.csv',help='gene length file path')
    parser.add_argument('--ref', nargs='?', default='hg38', choices=('hg19', 'hg38'), help='genome reference file')
    parser.add_argument('--savepath', nargs='?', default='./',help='save path for output')
    parser.add_argument('--maxaf', type=float, help='maximum allele frequency cutoff')
    parser.add_argument('--Ann',default='VEP', choices=('VEP', 'ANNOVAR'),help='EA annotation method')
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

    samples =pd.read_csv(args.samples, sep='\t', header=None,index_col=0)
    controls = samples[samples.iloc[:,0]==0].index.astype(str).tolist()
    cases = samples[samples.iloc[:,0]==1].index.astype(str).tolist()
    total_samples = samples.index.astype(str).tolist()
    if args.Ann=='ANNOVAR':
        matrix = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)(args.VCF, gene, ref.loc[gene], total_samples, min_af=0, max_af=args.maxaf,af_field='AF',EA_parser='canonical') for gene in tqdm(ref.index.unique()))
    if args.Ann=='VEP':
        matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP)(args.VCF, gene, ref.loc[gene], total_samples, max_af=args.maxaf, min_af=0) for gene in tqdm(ref.index.unique()))
    design_matrix = pd.concat(matrix, axis=1)
        ## reading gene length file
    gene_length = pd.read_csv(args.GeneLength, index_col=0)
    genes = set(design_matrix.columns.tolist()).intersection(set(gene_length.index.tolist()))
    gene_length = gene_length.loc[genes]
    
    design_matrix_case = design_matrix.loc[cases,genes]
    design_matrix_control = design_matrix.loc[controls,genes]
    
    SumEA_genes_case = np.sum(design_matrix_case,axis=0)
    SumEA_genes_control = np.sum(design_matrix_control,axis=0)


    ## mu calculation
    expected_energy_case = np.sum(SumEA_genes_case)/np.sum(gene_length['gene_length'])
    mu_matrix_case = pd.DataFrame(np.zeros((len(genes), 3)), index=genes, columns=['sumEA', 'gene_length', 'mu'])
    mu_matrix_case['sumEA']= SumEA_genes_case
    mu_matrix_case['gene_length'] = gene_length
    mu_matrix_case['mu'] = expected_energy_case/(mu_matrix_case['sumEA']/mu_matrix_case['gene_length'])
    
    ## mu calculation
    expected_energy_control = np.sum(SumEA_genes_control)/np.sum(gene_length['gene_length'])
    mu_matrix_control = pd.DataFrame(np.zeros((len(genes), 3)), index=genes, columns=['sumEA', 'gene_length', 'mu'])
    mu_matrix_control['sumEA']= SumEA_genes_control
    mu_matrix_control['gene_length'] = gene_length
    mu_matrix_control['mu'] = expected_energy_control/(mu_matrix_control['sumEA']/mu_matrix_control['gene_length'])
    mu_matrix = pd.DataFrame(np.zeros((len(genes), 2)), index=genes, columns=['mu_case', 'mu_control'])
    mu_matrix['mu_case'] = mu_matrix_case['mu'].copy()
    mu_matrix['mu_control'] = mu_matrix_control['mu'].copy()
    mu_matrix.to_csv(args.savepath+'mu.tsv', sep='\t', header=True, index=True) 
    
    
    distance_matrix = pd.DataFrame(np.zeros((len(genes), 1)), index=genes, columns=['distance'])
    distance_matrix['distance'] = mu_matrix_control['mu']-mu_matrix_case['mu']
    
    for i in range(1000):
        cases1 = random.sample(total_samples, len(cases))
        controls1 = list(set(total_samples) - set(cases1))
        design_matrix_case = design_matrix.loc[cases1,genes]
        design_matrix_control = design_matrix.loc[controls1,genes]
    
        SumEA_genes_case = np.sum(design_matrix_case,axis=0)
        SumEA_genes_control = np.sum(design_matrix_control,axis=0)
        
        expected_energy_case = np.sum(SumEA_genes_case)/np.sum(gene_length['gene_length'])
        mu_matrix_case['sumEA']= SumEA_genes_case
        mu_matrix_case[str(i)] = expected_energy_case/(mu_matrix_case['sumEA']/mu_matrix_case['gene_length'])
        
        expected_energy_control = np.sum(SumEA_genes_control)/np.sum(gene_length['gene_length'])
        mu_matrix_control['sumEA']= SumEA_genes_control
        mu_matrix_control[str(i)] = expected_energy_control/(mu_matrix_control['sumEA']/mu_matrix_control['gene_length'])
        
        distance_matrix[str(i)] = mu_matrix_control[str(i)]-mu_matrix_case[str(i)]
  
    distance_matrix.to_csv(args.savepath+'distance_matrix.tsv', sep='\t', header=True, index=True) 
    distance_matrix1 = distance_matrix.drop(columns='distance')
    mean = np.mean(distance_matrix1, axis=1)
    std = np.std(distance_matrix1, axis=1)
    distance_matrix['zscore'] = (distance_matrix['distance']-mean)/std
    distance_matrix.to_csv(args.savepath+'distance_matrix.tsv', sep='\t', header=True, index=True) 
    
    
    

if __name__ == "__main__":
    args = parse_args()
    main(args)        
