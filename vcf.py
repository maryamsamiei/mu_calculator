#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd 
import numpy as np 
from pysam import VariantFile


def validate_EA(ea):
    """
    Checks for valid EA score
    Args:
        ea (str/float/None): EA score as string
    Returns:
        float: EA score between 0-100 if valid, otherwise returns NaN
    """
    try:
        ea = float(ea)
    except ValueError:
        if type(ea) == str and (ea == 'fs-indel' or 'STOP' in ea):
            ea = 100
        else:
            ea = np.nan
    except TypeError:
        ea = np.nan
    return ea


def convert_zygo(genotype):
    """
    Convert a genotype tuple to a zygosity integer
    Args:
        genotype (tuple): The genotype of a variant for a sample
    Returns:
        int: The zygosity of the variant (0/1/2)
    """
    if genotype in [(1, 0), (0, 1)]:
        zygo = 1
    elif genotype == (1, 1):
        zygo = 2
    else:
        zygo = 0
    return zygo
### functions for VEP annottated vcf
def fetch_EA_VEP(EA, canon_ensp, all_ensp, csq):
    if 'stop_gained' in csq or 'frameshift_variant' in csq or 'stop_lost' in csq:# or 'splice_donor_variant' in csq or 'splice_acceptor_variant' in csq:
        return 100
    try:
        canon_idx = all_ensp.index(canon_ensp)
    except ValueError:
        return np.nan
    else:
        return validate_EA(EA[canon_idx])

def af_check(rec, max_af, min_af):
    """
    Check if variant allele frequency passes filters
    Args:
        rec (VariantRecord)
        af_field (str): Name of INFO field containing allele frequency information
        max_af (float): Maximum allele frequency for variant
        min_af (float): Minimum allele frequency for variant
    Returns:
        bool: True of AF passes filters, otherwise False
    """

    af = rec.info['AF']
    if type(af) == tuple:
        af = af[0]
    try:
        return min_af < float(af) < max_af
    except:
        return False


def parse_VEP(vcf_fn, gene, gene_ref, samples, max_af, min_af):
   
    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    dmatrix = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, columns=[gene])
    
    def _fetch_anno(anno):
        # for fields that could return either direct value or tuple depending on header
        if type(anno) == tuple:
            return anno[0]
        else:
            return anno
    for rec in vcf.fetch(contig=str(gene_ref.chrom), start=gene_ref.start, stop=gene_ref.end):
        all_ea = rec.info.get('EA', (None,))
        all_ensp = rec.info.get('Ensembl_proteinid', (rec.info['ENSP'][0],))
        canon_ensp = _fetch_anno(rec.info['ENSP'])
        csq = _fetch_anno(rec.info['Consequence'])
        rec_gene = _fetch_anno(rec.info['SYMBOL'])
        ea = fetch_EA_VEP(all_ea, canon_ensp, all_ensp, csq)
        pass_af_check = af_check(rec, max_af, min_af)
        if not np.isnan(ea).all() and gene == rec_gene and pass_af_check:
            gts = pd.Series([convert_zygo(rec.samples[sample]['GT']) for sample in samples], index=samples, dtype=int)
            dmatrix[gene] += ea*gts  
    return dmatrix   
def split_genes(rec):
    """
    If a variant has overlapping gene annotations, it will be split into separate records with correct corresponding
    transcripts, substitutions, and EA scores
    Args:
        rec (VariantRecord)
    Yields:
        VariantRecord
    """
    def _gene_map(gene_idxs_dict, values):
        genes = gene_idxs_dict.keys()
        if len(values) == len([i for gene_idxs in gene_idxs_dict.values() for i in gene_idxs]):
            val_d = {g: [values[i] for i in gene_idxs_dict[g]] for g in genes}
        elif len(values) == 1:
            val_d = {g: values for g in genes}
        else:
            raise ValueError('Size of values list does not match expected case.')
        return val_d

    ea = rec.info['EA']
    gene = rec.info['gene']
    nm = rec.info['NM']
    geneset = set(gene)
    idxs = {genekey: [i for i, g in enumerate(gene) if g == genekey] for genekey in geneset}
    ea_dict = _gene_map(idxs, ea)
    nm_dict = _gene_map(idxs, nm)
    for g in geneset:
        var = rec.copy()
        var.info['gene'] = g
        var.info['EA'] = tuple(ea_dict[g])
        var.info['NM'] = tuple(nm_dict[g])
        yield var

def fetch_variants(vcf, contig=None, start=None, stop=None):
    """
    Variant iterator
    Args:
        vcf (VariantFile): pysam VariantFile
        contig (str): Chromosome of interest
        start (int): The starting nucleotide position of the region of interest
        stop (int): The ending nucleotide position of the region of interest
    Yields:
        VariantRecord
    """
    for rec in vcf.fetch(contig=contig, start=start, stop=stop):
        if type(rec.info['gene']) == tuple:  # check for overlapping gene annotations
            for var in split_genes(rec):
                yield var
        else:
            yield rec
def fetch_EA_ANNOVAR(EA, nm_ids, canon_nm, EA_parser='canonical'):
    """
    Parse EA scores for a given variant
    Args:
        EA (tuple): The EA scores parsed from a variant
        nm_ids (tuple): Transcript IDs corresponding to EA scores
        canon_nm (str): Canonical NM ID based on smallest numbered transcript
        EA_parser (str): How to aggregate multiple transcript scores
    Returns:
        float/list: Valid EA scores, refactored as floats
    Note: EA must be string type
    """
    if EA_parser == 'canonical' and canon_nm in nm_ids:
        if len(EA) > 1:
            return validate_EA(EA[nm_ids.index(canon_nm)])
        else:
            return validate_EA(EA[0])
    elif EA_parser != 'canonical':
        newEA = []
        for score in EA:
            newEA.append(validate_EA(score))
        if np.isnan(newEA).all():
            return np.nan
        elif EA_parser == 'mean':
            return np.nanmean(newEA)
        elif EA_parser == 'max':
            return np.nanmax(newEA)
        else:
            return newEA
    else:
        return np.nan
        
def parse_ANNOVAR(vcf_fn, gene, gene_ref, samples, min_af=None, max_af=None, af_field='AF', EA_parser='canonical'):
    """
    Parse EA scores and compute pEA design matrix for a given gene with custom ANNOVAR annotations
    Args:
        vcf_fn (Path-like): Filepath to VCF
        gene (str): HGSC gene symbol
        gene_ref (Series): Reference information for given gene's transcripts
        samples (list): sample IDs
        min_af (float): Minimum allele frequency for variants
        max_af (float): Maximum allele frequency for variants
        af_field (str): Name of INFO field containing allele frequency information
        EA_parser (str): How to parse EA scores from multiple transcripts
    Returns:
        DataFrame: sumEA design matrix
    """

    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    dmatrix = pd.DataFrame(np.ones((len(samples), 1)), index=samples, columns=[gene])
    for rec in fetch_variants(vcf, contig=str(gene_ref.chrom), start=gene_ref.start, stop=gene_ref.end):
        ea = fetch_EA_ANNOVAR(rec.info['EA'], rec.info['NM'], gene_ref.canonical, EA_parser=EA_parser)
        pass_af_check = af_check(rec, af_field=af_field, max_af=max_af, min_af=min_af)
        if not np.isnan(ea).all() and gene == rec.info['gene'] and pass_af_check:
            gts = pd.Series([convert_zygo(rec.samples[sample]['GT']) for sample in samples], index=samples, dtype=int)
            dmatrix[gene] += ea*gts  
           
    return dmatrix 

