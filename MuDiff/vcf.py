#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import namedtuple
import numpy as np 
import pandas as pd 
from pathlib import Path
from pysam import VariantFile, VariantRecord
import re
import scipy.sparse as sp
from typing import Tuple, NamedTuple

SpliceAI = namedtuple("SpliceAI",
    "DP_AG DP_AL DP_DG DP_GL DS_AG DS_AL DS_DG DS_DL SYMBOL")


def validate_EA(ea: float) -> float:
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
        if type(ea) == str and (ea == "fs-indel" or "STOP" in ea):
            ea = 100
        else:
            ea = np.nan
    except TypeError:
        ea = np.nan
    return ea

def convert_zygo(genotype: tuple) -> int:
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
def fetch_EA_VEP(
        EA: tuple, 
        canon_ensp: tuple, 
        all_ensp: tuple, 
        csq: tuple
        ) -> float:
    """
    Return EA for VEP canonical transcript (ENSEMBL)
    """
    if "stop_gained" in csq or\
       "frameshift_variant" in csq or\
       "stop_lost" in csq:
        return 100
    elif "splice_donor_variant" in csq or "splice_acceptor_variant" in csq:# \
         #or "splice_region_variant" in csq:
        return 100
    try:
        canon_idx = all_ensp.index(canon_ensp)
    except ValueError:
        return np.nan
    else:
        return validate_EA(EA[canon_idx])

def af_check(
        rec: VariantRecord,
        min_af: float,
        max_af: float
        ) -> bool:
    """
    Check if variant allele frequency passes filters
    Args:
        rec (VariantRecord)
        min_af (float): Minimum allele frequency for variant
        max_af (float): Maximum allele frequency for variant
    Returns:
        bool: True of AF passes filters, otherwise False
    """

    af = rec.info["AF"]
    if type(af) == tuple:
        af = af[0]
    try:
        return min_af < float(af) < max_af
    except:
        return False

def _fetch_VEP_anno(anno):
    # for fields that could return either direct value 
    # or tuple depending on header
    if type(anno) == tuple:
        return anno[0]
    else:
        return anno

def validate_spliceAI(args):
    scores = [ float(s) if s else np.nan for s in args[:8] ]
    return SpliceAI(*scores, args[8])

def fetch_VEP_spliceAI(csq, canon_ensp):
    """
    Return consequence data associated to canonical transcript
    """
    def get_ensp(c):
        ensp = c.split("|")[ensp_idx]
        return ensp
    ensp_idx = 7
    spliceAI_idx = 10
    csq_canonical = next((c for c in csq if get_ensp(c) == canon_ensp), None)
    if csq_canonical:
        csq_canonical = csq_canonical.split("|")
        spliceAI = validate_spliceAI(csq_canonical[spliceAI_idx:spliceAI_idx+9])
        return spliceAI
    else:
        return None
            
def parse_VEP_synonymous(
        vcf_fn: Path, 
        gene: str, 
        gene_ref: str, 
        samples: list, 
        min_af: float,
        max_af: float,
        ) -> pd.DataFrame:
    """
    Parse EA scores and compute pEA design matrix for a given gene 
    with custom VEP annotations
    Args:
        vcf_fn (Path-like): Filepath to VCF
        gene (str): HGSC gene symbol
        gene_ref (Series): Reference information for given gene's transcripts
        samples (list): sample IDs
        min_af (float): Minimum allele frequency for variants
        max_af (float): Maximum allele frequency for variants
    Returns:
        DataFrame: sumEA design matrix
    """
   
    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    dmatrix = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, 
                           columns=[gene])
    for var in vcf:
        if re.search(r"chr", var.chrom):
            contig = "chr"+str(gene_ref.chrom)
        else:
            contig = str(gene_ref.chrom)
        break

    spliceAI_thres = 0.5

    for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
        canon_ensp = _fetch_VEP_anno(rec.info["ENSP"])
        csq = _fetch_VEP_anno(rec.info["Consequence"])
        rec_gene = _fetch_VEP_anno(rec.info["SYMBOL"])
        spliceAI = fetch_VEP_spliceAI(rec.info["CSQ"], canon_ensp)

        pass_af_check = af_check(rec, min_af, max_af)
        if pass_af_check and gene == rec_gene and "synonymous_variant" in csq \
             and (
             spliceAI.DS_AG < spliceAI_thres and \
             spliceAI.DS_AL < spliceAI_thres and \
             spliceAI.DS_DG < spliceAI_thres and \
             spliceAI.DS_DL < spliceAI_thres):
            gts = pd.Series([convert_zygo(rec.samples[sample]["GT"]) \
                             for sample in samples], index=samples, dtype=int)
            dmatrix[gene] += 100*gts  
    return dmatrix

def parse_VEP_spliceAI(
        vcf_fn: Path, 
        gene: str, 
        gene_ref: str, 
        samples: list, 
        min_af: float,
        max_af: float,
        spliceAI_thres: float = 0.5,
        ) -> pd.DataFrame:
    """
    Parse EA scores and compute pEA design matrix for a given gene 
    with custom VEP annotations
    Args:
        vcf_fn (Path-like): Filepath to VCF
        gene (str): HGSC gene symbol
        gene_ref (Series): Reference information for given gene's transcripts
        samples (list): sample IDs
        min_af (float): Minimum allele frequency for variants
        max_af (float): Maximum allele frequency for variants
    Returns:
        DataFrame: sumEA design matrix
    """
   
    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    dmatrix = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, 
                           columns=[gene])
    for var in vcf:
        if re.search(r"chr", var.chrom):
            contig = "chr"+str(gene_ref.chrom)
        else:
            contig = str(gene_ref.chrom)
        break

    for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
        canon_ensp = _fetch_VEP_anno(rec.info["ENSP"])
        csq = _fetch_VEP_anno(rec.info["Consequence"])
        rec_gene = _fetch_VEP_anno(rec.info["SYMBOL"])
        spliceAI = fetch_VEP_spliceAI(rec.info["CSQ"], canon_ensp)

        pass_af_check = af_check(rec, min_af, max_af)
        if pass_af_check and gene == rec_gene and "stop_gained" not in csq and \
            spliceAI is not None and (
            spliceAI.DS_AG >= spliceAI_thres or \
            spliceAI.DS_AL >= spliceAI_thres or \
            spliceAI.DS_DG >= spliceAI_thres or \
            spliceAI.DS_DL >= spliceAI_thres):
            gts = pd.Series([convert_zygo(rec.samples[sample]["GT"]) \
                             for sample in samples], index=samples, dtype=int)
            # print(csq, rec_gene, spliceAI, sum(gts), rec.info["AF"], 
            #       rec.pos, rec.ref, rec.alts,
            #       rec.info["AC_0"],
            #       rec.info["AC_1"],
            #       rec.info["EXON"],
            #       canon_ensp, rec.info["Ensembl_proteinid"],
            #       rec.info.get("EA", (None,))
            #       ,_fetch_VEP_anno(rec.info["Consequence"]))
            dmatrix[gene] += gts  
    return dmatrix

def parse_VEP(
        vcf_fn: Path, 
        gene: str, 
        gene_ref: str, 
        samples: list, 
        min_af: float,
        max_af: float,
        include_spliceAI: bool = True,
        spliceAI_thres: float = 0.5,
        ) -> pd.DataFrame:
    """
    Parse EA scores and compute pEA design matrix for a given gene 
    with custom VEP annotations
    Args:
        vcf_fn (Path-like): Filepath to VCF
        gene (str): HGSC gene symbol
        gene_ref (Series): Reference information for given gene's transcripts
        samples (list): sample IDs
        min_af (float): Minimum allele frequency for variants
        max_af (float): Maximum allele frequency for variants
    Returns:
        DataFrame: sumEA design matrix
    """
   
    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    dmatrix = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, 
                           columns=[gene])
    for var in vcf:
        if re.search(r"chr", var.chrom):
            contig = "chr"+str(gene_ref.chrom)
        else:
            contig = str(gene_ref.chrom)
        break

    for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
        all_ea = rec.info.get("EA", (None,))
        all_ensp = rec.info.get("Ensembl_proteinid", (rec.info["ENSP"][0],))
        canon_ensp = _fetch_VEP_anno(rec.info["ENSP"])
        csq = _fetch_VEP_anno(rec.info["Consequence"])
        rec_gene = _fetch_VEP_anno(rec.info["SYMBOL"])
        ea = fetch_EA_VEP(all_ea, canon_ensp, all_ensp, csq)
        pass_af_check = af_check(rec, min_af, max_af)

        if not np.isnan(ea).all() and gene == rec_gene and pass_af_check:
            gts = pd.Series([convert_zygo(rec.samples[sample]["GT"]) \
                             for sample in samples], index=samples, dtype=int)
            dmatrix[gene] += ea*gts

        if include_spliceAI and ea < 100:
            # SpliceAI data
            spliceAI = fetch_VEP_spliceAI(rec.info["CSQ"], canon_ensp)
            if  spliceAI is not None and (
                spliceAI.DS_AG >= spliceAI_thres or \
                spliceAI.DS_AL >= spliceAI_thres or \
                spliceAI.DS_DG >= spliceAI_thres or \
                spliceAI.DS_DL >= spliceAI_thres):
                gts = pd.Series([convert_zygo(rec.samples[sample]["GT"]) \
                                for sample in samples], 
                                index=samples, dtype=int)
                dmatrix[gene] += 100*gts  

    return dmatrix   

def parse_VEP_degenerate(
        vcf_fn: Path, 
        gene: str, 
        gene_ref: pd.Series, 
        samples: list, 
        min_af: float,
        max_af: float
    ) -> Tuple[pd.DataFrame, sp.csc_matrix]: 
    """
    Parse EA scores and compute pEA design matrix for a given gene 
    with custom VEP annotations
    ### Parameters
        - vcf_fn (Path-like): Filepath to VCF
        - gene (str): HGSC gene symbol
        - gene_ref (Series): Reference information for given gene's transcripts
        - samples (list): sample IDs
        - min_af (float): Minimum allele frequency for variants
        - max_af (float): Maximum allele frequency for variants
    ---
    ### Return
        DataFrame: sumEA design matrix
    """

    vcf = VariantFile(vcf_fn)
    vcf.subset_samples(samples)
    ea_matrix, gt_matrix = [], []

    for var in vcf:
        if re.search(r"chr", var.chrom):
            contig = "chr"+str(gene_ref.chrom)
        else:
            contig = str(gene_ref.chrom)
        break

    for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
        all_ea = rec.info.get("EA", (None,))
        all_ensp = rec.info.get("Ensembl_proteinid", (rec.info["ENSP"][0],))
        canon_ensp = _fetch_VEP_anno(rec.info["ENSP"])
        csq = _fetch_VEP_anno(rec.info["Consequence"])
        rec_gene = _fetch_VEP_anno(rec.info["SYMBOL"])
        ea = fetch_EA_VEP(all_ea, canon_ensp, all_ensp, csq)
        pass_af_check = af_check(rec, min_af, max_af)
        if not np.isnan(ea).all() and gene == rec_gene and pass_af_check:
            ea_matrix.append([gene, rec.id, ea])
            gts = [convert_zygo(rec.samples[sample]["GT"]) \
                   for sample in samples]
            gt_matrix.append(gts)

    # Avoid error when concatenating matrices
    if len(ea_matrix) == 0 and len(gt_matrix) == 0:
        return None, None

    ea_matrix = pd.DataFrame(ea_matrix, columns=["gene", "variant", "EA"])
    gt_matrix = sp.csc_matrix(gt_matrix, dtype=np.int8)
    return ea_matrix, gt_matrix   

def split_genes(rec: VariantRecord) -> VariantRecord:
    """
    If a variant has overlapping gene annotations, it will be split into 
    separate records with correct corresponding transcripts, substitutions,
    and EA scores
    ### Parameters
        rec (VariantRecord)
    ---
    ### Return
        VariantRecord
    """
    def _gene_map(gene_idxs_dict, values):
        genes = gene_idxs_dict.keys()
        if len(values) == len([i for gene_idxs in gene_idxs_dict.values() \
                               for i in gene_idxs]):
            val_d = {g: [values[i] for i in gene_idxs_dict[g]] for g in genes}
        elif len(values) == 1:
            val_d = {g: values for g in genes}
        else:
            raise ValueError("Size of values list does not match expected case.")
        return val_d

    ea = rec.info["EA"]
    gene = rec.info["gene"]
    nm = rec.info["NM"]
    geneset = set(gene)
    idxs = { genekey: [i for i, g in enumerate(gene) if g == genekey] \
            for genekey in geneset }
    ea_dict = _gene_map(idxs, ea)
    nm_dict = _gene_map(idxs, nm)
    for g in geneset:
        var = rec.copy()
        var.info["gene"] = g
        var.info["EA"] = tuple(ea_dict[g])
        var.info["NM"] = tuple(nm_dict[g])
        yield var

def fetch_variants(
        vcf: VariantFile,
        contig: str = None, 
        start: int = None,
        stop: int = None
        ) -> VariantRecord:
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
        # check for overlapping gene annotations
        if type(rec.info["gene"]) == tuple:
            for var in split_genes(rec):
                yield var
        else:
            yield rec

def fetch_EA_ANNOVAR(
        EA: tuple, 
        nm_ids: tuple,
        canon_nm: str, 
        EA_parser: str = "canonical"
        ) -> float:
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
    if EA_parser == "canonical" and canon_nm in nm_ids:
        if len(EA) > 1:
            return validate_EA(EA[nm_ids.index(canon_nm)])
        else:
            return validate_EA(EA[0])
    elif EA_parser != "canonical":
        newEA = []
        for score in EA:
            newEA.append(validate_EA(score))
        if np.isnan(newEA).all():
            return np.nan
        elif EA_parser == "mean":
            return np.nanmean(newEA)
        elif EA_parser == "max":
            return np.nanmax(newEA)
        else:
            return newEA
    else:
        return np.nan
        
def parse_ANNOVAR(
        vcf_fn: Path, 
        gene: str, 
        gene_ref: pd.Series, 
        samples: list,
        min_af: float = None, 
        max_af: float = None, 
        af_field: str = "AF", 
        EA_parser: str = "canonical"
        ) -> pd.DataFrame:
    """
    Parse EA scores and compute pEA design matrix for a given gene 
    with custom ANNOVAR annotations
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
    dmatrix = pd.DataFrame(np.ones((len(samples), 1)), index=samples, 
                           columns=[gene])
    for rec in fetch_variants(vcf, contig=str(gene_ref.chrom), 
                              start=gene_ref.start, stop=gene_ref.end):
        ea = fetch_EA_ANNOVAR(rec.info["EA"], rec.info["NM"], 
                              gene_ref.canonical, EA_parser=EA_parser)
        pass_af_check = af_check(rec, af_field=af_field, 
                                 min_af=min_af, max_af=max_af)
        if not np.isnan(ea).all() and gene == rec.info["gene"] and pass_af_check:
            gts = pd.Series([convert_zygo(rec.samples[sample]["GT"])\
                             for sample in samples], index=samples, dtype=int)
            dmatrix[gene] += ea*gts  
           
    return dmatrix 
