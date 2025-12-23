"""
Sets up MuData object for use in MOFA. 

This file combines data from the chronic lymphocytic leukaemia (CLL) into one dataset for use in MOFA. 

Source: This file is adapted from the EMBL-EBI workshop
“Multi-Omics Integration for Personalised Medicine”.
Original materials by Florin Walter (EMBL-EBI):
https://github.com/florinwalter/ebi_mofa_workshop/blob/main/project.ipynb
"""
import os
import pandas as pd
import anndata as ad
import mudata as md


def get_mdata():
    """
    Loads data from all of the different cll datasets found in data/raw/

    Output:
         mdata : mudata object containing CLL datasets as AnnData objects 

    """

    data_dict = {}

    metadata = load_metadata()

    data_dict["mrna"] = load_mrna()
    data_dict["mutations"] = load_mutations()
    data_dict["methylations"] = load_methylation()
    data_dict['drugs'] = load_drugs()

    mdata = md.MuData(data_dict)
    mdata.obs = metadata.loc[mdata.obs_names]

    return mdata
    

def load_metadata():
    """
    Loads CLL metadata
    Returns pandas df of metadata
    """
    # CLL Metadata
    metadata = pd.read_csv("data/raw/cll_metadata.csv", index_col = "Sample")

    metadata.rename(columns={"Gender" : "Sex"}, inplace=True)
    metadata.replace({"Sex" : {"m" : 0, "f" : 1}}, inplace=True)
    metadata.replace({"IGHV" : {"U" : 0, "M" : 1}}, inplace=True)
    metadata.fillna(-1, inplace=True)

    # Remove prediction variable
    metadata = metadata.drop('died', axis=1)

    return metadata

def load_mrna():
    """
    Loads CLL mRNA data.
    Returns AnnData object of mRNA dataset.
    """
    mrna = pd.read_csv("data/raw/cll_mrna.csv", index_col=0).T
    mrna_adata = ad.AnnData(mrna)

    gene_ids = pd.read_csv("data/raw/cll_geneids.csv", index_col=0)
    cols = list(mrna_adata.var_names)
    cols = [gene_ids.loc[gene_ids["GENEID"] == gene, "SYMBOL"].item() for gene in cols]
    mrna_adata.var_names = cols

    return mrna_adata


def load_mutations():
    """
    Loads CLL mutations data.
    Returns AnnData object of mutations dataset.
    """
    mutations = pd.read_csv("data/raw/cll_mutations.csv", index_col=0).T
    mutations_adata = ad.AnnData(mutations)
    mutations_adata.var_names = [f"m_{x}" for x in mutations_adata.var_names] #Need to change var names for mdata functionality
    
    return mutations_adata

def load_methylation():
    """
    Loads CLL methylation data.
    Returns AnnData object of methylation dataset.
    """
    methylation = pd.read_csv("data/raw/cll_methylation.csv", index_col=0).T
    methylations_adata = ad.AnnData(methylation)

    return methylations_adata

def load_drugs():
    """
    Loads CLL drugs data. 
    Returns AnnData object of drugs dataset.
    """
    drugs = pd.read_csv("data/raw/cll_drugs.csv", index_col=0).T
    drugs_adata = ad.AnnData(drugs)

    # assuming you called the AnnData object `adata_drugs`
    drug_names = pd.read_csv("data/raw/drugs.txt", sep=",", index_col=0)
    mapping = drug_names["name"].to_dict()
    cols = []
    for k in drugs_adata.var_names:
        for v in mapping.keys():
            if v in k:
                cols.append(k.replace(v, mapping[v]))
                break

    drugs_adata.var_names = cols

    return drugs_adata


