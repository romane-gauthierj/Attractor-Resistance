# Cell model passport - RNA seq expression

import numpy as np
import pandas as pd
import os
import diptest

from sklearn.mixture import GaussianMixture
from scipy import stats
from scipy.stats import kurtosis
from scipy.stats import chi2

import matplotlib.pyplot as plt

from sklearn.model_selection import KFold
import logging

logger = logging.getLogger(__name__)

### =========  Assess the gene's distribution ==========  
# Bimodal patterns: (Hartigan's Dip test < 0.05, Bimodality Index > 1.5, and Kurtosis < 1) 

# BI helps determine whether a dataset has two clear subpopulations (modes) 
# that can be statistically modeled as two bell-shaped curves (normal distributions) 
# with the same spread but different centers.

def identify_genes_distribution(rna_seq_data_filtered_analysis_pivoted):
    """
    Classifies each gene as 'bimodal', 'zero_infl', or 'unimodal' based on distribution tests.
    Adds a 'genes_distribution' column to the input DataFrame.
    """

    
    rna_seq_data_filtered_analysis_pivoted = rna_seq_data_filtered_analysis_pivoted.drop('Entrez_Gene_Id', axis=1, errors='ignore')
    
    # remove duplicates

    rna_seq_data_filtered_analysis_pivoted = (
    rna_seq_data_filtered_analysis_pivoted.groupby(rna_seq_data_filtered_analysis_pivoted.index)
    .mean()     
    )



    def hartigan_dip_test(gene_data):
        dip, pval = diptest.diptest(gene_data)
        return dip,pval 


    # Bimodal Index function (equivalent to R's bimodalIndex)
    def bimodal_index(data):
        """
        Calculate bimodality index using Gaussian Mixture Models
        BI = √[π(1-π)] × δ where π is mixing proportion and δ is standardized distance
        """
        # Remove NaN values
        data = data[~np.isnan(data)]
        # Reshape for sklearn
        data_reshaped = data.reshape(-1, 1)
        
        # Fit 1-component and 2-component Gaussian mixture models
        gmm1 = GaussianMixture(n_components=1, random_state=42)
        gmm2 = GaussianMixture(n_components=2, random_state=42)

        gmm1.fit(data_reshaped)
        gmm2.fit(data_reshaped)
            
        # Calculate bimodality index using 2-component model parameters
        mu1, mu2 = gmm2.means_.flatten()
        sigma = np.sqrt(gmm2.covariances_.mean())
        pi = gmm2.weights_[0]
        
        delta = abs(mu1 - mu2) / sigma
        bi = np.sqrt(pi * (1 - pi)) * delta
        return bi


    def compute_bimodial_genes(genes_data_filtered):
    # identification of bimodial genes - Hartigan's Dip test < 0.05, Bimodality Index > 1.5, and Kurtosis < 1
    # identification of bimodal genes - Hartigan's Dip test < 0.05, Bimodality Index > 1.5, and Kurtosis < 1
        genes_names = list(set(genes_data_filtered.index))
        results_bimodal = pd.DataFrame(index=genes_names, columns = ['dip', 'pval_dip_test', 'BI', 'Kurtosis'])
        for gene in genes_names:
            gene_data = genes_data_filtered.loc[gene].values  
            dip,pval  = hartigan_dip_test(gene_data)
            results_bimodal.loc[gene, 'dip'] = dip
            results_bimodal.loc[gene, 'pval_dip_test'] = pval

            bi = bimodal_index(gene_data)
            results_bimodal.loc[gene, 'BI'] = bi

            kurt = stats.kurtosis(gene_data, fisher=True)
            results_bimodal.loc[gene, 'Kurtosis'] = kurt

            results_genes_bimodal = results_bimodal[
            (results_bimodal['pval_dip_test'] < 0.05) & 
            (results_bimodal['BI'] > 1.5) 
            & (results_bimodal['Kurtosis'] < 1)
        ]
            
        # if results_genes_bimodal.empty:
        #     bimodal_genes = []
        # else:
        #     
        bimodal_genes = list(set(results_genes_bimodal.index))

        return bimodal_genes


    
    def assess_zero_inflation_comprehensive(genes_data_filtered):
        """
        Comprehensive zero-inflation assessment using multiple approaches
        """
        genes_names = list(genes_data_filtered.index)
        results_zero_inflation = pd.DataFrame(
            index=genes_names, 
            columns=['proportion_zero', 'expected_zero_poisson', 'expected_zero_nb', 
                    'excess_zeros_poisson', 'excess_zeros_nb', 'zi_score', 'is_zero_inflated']
        )
        
        for gene in genes_names:
            gene_data = genes_data_filtered.loc[gene].values
            gene_data = gene_data[~np.isnan(gene_data)]  # Remove NaN
            
            n_total = len(gene_data)
            n_zeros = np.sum(gene_data == 0)  # Exact zeros
            prop_zeros = n_zeros / n_total
            
            # For RNA-seq, also consider very low expression as "zeros"
            n_low_expr = np.sum(gene_data < 1)  # TPM < 1
            prop_low_expr = n_low_expr / n_total
            
            results_zero_inflation.loc[gene, 'proportion_zero'] = prop_zeros
            results_zero_inflation.loc[gene, 'proportion_low_expr'] = prop_low_expr
            
            # Only analyze non-zero values for distribution fitting
            non_zero_data = gene_data[gene_data > 0]
            
            if len(non_zero_data) > 5:  # Need enough data points
                # Fit Poisson distribution to non-zero data
                lambda_poisson = np.mean(non_zero_data)
                expected_zero_poisson = stats.poisson.pmf(0, lambda_poisson)
                
                # Fit Negative Binomial to non-zero data
                # Convert to discrete for NB fitting
                discrete_data = np.round(non_zero_data).astype(int)
                try:
                    # Estimate NB parameters
                    mean_data = np.mean(discrete_data)
                    var_data = np.var(discrete_data)
                    
                    if var_data > mean_data:  # Overdispersed
                        p = mean_data / var_data
                        r = mean_data * p / (1 - p)
                        expected_zero_nb = stats.nbinom.pmf(0, r, p)
                    else:
                        expected_zero_nb = expected_zero_poisson
                        
                except:
                    expected_zero_nb = expected_zero_poisson
                
                # Calculate excess zeros
                excess_zeros_poisson = prop_zeros - expected_zero_poisson
                excess_zeros_nb = prop_zeros - expected_zero_nb
                
                # Zero-inflation score (how much more zeros than expected)
                zi_score = excess_zeros_poisson / (expected_zero_poisson + 1e-10)
                
                results_zero_inflation.loc[gene, 'expected_zero_poisson'] = expected_zero_poisson
                results_zero_inflation.loc[gene, 'expected_zero_nb'] = expected_zero_nb
                results_zero_inflation.loc[gene, 'excess_zeros_poisson'] = excess_zeros_poisson
                results_zero_inflation.loc[gene, 'excess_zeros_nb'] = excess_zeros_nb
                results_zero_inflation.loc[gene, 'zi_score'] = zi_score
                
                # Criteria for zero-inflation
                is_zero_inflated = (
                    (prop_zeros > 0.3) &  # At least 30% zeros
                    (excess_zeros_poisson > 0.1) &  # Significant excess
                    (zi_score > 0.5)  # 50% more zeros than expected
                )
                

                results_zero_inflation.loc[gene, 'is_zero_inflated'] = is_zero_inflated

                # list_zero_infla_genes = list(results_zero_inflation[
                #     results_zero_inflation['is_zero_inflated'] == True
                # ].index)
        if 'is_zero_inflated' in results_zero_inflation.columns:
            list_zero_infla_genes = list(
            results_zero_inflation[results_zero_inflation['is_zero_inflated'].fillna(False) == True].index
        )
        else:
            list_zero_infla_genes = []
        
        return list_zero_infla_genes
    

    bimodal_genes = compute_bimodial_genes(rna_seq_data_filtered_analysis_pivoted)
    zero_infl_genes = assess_zero_inflation_comprehensive(rna_seq_data_filtered_analysis_pivoted)

    def assign_distribution(gene):
        if gene in bimodal_genes:
            return 'bimodal'
        elif gene in zero_infl_genes:
            return 'zero_infl'
        else:
            return 'unimodal'



    rna_seq_data_filtered_analysis_pivoted['genes_distribution'] = [
        assign_distribution(gene) for gene in rna_seq_data_filtered_analysis_pivoted.index
    ]

    return rna_seq_data_filtered_analysis_pivoted









## Compute the normalization of each gene 
def compute_multi_distrib_normalization(expression_genes):
    # expect expression_genes with index the genes names (Hugo_Symbol) and colnames the patients id
    
    # expression_genes needs to have the column for the genes_distribution (unimodal, bimodal, zero_inf)

    # --- Zero-inflated genes normalization ---
    def compute_zero_inflat_genes_normalization(genes_expression, list_zero_infla_genes_test):

        patients_ids = genes_expression.columns.drop('genes_distribution')


        data_zero_infla_genes_normalized = pd.DataFrame(index=list_zero_infla_genes_test, columns=patients_ids)

        for gene in list_zero_infla_genes_test:
            # Get the expression values for the gene (as a 1D array, not reshaped)
            gene_expression = genes_expression.loc[gene].drop('genes_distribution').values.astype(float)
        

            q1 = np.quantile(gene_expression, 0.01)
            q99 = np.quantile(gene_expression, 0.99)


            interval = q99 - q1
            if interval > 0:
                genes_exp_normalized = (gene_expression - q1) / interval
                genes_exp_normalized[genes_exp_normalized < 0] = 0
                genes_exp_normalized[genes_exp_normalized > 1] = 1

            else:
                genes_exp_normalized = np.full_like(gene_expression, 0.5)  # If all values are the same

            # Assign normalized values to the DataFrame
            data_zero_infla_genes_normalized.loc[gene, patients_ids] = genes_exp_normalized
        results_proba_zero_inflated_genes_long = (
        data_zero_infla_genes_normalized
        .reset_index()
        .melt(id_vars='index', var_name='model_id', value_name='rsem_tpm_normalized')
        .rename(columns={'index': 'gene_symbol'})
    )
        return results_proba_zero_inflated_genes_long



    # unimodal genes -> sigmoid normalization 

    # --- Unimodal genes normalization ---
    def compute_unimodal_genes_normalization(group):
            """
            Apply sigmoid normalization to a group (gene) -> when data is unimodal
            1. Subtract median from all values
            2. Calculate lambda using MAD (Median Absolute Deviation)
            3. Apply sigmoid function: 1 / (1 + exp(-λ * x))
            """
            # Step 1: Calculate median and subtract from all values
            median_value = group['rsem_tpm'].median()
            centered_values = group['rsem_tpm'] - median_value
            
            # Step 2: Calculate MAD (Median Absolute Deviation)
            # MAD = median(|xi - median(X)|)
            mad_value = np.median(np.abs(centered_values))
            
            # Handle edge case where MAD is 0 (all values are identical)
            if mad_value == 0:
                # If all values are the same, set them all to 0.5 (median mapping)
                sigmoid_values = np.full(len(group), 0.5)
            else:
                # Step 3: Calculate lambda (λ) = log(3) / MAD
                lambda_param = np.log(3) / mad_value
                
                # Step 4: Apply sigmoid function with lambda scaling
                # sigmoid(x) = 1 / (1 + exp(-λ * x))
                sigmoid_values = 1 / (1 + np.exp(-lambda_param * centered_values))
            
            # Store the normalized values
            group['rsem_tpm_normalized'] = sigmoid_values
            
            return group


    # --- Bimodal genes normalization ---
    def compute_bimodal_genes_gaussian_mixture_model(genes_expression, bimodal_genes):
    # genes_expression: has genes as index names (Hugo_Symbol),and patients id as columns 
    #  
        patients_ids = genes_expression.columns.drop('genes_distribution')

        results_proba_bimodal_genes = pd.DataFrame(index = bimodal_genes, columns=patients_ids)
        for bimodal_gene in bimodal_genes:
            gmm = GaussianMixture(n_components=2, covariance_type='full', random_state=0)
            gene_expression = genes_expression.loc[bimodal_gene, patients_ids].values.reshape(-1, 1)
            gmm.fit(gene_expression)
            probabilities = gmm.predict_proba(gene_expression)

            if gmm.means_[0][0] < gmm.means_[1][0]:
                # Column 0 is M0, Column 1 is M1
                p_M0 = probabilities[:, 0]
                p_M1 = probabilities[:, 1]
            else:
                # Column 1 is M0, Column 0 is M1
                p_M0 = probabilities[:, 1]
                p_M1 = probabilities[:, 0]
            results_proba_bimodal_genes.loc[bimodal_gene, patients_ids] = p_M1

        results_proba_bimodal_genes_long = (
        results_proba_bimodal_genes
        .reset_index()
        .melt(id_vars='index', var_name='model_id', value_name='rsem_tpm_normalized')
        .rename(columns={'index': 'gene_symbol'})
    )

        return results_proba_bimodal_genes_long
    

    # --- Split by distribution type ---
    zero_infl_genes = expression_genes[expression_genes['genes_distribution'] == 'zero_infl']
    bimodal_genes = expression_genes[expression_genes['genes_distribution'] == 'bimodal']
    unimodal_genes = expression_genes[expression_genes['genes_distribution'] == 'unimodal']

    zero_infl_list = zero_infl_genes.index.tolist()
    bimodal_list = bimodal_genes.index.tolist()

    # --- Apply normalization for each group ---
    if len(zero_infl_list) != 0:
        results_zero = compute_zero_inflat_genes_normalization(zero_infl_genes, zero_infl_list)
    else:
        logger.warning('no genes with zero distribution')

        results_zero = pd.DataFrame()

    if len(bimodal_list) != 0:
        results_bimodal = compute_bimodal_genes_gaussian_mixture_model(bimodal_genes, bimodal_list)
    else:
        logger.warning('no genes with bimodal distribution')
        results_bimodal = pd.DataFrame()

    # Unimodal genes: convert to long format and apply sigmoid normalization
    unimodal_long = (
        unimodal_genes
        .drop(columns=['genes_distribution'])
        .reset_index()
        .melt(id_vars='Hugo_Symbol', var_name='model_id', value_name='rsem_tpm')
        .rename(columns={'Hugo_Symbol': 'gene_symbol'})
    )

    results_unimodal = (
        unimodal_long
        .groupby('gene_symbol', group_keys=False)
        .apply(compute_unimodal_genes_normalization)
    )
    results_unimodal = results_unimodal[['gene_symbol', 'model_id', 'rsem_tpm_normalized']]

    # --- Combine all results ---
    combined_results = pd.concat([results_zero, results_bimodal, results_unimodal], ignore_index=True)
    return combined_results

    










def process_genes(patients_ids, rna_seq_data, all_montagud_nodes, synonyms_to_nodes_dict):
    rna_seq_data = rna_seq_data[rna_seq_data["model_id"].isin(patients_ids)]
    rna_seq_data_filtered = rna_seq_data[rna_seq_data['gene_symbol'].isin(all_montagud_nodes)]


  # Special cases handling: Create both mTORC1 and mTORC2 from MTOR data
    syn_dict = {'MTOR': ['mTORC1', 'mTORC2'], 'MYC': ['MYC', 'MYC_MAX'], 'PIK3CA': ['PI3K', 'PIP3'], 'LDHA': ['LDHA', 'Lactic_acid'], 'ERG': ['AR_ERG', 'ERG']}
    # Get all keys
    list_genes_duplicates = syn_dict.keys()

    # apply the synonym mapping 
    for gene in rna_seq_data_filtered['gene_symbol'].unique():
        if gene in synonyms_to_nodes_dict and gene not in list_genes_duplicates:
            rna_seq_data_filtered.loc[rna_seq_data_filtered['gene_symbol'] == gene, 'gene_symbol'] = synonyms_to_nodes_dict[gene]


    for duplicate_gene in list_genes_duplicates:
        gene_duplicate_data = rna_seq_data_filtered[rna_seq_data_filtered['gene_symbol'] == duplicate_gene].copy()
        if not gene_duplicate_data.empty:
            # Create mTORC1 data
            gene_duplicate_1_data = gene_duplicate_data.copy()
            gene_duplicate_1_data['gene_symbol'] = syn_dict[duplicate_gene][0]
            
            # Create mTORC2 data  
            gene_duplicate_2_data = gene_duplicate_data.copy()
            gene_duplicate_2_data['gene_symbol'] = syn_dict[duplicate_gene][1]
            
            # Remove original MTOR and add both complexes
            rna_seq_data_filtered = rna_seq_data_filtered[rna_seq_data_filtered['gene_symbol'] != duplicate_gene]
            rna_seq_data_filtered = pd.concat([rna_seq_data_filtered, gene_duplicate_1_data, gene_duplicate_2_data], ignore_index=True)
        

    rna_seq_data_models_filtered = rna_seq_data_filtered[["model_id", "gene_symbol", "rsem_tpm"]]
    return rna_seq_data_models_filtered



# if combined with proteins

def process_genes_proteins(df_melted_protein,rna_seq_data_models_filtered):
    df_melted_protein = df_melted_protein.groupby(['model_id', 'protein_symbol'], as_index=False).agg({'rsem_tpm': 'mean'})
    protein_pairs = set(zip(df_melted_protein['model_id'], df_melted_protein['protein_symbol']))

    # Filter the genes dataframe to exclude those pairs
    filtered_genes_df = rna_seq_data_models_filtered[~rna_seq_data_models_filtered.apply(lambda row: (row['model_id'], row['gene_symbol']) in protein_pairs, axis=1)]


    return filtered_genes_df