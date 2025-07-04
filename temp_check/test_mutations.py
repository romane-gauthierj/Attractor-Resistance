# # add mutations info

# # Add somatic mutations data
# # keep only the one we had for the main pipeline

# mutations_data = pd.read_csv('data/TCGA_data/TCGA_mutations_mutect2_GDC-PANCAN.csv')
# mutations_data_filtered = mutations_data[mutations_data['Sample_ID'].isin(patients_id)]
# mutations_data_filtered = mutations_data_filtered[mutations_data_filtered['gene'].isin(montagud_nodes)]

# # check if genes are TSG/ Oncogenes
# onco_tsg_gene = pd.read_csv('data/unknown_origin/oncogenes_tsg.csv')
# onco_tsg_gene = onco_tsg_gene[['Hugo Symbol', 'Is Oncogene', 'Is Tumor Suppressor Gene']]
# onco_tsg_gene_filtered = onco_tsg_gene[onco_tsg_gene['Hugo Symbol'].isin(montagud_nodes)]
# onco_tsg_gene_filtered = onco_tsg_gene_filtered.rename(columns={'Hugo Symbol': 'gene'})
# # oncogenes = onco_tsg_gene_filtered[onco_tsg_gene_filtered['Is Oncogene'] == 'Yes']
# # tsg_genes = onco_tsg_gene_filtered[onco_tsg_gene_filtered['Is Tumor Suppressor Gene'] == 'Yes']


# mutations_annotated = mutations_data_filtered.merge(
#     onco_tsg_gene_filtered[['gene', 'Is Oncogene', 'Is Tumor Suppressor Gene']],
#     on='gene',
#     how='left'
# )
# mutations_annotated = mutations_annotated.rename(columns={'Is Oncogene': 'oncogene', 'Is Tumor Suppressor Gene': 'tsg'})

# mutations_annotated = mutations_annotated[
#     mutations_annotated['oncogene'].notna() | mutations_annotated['tsg'].notna()
# ]

# # loss function mutation assumption -> TSG and
# lof_effects = ["frameshift_variant", "stop_gained", "start_lost", "splice_region_variant"]
# lof_mutations = mutations_annotated[mutations_annotated['effect'].isin(lof_effects)]
# lof_mutations_tsg = lof_mutations[(lof_mutations['tsg'] == 'Yes') & (lof_mutations['oncogene'] == 'No')]
# lof_mutations_tsg_filtered = lof_mutations_tsg[['Sample_ID', 'gene']]


# mutations_onco = mutations_annotated[(mutations_annotated['tsg'] == 'No') & (mutations_annotated['oncogene'] == 'Yes')]
# print(mutations_onco.head())
# # dna_vaf > 0.5 -> clonal mutation (mutation probably in the early tumor cells)
# # gof_effects = ['p.G12D','p.S249C', 'p.Y373C']
# # gof_mutations = mutations_onco[mutations_onco['Amino_Acid_Change'].isin(gof_effects)]
# # gof_mutations_filtered = gof_mutations[['Sample_ID', 'gene']]






# tailor_bnd_mutat_validation(lof_mutations_tsg_filtered,gof_mutations_filtered,folder_pers_models, tissue)