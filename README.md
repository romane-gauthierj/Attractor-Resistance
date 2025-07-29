# Use attractors to explore drug resistance

## Set up

1. download the data and saved it in the data folder (TCGA datasets, Cell model passport)
2. create virtual environment and install requirements txt

```
pip install -r requirements.txt
```


## User parameters 
1. select the drug, the drug target and the type of models (proteins and genes)



## Pipelines

1. pipeline generic: Growth Truth of the personalization, run maboss simulation on generic model (no personalization) 
2. Pipeline 0: Validate the personalization with the TCGA dataset: create personalized boolean networks with TCGA prostate data, and compare to gleason score (Montagud paper)(approx. 20 min)
3. pipeline 1: create personalized networks (gene expression, cnv) (approx. 4 min) for resistant, sensitive and healthy patients. 
4. pipeline 2: run maboss simulation and gene diff expressed analysis
5. pipeline 3: simulate a KO of genes identified and check what lead to reduction of proliferation (approx 3 hours for all the genes)






step 1
- import data montagud cfg et bnd et enlever tous les fused event et ar_erg dans le modele generic 
- run pipeline 1










Run the Pipeline pers to create the personalized models (by modifying the models and results folder with 'interventions')-> until step 6 (included)
Run downstream analysis gene intervention 


check how to create the prolif_egf enrichment file and then analyze it (should i have three pipelines ?- one for creating the models, one for computing the gene enrichment and one for identify genes to block?)



