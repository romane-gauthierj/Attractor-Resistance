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

Pipeline generic: run MaBoSS on generic boolean model (Montagud)
Main Pipeline:  create personalized cell models boolean networks from the Montagud boolean networks based on expression data, run the MaBoSS stochastic simulation, compute correlation of probability phenotype and genes signature, perform genes differential expression and identify potential genes drivers of resistance.
Simulate KO of some potential candidates and compute invasion mean value as well as stats test between before KO and after

Pipeline Validation: Create personalized Boolean networks from Fumia Boolean network with METABRIC datasets, run MaBoSS, compute correlation between phenotype probability and Proliferation/ Apoptosis as well as to MIP 




step 1
- import data montagud cfg et bnd et enlever tous les fused event et ar_erg dans le modele generic 
- run pipeline 1










Run the Pipeline pers to create the personalized models (by modifying the models and results folder with 'interventions')-> until step 6 (included)
Run downstream analysis gene intervention 


check how to create the prolif_egf enrichment file and then analyze it (should i have three pipelines ?- one for creating the models, one for computing the gene enrichment and one for identify genes to block?)



