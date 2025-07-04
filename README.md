# Use attractors to explore drug resistance

## Set up

1. download the data and saved it in the data folder 
2. create virtual environment and install requirements txt

```
pip install -r requirements.txt
```

3. Pipeline_validation_prostate: 
create personalized boolean networks with TCGA prostate data, and compare to gleason score (Montagud paper)




pipeline generic: run maboss simulation on generic model (no personalization) 
Pipeline 0: Validate the personalization with the TCGA dataset and prostate model of Montagud Analysis(approx. 20 min)
pipeline 1: create personalized networks (gene expression, cnv) (approx. 4 min)
pipeline 2: run maboss simulation and gene diff expressed analysis
pipeline 3: simulate a KO of genes identified and check what lead to reduction of proliferation (approx 3 hours)













Run the Pipeline pers to create the personalized models (by modifying the models and results folder with 'interventions')-> until step 6 (included)
Run downstream analysis gene intervention 


check how to create the prolif_egf enrichment file and then analyze it (should i have three pipelines ?- one for creating the models, one for computing the gene enrichment and one for identify genes to block?)



