# Use attractors to explore drug resistance

## Set up

1. download the data and saved it in the data folder 
2. create virtual environment and install requirements txt

```
pip install -r requirements.txt
```

3. Pipeline_validation_prostate: 
create personalized boolean networks with TCGA prostate data, and compare to gleason score (Montagud paper)

Pipeline_generic:
Run MaBoSS simulation on generic model (no personalization) 

Pipeline_pers:
Create personalized Boolean Netorks with Cell Model Passport data of PAN CANCER (Except Hematopoietic tissue)