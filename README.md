# A Pan-Cancer Computational Pipeline for Modeling Adaptive Resistance to PI3K Inhibitor Pictilisib through Markovian Boolean Network Simulations


## Project Pipelines
This project uses a three-step pipeline to analyze and model cell behavior, focusing on gene expression and resistance pathways.

- ### Main Analysis Pipeline
  This is the core pipeline for understanding cell resistance. It focuses on identifying key genes and pathways that drive resistance.

  * #### Personalize Cell Models
    We create custom Boolean network models for individual cells using the Montagud network as a base and incorporating specific gene expression data.

  * #### Stochastic Simulation
    The MaBoSS stochastic simulation is run on these personalized models to predict cell behavior.
  * #### Genes differential analysis and simulation of gene KO
    Identify Drivers of Resistance: We analyze the correlation between phenotype probabilities and gene signatures. We also perform gene differential expression analysis to pinpoint potential "driver genes" of resistance.

    Validate Candidates: We simulate the knockout (KO) of these potential candidate genes and analyze how their removal affects invasion rates, using statistical tests to confirm the impact.

- ### Validation Pipeline
  This pipeline validates the core model using a different dataset to ensure its accuracy and robustness.

  * #### Personalize Cell Models
    Personalize from METABRIC: We create personalized Boolean networks using the Fumia network and METABRIC datasets.

  * #### Stochastic Simulation
    Simulate with MaBoSS: We run MaBoSS simulations on these networks.

  * #### Compute NIP correlation, and pathway genes signature correlation

    Compute Correlation: We calculate the correlation between the phenotype probabilities and key biological processes like Proliferation and Apoptosis, as well as the MIP (Most Important Pathogen).



- ### Baseline(Generic) Pipeline
This is the baseline pipeline used for initial model setup and validation.

  * #### Stochastic Simulation
    Baseline Simulation: We run a MaBoSS simulation on the generic Montagud Boolean model before any personalization is performed. This provides a baseline for comparison with the personalized models.





## Set up

1. download the data and saved it in the folder data and in the subfolder corresponding(TCGA datasets, Cell model passport)
2. create virtual environment 


3. activate it, in the terminal type the following command. If properly activated, you should see (.env) on left to the terminal window.

```
.env/bin/activate
```


3. Install requirements txt file that contains all the required library to run the code. In the terminal run:

```
pip install -r requirements.txt
```


## Configuration- User parameters

This pipeline is designed so that users can easily change variables to fit their needs.  

You can configure:
- **Number of cell models/patients**
- **Drug of interest**
- **Normalization technique**
- **Other parameters**

All variables are stored in the [`config.env`](./config.env) file.  
To update a value, edit the line corresponding to the variable.  
You can also change dataset paths in this file.  

Note: All variables for the **three pipelines** are in the same `config.env`.  
Be careful to edit the variables corresponding to the pipeline you are running.

---

### Available options

```env
genetic_intervention = KO | KI
drug_name            = Pictilisib | AZD7762
continuous_variable  = genes | proteins | genes_proteins
discrete_variable    = mutations | cnv
normalization_techniques = sigmoid | min-max | ...
```



## Run the pipelines

1. Validation pipeline with the command:

```
python run_validation_pipeline.py
```

2. Main pipeline with the command:

```
python run_main_pipeline.py
```


## Output - Results
The results and output figures can be found in the analysis folder. 



