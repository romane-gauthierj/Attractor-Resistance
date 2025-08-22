# A Pan-Cancer Computational Pipeline for Modeling Adaptive Resistance to PI3K Inhibitor Pictilisib through Markovian Boolean Network Simulations



## Project Overview

This project provides a computational framework to analyze and model cell behavior, with a focus on understanding gene expression and resistance mechanisms. The pipeline integrates Boolean network modeling and stochastic simulations to identify key drivers of resistance and validate findings across datasets. It includes:

- A main analysis pipeline for discovering resistance pathways and candidate genes.
- A validation pipeline to ensure robustness using independent data.
- A baseline pipeline for initial model setup and comparison.

This approach enables systematic exploration of resistance mechanisms and supports the development of personalized cell models.


## Dependency
1. Install miniconda (https://www.anaconda.com/docs/getting-started/miniconda/install)

2. Install MaBoSS dependency
```
conda install -c colomoto::maboss
``` 


## Set up

1. create directory

```
mkdir attractors_project
``` 

2. go inside the directory
```
cd attractors_project
``` 

3. clone the github repo (with all the code) inside the directory 

```
git clone https://github.com/romane-gauthierj/Attractor-Resistance.git
``` 

4. create virtual environment 

```
python -m venv .env
``` 

5. activate virtual environment 
```
source .env/bin/activate
``` 


6. Install requirements txt file that contains all the required library to run the code. In the terminal run:

```
pip install -r requirements.txt
```


7. Load the datasets (original datasets are saved in a drive folder)

```
gdown --folder https://drive.google.com/drive/folders/1Tp_wRLTVEFLWm_mwrubFYBSsbt6xTnQG

```

8. Run the pipelines with the default parameters 

Validation pipeline with the command (100 cell models, personalization using mutations and genes):

```
python run_validation_pipeline.py
```

Main pipeline with the command (35 cell models in each group, Pictilisib drug, personalization using mutations and genes, sigmoid normalization):

```
python run_main_pipeline.py
```


## Output - Results
The results and output figures can be found in the analysis folder in the respective folder. 






## Configuration- User parameters

It is also possible to change some default parameters.  

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













