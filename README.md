_# Basic info

## Setting up the repository and environment locally
### Clone the repository
* `git clone https://github.com/ivanajanickova/DARTpaths-vis.git`
### Setting up the python environment
* `python3 -m venv dart-env`
* `source dart-env/bin/activate` you need to run this everytime
* `pip install -r requirements.txt` to install all necessary packages

## Running `run_enrichment_for_pathways.sh`

* you do not need to do this, since it was already done for you
* however, you can run this if you want study other pathways
* `cd` to the `Phenotype_Enrichment/backend` folder
* in the terminal window: `run_enrichment_for_pathways.sh`
* where: 
  * R-HSA-69239 is reactome pathway
  * "AHR" is the name of the pathway
  * "Phase1CompoundFunctionalization" is upper pathway name
  * `$HOME/DARTpaths-vis/` is path to the root 
  * you need to change root pathway in the bash script in case you run this on Windows machine

## Running `app.py`

* go to `Phenotype_Enrichment/visualisation/`
* run `app.py`
* click on the link that appears
