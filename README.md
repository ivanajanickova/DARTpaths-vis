## Setting up the repository and environment locally
### Clone the repository
* `git clone https://github.com/ivanajanickova/DARTpaths-vis.git`
* preferably, start this project in PyCharm or other python IDE
*
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
  * "Phase1CompoundFunctionalization" is upper level pathway name
  * `$HOME/DARTpaths-vis/` is path to the root 
  * you need to change root pathway in the bash script in case you run this on Windows machine

## Running `app.py`

* go to `Phenotype_Enrichment/visualisation/`
* run `app.py`
* on an average laptop, this step will take ~ 15-20 minutes
* click on the link that appears


## Troubleshooting
* Windows users can experience issues with the `run_enrichment_for_pathways.sh` since the `Python_APRIL_2021_phenotype_enrichment.py` is not optimised for Windows machines
* We recommend running this step on Linux machine
* You may have to install some of the packages from requirements.txt manually

