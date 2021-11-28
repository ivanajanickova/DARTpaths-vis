_# Basic info

## Setting up the repository and environment locally
### Clone repo
* `git clone https://github.com/ivanajanickova/DARTpaths-vis.git`
### Set-up the python evironment
* `python3 -m venv dart-env`
* `source dart-env/bin/activate`   #run every time
* `pip install -r requirements.txt`

## Running `Python_APRIL_2021_phenotype_enrichment.py`
* in /DARTpaths-vis/databases folder create folders: phenotype, orthologs, ontology
* `cd` to the `Python_Phenotype_Enrichment-version_April_2021` folder
* in the terminal window: `python3 Python_APRIL_2021_phenotype_enrichment.py R-HSA-8937144 "AHR" $HOME/DARTpaths-vis/`
* where: 
  * R-HSA-69239 = reactome pathway
  * "AHR" = pathay name
  * $HOME/DARTpaths-vis/ = path to the root (WILL BE DIFFERENT IN WINDOWS)

## Extending code
* activate env
* `git pull` (in the `master` branch the last changes)
*  always commit to a branch different from the `master` 
*  create a pull request and assign reviews 
