# Devono-Mississippian winds over the greater US Midcontinent 

## purpose
Analysis of the simulated wind fields and precipitation in the Devono-Mississippian simulations of [Valdes et al. (2021)](https://cp.copernicus.org/articles/17/1483/2021/). Particular focus is on a possible eolian dust source for strata of the North American midcontinent. Paleolocations for the study sites and potential dust sources are reconstructed with [pyGPlates](https://www.gplates.org/docs/pygplates/) using my [template repo](https://github.com/sebsteinig/pyGPlates-reconstructions-template).

## source data
HadCM3BL model output from the BRIDGE webpage processing. Simulations are from the `Scotese_02` set as described in [Valdes et al. (2021)](https://cp.copernicus.org/articles/17/1483/2021/). Data can be accessed at the [BRIDGE website](https://www.paleo.bristol.ac.uk/ummodel/scripts/papers/Valdes_et_al_2021.html). 

## running the notebooks
Notebooks can either be run on [Google Colab](https://colab.research.google.com/) (online, Google account required) or locally. Notebooks should include a button to open it in Colab, otherwise you can also directly load a GitHub repo within Colab. Easiest way to run locally is to first download the repo with


```
git clone https://github.com/USERNAME/REPOSITORY
``` 

and then install [conda](https://conda.io/projects/conda/en/latest/index.html) (if not installed already). Then create an environment `env_name` with 

```
conda env create --name env_name --file=environment.yml
``` 

using the `environment.yml` file from this repository to install all necessary python packages. The notebooks can then be run interactively by typing

```
jupyter lab
```

We further use the [jupytext](https://jupytext.readthedocs.io/en/latest/index.html) package to automatically save a pure python script (*.py) each time the notebook is saved. This is very helpful for clean diffs in the version control and allows you to run the analysis in your local terminal with:

```
python notebook_name.py
```
The python file can also be shared with others to work on the code together using all the version control benefits (branches, pull requests, ...). You can edit it with any tex editor/IDE and it can also be converted back to a jupyter notebook (with no output) via
```
jupytext --to notebook notebook_name.py
```

