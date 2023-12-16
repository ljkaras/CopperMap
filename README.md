![CopperMap Logo](CopperMap.png)

**CopperMap:** *Enhancing Success in Ullmann Couplings*

CopperMap is a Jupyter notebook developed with the goal of enhacing success in Cu-catalyzed C–N couplings (Ullmann couplings). The machine learning algorithm behind CopperMap was trained and validated using approximately 1000 experimental data points and was found to predict the yield outcome of Ullmann couplings (below or above 20% yield) between primary amines and aryl-bromides with an average accuracy of 89%. CopperMap also includes a ligand suggestion tool capable of finding the most probable working ligands for a given pair of coupling partners based on experimental results for similar compounds. 

## Manuscript Reference

CopperMap is derived from our manuscript posted on [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/6532eb5cc3693ca993c1ce40). For a more in-depth understanding of the model, its training, and validation, please refer to the manuscript.

## Features

- **Prediction Tool:** Find out the likelihood of success for a given pair of substrates.
- **Ligand Suggestion Tool:** Receive information about the ligands most likely to work for the specific coupling partners.

## Getting-Started

1. Clone this repository.
```shell
git clone https://github.com/ljkaras/CopperMap.git
```
2. Create and activate the Python environment.
```shell 
conda env create -f environments/coppermap_mac.yml  # for MacBook
```
```shell 
conda env create -f environments/coppermap_windows.yml  # for Windows
```
```shell 
conda activate coppermap
```
3. Install IPykernel and create a Jupyter kernel for the Python environment.
```shell 
conda install ipykernel
```
```shell 
python -m ipykernel install --user --name=coppermap --display-name "coppermap"
```
4. Open the notebook: `CopperMap.ipynb`.
5. Follow the instructions within the notebook for predictions and ligand suggestions.

## Citation

If you find CopperMap useful, please cite our manuscript:

> M. Samha, L. Karas, D. Vogt, E. Odogwu, J. Elward, J. Crawford, J. Steves, M. Sigman. Predicting Success in Cu-Catalyzed C–N Coupling Reactions using Data Science. ChemRxiv, 2023. DOI: 10.26434/chemrxiv-2023-f50w6.

## Credits

Developed by [The Sigman Lab](https://www.sigmanlab.com)  
Developer: Lucas Karas
