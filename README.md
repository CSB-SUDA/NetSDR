# NetSDR

# Drug repurposing for cancers based on subtype-specific network modularization and drug perturbation analysis

## Framework
![Framework](Picture/Framework.png)

## Highlights  
* We proposed the NetSDR framework for drug repurposing for different cancer subtypes.
*	NetSDR integrated, proteomics and drug sensitivity data by network topology and biophysics-based analysis.
*	Apply NetSDR to gastric cancer revealed GSK1904529A as a repurposable drug for G-IV by targeting LAMB2.

## Installation
### Requirements
* R
* Anaconda3

This package runs in R, but its functionality depends on Python scripts. Therefore, before using it, you need to create a virtual environment named `py3.9` using anaconda3, with the Python version no lower than 3.7. Please open the Anaconda Prompt and enter the following command:
### Install Python
```
conda create -n py3.9 python=3.9
conda activate py3.9
pip install rdkit
pip install descriptastorus 
pip install DeepPurpose
pip install seaborn
pip install goatools
pip install prody
conda deactivate
```
The Python environment is created to use the DeepPurpose tool and PRS analysis. For more details about DeepPurpose, please visit https://github.com/kexinhuang12345/DeepPurpose.

> **Note that the path to the py3.9 environment should be obtained using `conda env list`. In subsequent analyses, this path will be used as a parameter in the functions of NetSDR.**


### Install NetSDR package in R
Use the install_github function from the devtools library to download and install the NetSDR package.
```
install.packages("devtools")
library(devtools)
install_github("CSB-SUDA/NetSDR")
library(NetSDR)
```

## Usage

### Package organization
```
NetSDR/
├── LICENSE
├── README.md
├── R/           <- Contains R scripts for the above framework calculations.
├── data/        <- Contains data files used by the package.
├── inst/        <- Contains example data and Python scripts.
│     ├── python/           <- Contains some code from enm_package.For more details, please visit https://github.com/oacar/enm_package.
│            └── enm/       <- Source code for enm.
│     ├── extdata/          <- Some example data.
├── man/         <- Contains .Rd help documentation.
├── DESCRIPTION  <- Package metadata.
├── NAMESPACE    <- Package namespace information.
└── NetSDR.Rproj <- Rproject.
```
If installing NetSDR via `devtools` fails, its source code can be directly downloaded. Then, open the `NetSDR.Rproj` file in R and run `devtools::build()` to create a local installation package for local installation.

### Usage Example

#### Step1: Identify subtype-specific modules.
The `getModule` function can identify subtype-specific modules, which accepts three parameters:
* _expr_file_: The path of TXT file storing expression profile. The specific format is referred to in the 'inst/extdata/expression. txt' file or as follow:
  
  || sample1 | sample2 |... |
  | --- | --- | --- | --- |
  | protein1 | 2.345 | 6.480 | ... |
  | protein2 | 7.985 | 4.621  | ... |
  | ... | ... | ...  |...|
  
* _group_file_:  The path of TXT file storing subtype information. In this file, the first column is samples and the second column is subtype grouping.
  
  |Sample| Group |
  | --- | --- |
  | sample1 | subtype1 |
  | sample2 | subtype1 |
  | sample3 | subtype2 |
  | ... | ... |
  
* _subtype_: A vector representing the analysis of a specific subtype. For example, subtype="subtype1".

* _ppi_file_: The path of TXT file storing interactome data.

  |node1| node2 |
  | --- | --- |
  | protein1 | protein2 |
  | protein1 | protein3 |
  | protein2 | protein4 |
  | ... | ... |

```
expr_file <- "inst/extdata/expression.txt"
group_file <- "inst/extdata/group.txt"
subtype <- "G-IV"
ppi_file <- "inst/extdata/ppi_String400.txt"
getModule(expr_file,group_file,subtype,ppi_file)
```

It generates `NetSDR_results/Modules` result file, including:
* The signature proteins are saved in the _'signature.csv'_ file.

  |protein| log2FoldChange | p.value | frequency | label | 
  | --- | --- | --- | --- | --- |
  | pro1 | -1.41 | 0.00046 | 20 | down |
  | pro2 | 4.26 | 0.003 | 18 | up |
  | ... | ... | ... | ... | ... |
* The _'network.txt'_ stores a network composed of signature proteins.
  
  |node1| node2 |
  | --- | --- |
  | pro1 | pro2 |
  | pro1 | pro3 |
  | ... | ... |
  
* The _'edges.txt'_ file stores the subtype-specific network.

  |node1| node2 | module |
  | --- | --- | --- |
  | pro1 | pro2 | M1 |
  | pro1 | pro3 | M2 |
  | ... | ... | ... |
  
* The _'node_Module.txt'_ and _'edge_Module.txt'_ files provide information on the nodes and edges of robust modules, respectively.

  |node| module |
  | --- | --- |
  | pro1 | M1 |
  | pro2 | M1 |
  | ... | ... |
  
* The _'node_Module_select.txt_' provides information on modules with more than 9 nodes.

  |node| module |
  | --- | --- |
  | pro1 | M1 |
  | pro2 | M1 |
  | ... | ... |

Then, the module associated with drug response can be identified by `getModuleResponse` function . A drug response network is constructed on this module, and drug repositioning based on PRS is performed.

#### Step2: Build a drug response network.
The `getDRN` function is used to build a drug response network, with _"expr_file"_ as input, namely:
* _expr_file_: The path of TXT file storing expression profile of a module, It contains the following information：

  |  | sample1 | sample2 | ... |
  | --- | --- | --- | --- |
  | protein1 | 7.335 | 1.345 | ... |
  | protein2 | 6.782 | 8.481  | ... |
  | ... | ... | ...  |...|
  
```
# For example
expr_file <- "inst/extdata/expression_module.txt"
getDRN(expr_file)
```
It exports `NetSDR_results/DRN` result file, including:
* _calcPhenotype_Output_: The results predicted by DeepPurpose tools. The predicted clinical drug reactions of patients are saved in the "DrugPredictions.csv" file.

  |  | drug1 | drug2 | ... |
  | --- | --- | --- | --- |
  | sample1 | 12.451 | 0.567 | ... |
  | sample2 | 8.392 | 6.179  | ... |
  | ... | ... | ...  |...|
  
* _DRN.txt_: The predicted drug response network.

  |protein| drug | pvalue | label | 
  | --- | --- | --- | --- |
  | pro1 | drug1 | 0.0001 | sig |
  | pro2 | drug1 | 0.0023 | 18 |
  | ... | ... | ... | ... |
  
* _DRN_info.txt_: Annotation information of nodes in the drug response network.

  |node| type |
  | --- | --- |
  | pro1 | protein |
  | drug1 | drug |
  | ... | ... |

#### Step3: Perform PRS-based Drug Repurposing.
It initially predicts the binding affinity between drugs and proteins in the DRN, and then performs a drug repositioning method based on PRS to rank the DPIs.
The `getAffinity` function can predicts the binding affinity, with five parameters as input:
* _smiles_file_: The path of TXT file storing drug smiles. The specific content is:

  | drug | SMILES |
  | --- | --- |
  | drug1 |	C1=CC=C(C=C1)NC(=O)CCCCCCC(=O)NO |
  | drug2 |	C1=CC(=CC=C1/C=C\2/C(=O)NC(=N)S2)O |
  | drug3 | C1=C(C(=O)NC(=O)N1)F |
  | ... | ... |

* _seq_file_: The path of TXT file storing protein sequence. The specific content is:

  | protein | Sequence |
  | --- | --- |
  | protein1 | MAADISESSGADCKGDPRNSAKLDADYPLRVLYCGVCSLPTEYCEYMPDVAKCRQWLEKNFPNEFAKLTVENSPKQEAGISEGQGTAGEEEEKKKQKRGGRGQIKQKKKTVPQKVTIAKIPRAKKKYVTRVCGLATFEIDLKEAQRFFAQKFSCGASVTGEDEIIIQGDFTDDIIDVIQEKWPEVDDDSIEDLGEVKK |
  | protein2 | MPLAKDLLHPSPEEEKRKHKKKRLVQSPNSYFMDVKCPGCYKITTVFSHAQTVVLCVGCSTVLCQPTGGKARLTEGCSFRRKQH |
  | protein3 | MPGPTPSGTNVGSSGRSPSKAVAARAAGSTVRQRKNASCGTRSAGRTTSAGTGGMWRFYTEDSPGLKVGPVPVLVMSLLFIASVFMLHIWGKYTRS |
  | ... | ... |
  
* _DRN_file_: The path of the TXT file storing the predicted drug response network by `getDRN` function, such as: DRN_file = "NetSDR_results/DRN/DRN.txt".
* _py_env_: The virtual environment path for py3.9 installed using anaconda3, which can be obtained by running `conda env list` in `Anaconda Prompt`.
* _pretrained_model_: The DeepPurpose model type used, which defaults to 'MPNN_CNN_BindingDB_IC50'.

```
# Predict the binding affinity between drugs and proteins
smiles_file <- "inst/extdata/smiles.txt"
seq_file <- "inst/extdata/sequences.txt"
DRN_file <- "NetSDR_results/DRN/DRN.txt"  #"inst/extdata/DRN.txt"
py_env <- "D:/Software/anaconda3/envs/py3.9"
getAffinity(smiles_file,seq_file,DRN_file,py_env)
```
It outputs the result file "result/virtual_screening.txt" in "NetSDR_results/PS" dir, which is used for ps calculation.The content of "virtual_screening.txt" is as follows:

  +------+-------------------------------+-------------+---------------+
  
  | Rank |           Drug Name           | Target Name | Binding Score |
  
  +------+-------------------------------+-------------+---------------+
  
  |  1   |             drug4             |    pro11    |     61.69     |
  
  |  2   |             drug1             |    pro2     |     38.48     |
  
  |  3   |             drug10            |    pro1     |     37.10     |
  

  
The `getPS` function accepts four parameters:
* _BA_file_: The path to the binding affinity file (.txt) predicted by the `getAffinity` function, such as "NetSDR_results/PS/result/virtual_screening.txt", whose format is shown above.
* _PPIN_file_: The path to the PPI network file (.txt) . This file contains protein-protein interactions with columns named gene1 and gene2.

  |gene1| gene2 |
  | --- | --- |
  | pro1 | pro2 |
  | pro1 | pro3 |
  | ... | ... |
  
* _py_env_: The virtual environment path for py3.9 installed using anaconda3, same as  `getAffinity`.

```
# Calculat drug score(ps) using PRS methods for DRN.
py_env <- "D:/Software/anaconda3/envs/py3.9"
BA_file <- "NetSDR_results/PS/result/virtual_screening.txt"  #"inst/extdata/virtual_screening.txt"
PPIN_file <- "inst/extdata/edge_module.txt"
getPS(BA_file,PPIN_file,py_env)
```
It returns a data frame contains the ps and saves it in the "NetSDR_results/PS/prs_dti_score.csv" file. The format is as follows:

  | Drug.Name | Target.Name | Binding.Score | sens | ps |
  | --- | --- | --- | --- | --- |
  | drug1 | pro1 | 21.9 | 5.42 | 118.89 |
  | drug1 | pro2 | 16.21 | 3.21 | 52.03 |
  | ... | ... | ... | ... | ... |

***Run ROC analysis***
```
# Read the file containing the ps scores obtained in the previous step
DPI_ps_file <- "inst/extdata/dpi_ps_score.txt"
runROC(DPI_ps_file)
```
The parameter _DPI_ps_file_ of function `runROC` is the path of the file containing DPIs and their ps. It perform the ROC analysis for different ps cut-off in sequence and returns their AUCs.

> ***Note:*** _For more detailed information about each function, please refer to the function comments in the respective script._
