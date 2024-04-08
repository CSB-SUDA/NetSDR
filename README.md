# DRNMDRN

# Drug Repositioning for Cancers - based on network modularization and drug response network

## Graph abstract
![Framework](Picture/Framework.png)

## Highlights  
* A comprehensive computational framework integrating expression profiles, network topology, and biophysics is developed.
* The ECM module of G-IV characterized by poor prognosis can also serve as a therapeutic target.
* Combining protein expression with drug sensitivity data guides the construction of drug response networks.
* GSK1904529A displays prospect as a repurposable drug for G-IV by targeting LAMB2.

## Installation
### Requirements
* R
* Anaconda3

### Installation of DRNMDRN
Use the install_github function from the devtools library to download and install the DRNMDRN package.
```
install.packages("devtools")
library(devtools)
install_github("Bin-suda/DRNMDRN")
library(DRNMDRN)
```

### Installation of DeepPurpose and ProDy
Use anaconda3 to install the DeepPurpose tool. Open the command line terminal and enter the following command:
```
git clone https://github.com/kexinhuang12345/DeepPurpose.git
cd DeepPurpose
conda env create -f environment.yml
conda activate DeepPurpose
pip install prody
conda deactivate 
```
For more details about DeepPurpose, please visit https://github.com/kexinhuang12345/DeepPurpose

## Usage

### Package organization
```
DRNMDRN/
├── LICENSE
├── README.md
├── R/           <- Contains R scripts for the above framework calculations.
├── data/        <- Contains data files used by the package.
├── man/         <- Contains .Rd help documentation.
├── DESCRIPTION  <- Package metadata.
└── NAMESPACE    <- Package namespace information.
```
[pipline]

### Usage Example

#### Step1：Identify subtype-specific modules.
```
# Load the proteome expression profiles and grouping information of subtypes.
data(Expr_Group)
getModule(exprDF=Expr,groupDF=Group,subtype="G-IV")
```
The `getModule` function accepts three parameters:
* _exprDF_: A data frame storing expression values, with rows representing proteins and columns representing samples.
* _groupDF_: A data frame storing subtype information, with the first column being samples and the second column being subtype grouping.
* _subtype_: A vector representing the analysis of a specific subtype.

It generates several result files, as follows:
* The signature proteins are saved in the _'signature.csv'_ file.
* The _'edges.txt'_ file stores the subtype-specific network
* The _'node_Module.txt'_ and _'edge_Module.txt'_ files provide information on the nodes and edges of robust modules, respectively.


#### Step2：Build a drug response network.
```
# for example
moduleDF <- Expr[1:50,1:20]
data(Expr_Group)
sig_count <- getModuleResponse(moduleDF)
drugDF <- read.csv("calcPhenotype_Output/DrugPredictions.csv",row.names = 1,check.names = F)
DRN_ls <- getDRN(moduleDF,drugDF)
```
At first, it predicts the clinical drug response of patients within a specific subtype based on module expression `moduleDF` specific to that subtype. Then, They are utilized to predict potential interactions between the drug and the protein.
It returns a list containing the drug response network and its node information.

#### Step3：Perform PRS-based Drug Repurposing.





#### Step4：Run ROC analysis.



