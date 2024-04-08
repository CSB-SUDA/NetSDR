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

It generated several result files, as follows:
* The signature proteins are saved in the _'signature.csv'_ file.
* The _'edges.txt'_ file stores the subtype-specific network
* The _'node_Module.txt'_ and _'edge_Module.txt'_ files provide information on the nodes and edges of robust modules, respectively.


#### Step2：Build a drug response network.





#### Step3：Perform PRS-based Drug Repurposing.





#### Step4：Run ROC analysis.






一、标题写法：  
第一种方法：  
1、在文本下面加上 等于号 = ，那么上方的文本就变成了大标题。等于号的个数无限制，但一定要大于0个哦。  
2、在文本下面加上 下划线 - ，那么上方的文本就变成了中标题，同样的 下划线个数无限制。  
3、要想输入=号，上面有文本而不让其转化为大标题，则需要在两者之间加一个空行。  
另一种方法：（推荐这种方法；注意⚠️中间需要有一个空格）
关于标题还有等级表示法，分为六个等级，显示的文本大小依次减小。不同等级之间是以井号  #  的个数来标识的。一级标题有一个 #，二级标题有两个# ，以此类推。
例如：
# 一级标题  
## 二级标题  
### 三级标题  
#### 四级标题  
##### 五级标题  
###### 六级标题 
二、编辑基本语法  
1、字体格式强调
 我们可以使用下面的方式给我们的文本添加强调的效果  
*强调*  (示例：斜体)  
 _强调_  (示例：斜体)  
**加重强调**  (示例：粗体)  
 __加重强调__ (示例：粗体)  
***特别强调*** (示例：粗斜体)  
___特别强调___  (示例：粗斜体)  
2、代码  
`<hello world>`  
3、代码块高亮  
```
@Override
protected void onDestroy() {
    EventBus.getDefault().unregister(this);
    super.onDestroy();
}
```  
4、表格 （建议在表格前空一行，否则可能影响表格无法显示）
 
 表头  | 表头  | 表头
 ---- | ----- | ------  
 单元格内容  | 单元格内容 | 单元格内容 
 单元格内容  | 单元格内容 | 单元格内容  
 
5、其他引用
图片  
![图片名称](https://www.baidu.com/img/bd_logo1.png)  
链接  
[链接名称](https://www.baidu.com/)    
6、列表 
1. 项目1  
2. 项目2  
3. 项目3  
   * 项目1 （一个*号会显示为一个黑点，注意⚠️有空格，否则直接显示为*项目1） 
   * 项目2   
 
7、换行（建议直接在前一行后面补两个空格）
直接回车不能换行，  
可以在上一行文本后面补两个空格，  
这样下一行的文本就换行了。
或者就是在两行文本直接加一个空行。
也能实现换行效果，不过这个行间距有点大。  
 
8、引用
> 第一行引用文字  
> 第二行引用文字 

