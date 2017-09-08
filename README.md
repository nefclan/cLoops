# Clustering based loops calling for ChIA-PET, HiChIP and high resolution Hi-C data with cLoops

# not available yet 

## Introduction
By taking the mapped paired-end tags from ChIA-PET or HiChIP as 2D points, the problem for calling loops is converted to draw significant clusters from sparse points with noise. After classifying the detected clusters into self-ligation and inter-ligation clusters, the significances of the inter-ligation clusters are estimated using permuted local backgrounds. We implemented the approach in the “cLoops (see loops)” package. Although without the peak calling step, the anchors determined by cLoops shows a high overlap with the peaks. By comparing to peaks based loop calling tools, we show that cLoops can detect more interactions with better ranked p-values, better supported by Hi-C data, sharper anchors, higher enrichment for TF motifs, work well both for sharp and broad peak like ChIA-PET data.

If you find cLoops is useful, please cite our paper:    
**### Clustering based loops calling for ChIA-PET, HiChIP and high resolution Hi-C data with cLoops ###**

--------
## Install
[scipy](https://www.scipy.org/),[numpy](http://www.numpy.org/), [seaborn](https://seaborn.pydata.org/), [pandas](http://pandas.pydata.org/) and [joblib](https://pythonhosted.org/joblib/) are required. If you have problems for installing scipy, please refer to [Anaconda](https://docs.continuum.io/anaconda/) or [SAGE](http://www.sagemath.org/).
```
wget 
tar xvzf cLoops.tar.gz
cd cLoops
python setup.py install    
```

or just

```
pip install cLoops
```
Please refer to [here](https://docs.python.org/2/install/index.html) to install cLoops to customized path.

--------
## Usage
Run ***cLoops -h*** to see all options. Key parameters are ***eps*** and ***minPts***  ***minPts*** defines at least how many PETs are required for a candidate loop, ***eps*** defines the distance requried for two PETs being neighbors. For ChIA-PET data with sharp peaks, cLoops can auto estimate ***eps*** from the data as 2 fold of the fragment size, and ***minPts***=5 is good. For ChIA-PET data with broad peaks (like H3K4me1), empirical experience is set ***eps*** to 2000. For HiChIP, set a series ***eps***=2000,4000,6000,8000,10000 & ***minPts***=20 worth a first trial, if sequencing deep, increase ***minPts*** to 30 or 50. For practically usage, using the PETs in the smallest chromosome except chrY and chrM, then run a series of ***eps***, choose the smallest ***eps*** that can get well seperated inter-ligation and self-ligation PETs distance distributions. 

--------
### Input  
Mapped PETs in [BEDPE format](http://bedtools.readthedocs.io/en/latest/content/general-usage.html), compressed files with gzip are also accepected, first 6 columns as following are necessary: chrom1,start1,end1,chrom2,start2,end2.

--------
### Output
The main output is a loop file and a PDF file or PDFs for the plot of self-ligation and inter-ligation PETs distance distributions.
For the .loop file, columns and explaination are as follwing:

column | name | explaination
------ | ---- | ------------
0th | loopId | Id for a loop, like chr1-chr1-1
1th | ES | Enrichment score for the loop, caculated by observed PETs number divided by the mean PETs number of nearby permutated  regions
2th | FDR | false discovery rate for the loop, caculated as the number of permutated regions that there are more observed PETs than the region  
3th | binomal\_p-value | binomal test p-value for the loop
4th | distance | distance (bp) between the centers of the anchors for the loop
5th | hypergeometric\_local\_FDR | FDR for the hypergeometric test p-value compared to permutated regions
6th | hypergeometric\_p-value | hypergeometric test p-value for the loop
7th | iva | genomic coordinates for the left anchor, for example, chr13:50943050-50973634
8th | ivb | genomic coordinates for the right anchor
9th | poisson_p-value | poisson test p-value for the loop
10th | ra | observed PETs number for the left anchor
11th | rab | observed PETs number linking the left and right anchors
12th | rb | observed PETs number for the right anchor
13th | poisson\_p-value\_corrected | Bonferroni corrected poisson p-value according to number of loops for each chromosome
14th | binomal\_p-value\_corrected | Bonferroni corrected binomal p-value according to number of loops for each chromosome
15th | hypergeometric\_p-value\_corrected | Bonferroni corrected hypergeometric p-value according to number of loops for each chromosome
16th | significant | 1 or 0, 1 means we think the loop is significant compared to permutated regions. For ChIA-PET data, significant requiring ES >=1.0, FDR <=0.05, hypergeometric\_local\_FDR <=0.05 and all uncorrected p-values <= 1e-5; For HiChIP and high resolution Hi-C data, significant requiring ES >= 2.0, FDR <=0.05, sqrt(binomal\_p-value * poisson\_p-value) <=1e-5. You can ignore this and customize your cutoffs by visualization such like in the [Juicebox](https://github.com/theaidenlab/juicebox) to determine your cutoffs.

--------
## Examples
All following examples source data, result and log file can be found in the [examples](https://github.com/YaqiangCao/cLoops/tree/master/examples).

1. ChIA-PET data    
We provide a test data from GM12878 CTCF ChIA-PET ([GSM1872886](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1872886)), just the chromosome 21 mapped to hg38. Run the command as following then you will get the result if cLoops is successfuly installed. The ***eps*** is auto estimated and default ***minPts*** is 5,**-w** option will generate loops for visualization in [washU browser](http://epigenomegateway.wustl.edu/browser/),**-j** option will generate loops for visualization in [Juicebox](https://github.com/theaidenlab/juicebox) .
```
wget https://github.com/YaqiangCao/cLoops/blob/master/examples/GSM1872886_GM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz
cLoops -f GSM1872886_GM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz -o chiapet -w 1 -j 1
```      
For ChIA-PET data with sharp peak, like the CTCF here, you will get the inter-ligation and self-ligation PETs distance distribution like [this](https://github.com/YaqiangCao/cLoops/blob/master/examples/chiapet_disCutoff.pdf). If your experimental data doesn't look like this by auto estimated ***eps***, which could be true for some ChIA-PET data with broad peak (like H3K27ac), please use the small chromosome (chr21 in human and chr19 in mouse) run a series of ***eps***, then chose the smallest one that generate the well seperated distance distribution to run cLoops, or just using the series. 

2. HiChIP data   
We provide two data of from GM12878 cohesin HiChIP of two biological replicates, just the chromosome 21 mapped to hg38. Run the command as following to call merged loops. ***-s*** option is used to keep working directory and temp files, which could be used by deLoops,jd2washU (BEDTOOLS needed) and jd2juice (Juicer needed).***-hic*** option means using cutoffs design for Hi-C like data, see above. 
```
wget https://github.com/YaqiangCao/cLoops/blob/master/examples/GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio1.bedpe.gz 
wget https://github.com/YaqiangCao/cLoops/blob/master/examples/GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio2.bedpe.gz 
cLoops -f GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio1.bedpe.gz,GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio2.bedpe.gz -o hichip -eps 1000,2000,4000,6000,8000,10000 -minPts 50 -s 1 -hic 1 -w 1 -j 1
```    

--------
## Questions & Answers  
Please address questions and bugs to Yaqiang Cao (caoyaqiang@picb.ac.cn) or Xingwei Chen (chenxingwei@picb.ac.cn) or Daosheng Ai (aidaosheng@picb.ac.cn), using the subject as "cLoops: questions about" to escape misjudged as spams.  
