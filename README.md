# **MethAgingDB: a comprehensive DNA methylation database for aging biology**

![MINGLE](./MethAgingDB.png)

## Web interface

https://methagingdb.biox-nku.cn/

## Script description

| Script               | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| getDMS.R             | Get differentially methylation sites (DMSs) between age groups. |
| getDMR_Bumphunter.R  | Get differentially methylation regions (DMRs) between age groups using Bumphunter method. |
| getDMS_continuous.R  | Get differentially methylation sites (DMSs) using age as continuous variable. |
| getDMR_DMRcate.R     | Get differentially methylation regions (DMRs) between age groups using DMRcate method. |
| run_tSNE.py          | Get tSNE scatter plots of each dataset.                      |
| run_PCA.py           | Get scatter plots of each dataset using PCA.                 |
| run_tSNE_removeXY.py | Get tSNE scatter plots of each dataset without sex chromosomes. |
| run_PCA_removeXY.py  | Get scatter plots of each dataset using PCA without sex chromosomes. |
