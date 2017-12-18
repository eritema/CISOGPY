## Usage
Please run `python autoCIS.py -h` for usage
i.e.
```
python autoCIS.py -i fileIn -o fileOut.csv -H
```

## Input format

In order to proceed with the analysis with autoCIS, the data must be in the following format (tsv,csv)

```
---------------------------------------------------------------
 Chr          Pos  Gene    Ann1  Ann2                  Entropy 
                                                       Label   
---------------------------------------------------------------
 chr1     1753243  GNB1          NM_001282538_Intron3  P1      
 chr12  125578270  AACS          NM_023928_Intron5     P1      
 chr12   32860361  DNM1L         NM_001278465_Intron4  P1      
 chr12    3931561  PARP11        NM_020367_Intron4     P1      
 chr15   93526100  CHD2          NM_001271_Intron24    P1      
 chr16   50811009  CYLD          NM_015247_Intron8     P1      
---------------------------------------------------------------
```

## Conda environement (osx)
To build the appropriate conda env for autoCIS run the conda command
```
conda create --name autoCIS --file spec-cisog.txt
source activate autoCIS
```
