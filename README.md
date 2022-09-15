# contig-extract
Python3  
Extract contigs from fasta files  

## Common Usage  
1. extract/remove single contig from one fasta file  
```
python extractor2.py -i example.fasta -e/-r 1
```
2. run blast and extract contigs containing target genes.  
```
python extractor2.py -i example.fasta -b PlasmidFinder.fasta
```
3. use [mlplasmids](https://gitlab.com/sirarredondo/mlplasmids) output to divide genome into chromosome and plasmids.  
```
python extractor2.py -i example.fasta -p example.tab
```

## Batch running  
1. use "-i" and run it in a loop.
2. use "-e" or "-r" with a list file (WITHOUT "-i") 
to extract/remove contigs from multiple fasta files
```
python extractor2.py -e list.txt
```  
> List file contains fasta file names and contig names to extract or remove, separated by tab, with no heading.  

  Column 1  | Column 2
  ----  | ----
  file name1  |  contig name1  
  file name1  |  contig name2  
  file name1  |  contig name3  
  file name2  |  contig name1  
  ......  |  
  
### Versions
#### V2.0（in master branch)
The output contigs were put into one file. 1	input -> 1 output.  

#### V2.1 (in THIS branch)
The output files were separated by contigs. 1 input -> n output.  
