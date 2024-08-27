# mappability-diff
Compare mappability tracks 



## Setup and Installation 

```
conda create -n mappability-diff

conda install bioconda::ucsc-wigtobigwig python=3.10 

python -m pip install -r requirements.txt
``` 


## Run 

Create `data` directory with `fasta` and `bigwig` subdirectories. 
```
mkdir -p data/fasta
mkdir -p data/bigwig 
```
Move your fasta files to `data/fasta` before proceding


Run genmap to create mappability tracks 
```
bash run_genmap.sh
```

Then assuming you have downloaded an annotation file for your genome run: 
```
python src/mappability_diff.py data/bigwig results --gtf annotation.gtf
```

