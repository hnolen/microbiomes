#QIIME2 Pipeline for microbiome analysis of soil samples


version qiime2-2020.2


### 1. Import demultiplexed reads, trim primer sequences using cutadapt, and create base call quality plots
view plots at https://view.qiime2.org/
`sbatch import_trim.slurm`


### 2. Error correction using DADA2
`sbatch qiime-dada.slurm` 

### 3. Visualize DADA2 outputs - look at dns.qzv to determine rarefaction level
`sbatch qiime-tables.slurm` 

### 4. Assign taxonomy


These scripts train the classifier using sklearn, 


`sbatch unite_train_classify.slurm` for ITS


or  


`sbatch silva_train_classify.slurm` for 16S

### 5. Build tree and perform core diversity metrics
`qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ./rep-seqs.qza --o-alignment ./aligned-rep-seqs.qza --o-masked-alignment ./masked-align-rep-seqs.qza --o-tree ./unrooted-tree.qza --o-rooted-tree ./rooted-tree.qza`


`qiime diversity core-metrics-phylogenetic --i-phylogeny ./rooted-tree.qza --i-table ./table.qza --p-sampling-depth 10100 --m-metadata-file ./metadata_its2020.tsv --output-dir ./core-div-metrics`


### 6. Export files for analysis in R
`qiime tools export --input-path rarefied_table.qza --output-path ./phyloseq/`

`biom convert -i feature-table.biom -o feat_table.txt --to-tsv`

`qiime tools export --input-path its_classification.qza --output-path core-div-metrics/phyloseq`

`qiime tools export --input-path unrooted-tree.qza --output-path core-div-metrics/phyloseq`


### 7. rest of pipeline in R
