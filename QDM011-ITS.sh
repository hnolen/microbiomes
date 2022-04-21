##QDM-011 ITS Pipeline

##using qiime2-2020.2

##importing reads from genome center
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./reads/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path ./imported_reads.qza

##getting quality plots
qiime demux summarize --i-data imported_reads.qza --o-visualization quality.qzv
##trim forward at 227

#running cutadapt and itsxpress
sbatch itsXpress.slurm #output file xpress_trimmed_reads.qza
##note: itsXpress actually worked with these reads, when I did ITS for 2019 samples I kept getting
## a non zero exit status error during the vsearch step - going to keep going but maybe I should
##not do itsxpress with these reads since I didn't do it on my 2019 reads

##running dada2 on xpress_trimmed_reads.qza
sbatch qiime-dada.slurm ##ran successfully, filtering removed all reads from some samples see dada2.output file

##now run qiime-tables.slurm to get the dns.qzv, table.qzv, and rep-seqs.qzv visualization files
sbatch qiime-tables.slurm ##looking at dns.qzv - decided to rarefy at 10100

##running taxonomy
sbatch unite_tax.slurm

##core diversity metrics
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ./rep-seqs.qza --o-alignment ./aligned-rep-seqs.qza --o-masked-alignment ./masked-align-rep-seqs.qza --o-tree ./unrooted-tree.qza --o-rooted-tree ./rooted-tree.qza 
qiime diversity core-metrics-phylogenetic --i-phylogeny ./rooted-tree.qza --i-table ./table.qza --p-sampling-depth 10100 --m-metadata-file ./metadata_its2020.tsv --output-dir ./core-div-metrics

