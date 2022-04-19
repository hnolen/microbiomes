##QDM-011 ITS Pipeline

##using qiime2-2020.2

##importing reads from genome center
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./reads/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path ./imported_reads.qza

##getting quality plots
qiime demux summarize --i-data imported_reads.qza --o-visualization quality.qzv
##trim forward at 227

#running cutadapt and itsxpress
sbatch itsXpress.slurm
