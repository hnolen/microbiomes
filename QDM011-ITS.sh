##QDM-011 ITS Pipeline

##using qiime2-2020.2

##importing reads from genome center
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./reads/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path ./imported_reads.qza



