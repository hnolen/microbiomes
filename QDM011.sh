###transferring QDM work to vscode - it seems easier to keep track of commands this way versus using notion


##logging in - working on QDM-011 16s data
ssh hbn1002@premise.sr.unh.edu

##based on interactive plots I will rarefy to 2300 

##next step is taxonomy
cp /mnt/home/poleatewich/hbn1002/QDM009/data/16s/silva-training-classifier/silva_train_classify.slurm /mnt/home/poleatewich/hbn1002/QDM011/16s

##edited silva_train_classify.slurm in the QDM011 folder for 2020 file names

pwd
sbatch silva_train_classify.slurm #just ran the classifier step output: silva_classifier.output
sbatch silva_train_classify.slurm #just ran the testing classifier step output: silva_classifier.output
##getting pop from empty list error
##running an interactive slurm session to look at the /tmp file on node117 (where the program ran)
srun --nodelist=node117 --pty bash -i 
exit
##changed the test classifier command to:
qiime feature-classifier classify-sklearn --i-classifier ./silva_trained_classifier.qza --i-reads ./rep-seqs.qza --o-classification silva_taxonomy.qza --p-confidence 1 --p-read-orientation same

sbatch silva_train_classify.slurm #just ran the taxa bar plots step output: silva_barplots.output

##got taxa barplots back and about 90% of each sample was unassigned - going to rerun taxonomy classifier step with a lower confidence
##running classify-sklearn command at --p-confidence 0
##got the pop from empty list error - changed to --p-confidence 0.5
##keep getting pop from empty list error - changed --p-confidence to 'disable'
##this worked, gave me a taxonomy file
