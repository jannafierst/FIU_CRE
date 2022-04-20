# FIU_CRE
Comparative genomics, from annotation onwards (2022)

*This is a work in progress, edited from a 2021 project started at the University of Alabama*

## PART 5: Annotation

We have assembled sequences (ACGT etc) but how do we work with them? The first step is annotating genome features including transposable elements and protein-coding genes.

We will be annotating sequences that were extracted, assembled and decontaminated by last year's CRE cohort. You can access them at the NCBI by searching for "Oscheius" and identifying the University of Alabama deposited sequences. Here, I will work you through one example - Oscheius dolichura.

Obtain the assembled sequences here https://www.ncbi.nlm.nih.gov/assembly/GCA_022343505.1 through the side menu that says 'FTP directory for GenBank assembly.' You can click on this link then right click on the file name that ends with .fna. This is the NCBI notation for assembled genome sequences. Go back to your FIU HPC window, type

	$ wget 
	
	And paste in the link you copied when you right-clicked on the file name. It should start a download and obtain the assembled sequence. It is transferred as a compressed (zipped) archive, unzip it like this
	
	$ gunzip GCA_022343505.1_ASM2234350v1_genomic.fna.gz 
	
	And you are ready to work.
    
### 5.1 Characterizing Repeats and Transposable Elements

### 5.1.1 [RepeatMasker](http://www.repeatmasker.org/)

RepeatMasker is a clearinghouse for repeat-associated sequences in a variety of organisms. Here, we will use the RepeatMasker tools to create 

	Use the queryRepeatDatabase.pl script inside RepeatMasker/util to extract Rhabditida repeats
	
	$ module load RepeatMasker-4.1.0

	$ queryRepeatDatabase.pl -species rhabditida | grep -v "Species:" > Rhabditida.repeatmasker

Finally, we will use RepeatMasker to create a masked version of our assembled genome sequence. There are two versions of masking. Hard-masking means we replace every nucleotide in a repeat region with 'N' and soft-masking means we replace the normally capitalized nucleotides with lower-case nucleotides in repeat regions. Here we soft-mask (-xsmall) and do not mask low complexity elements (-nolow). This takes some time and we will use a slurm script to send the computations to a node.


	Create a new file with vim
	
	$ vi repeat_slurm.sh
	
	Once inside vi type 'i' for insertion mode and enter the following:

	#!/bin/bash

	#SBATCH --qos=pq_bsc4934-5935
	#SBATCH --account=acc_bsc4934-5935

	# Number of nodes
	#SBATCH -N 1

	# Number of tasks
	#SBATCH -n 8

	#SBATCH --output=log

	##########################################################
	# Setup envrionmental variable. 
	##########################################################
	export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
	. $MODULESHOME/../global/profile.modules
	
	module load RepeatMasker-4.1.0

	RepeatMasker -lib Rhabditida.repeatmasker -pa 8 -xsmall -nolow [genome.fasta] 
	
	```
	Once you have finished this, type <escape> to enter command mode and then :wq to save your script and exit vi.
	
	Submit this to the FIU HPC with
	
	$ sbatch repeat_slurm.sh
	
	```
	End for today! Great work team :)
	
We can explore our RepeatMasker results a little, file by file. Use 'more' to view the contents and we can use grep (Global regular expression print) to search for text patterns. Here, we have asked RepeatMasker to produce a soft-masked file so we should have DNA characters (ACGT) in lower case (agct). We can, for example ask grep to count (-c)

	$ grep [acgt] <your_genome.masked>

Where <your_genome.masked> is the name of your masked .fasta file. We are asking grep to search for any of the characters within the brackets [] and we could search instead for any digit with [0-9] or any capital letter with [A-Z].  

grep has a number of different options, view them with

	$ man grep
	
To search for soft-masked sequence we can use a more compact option, -c. This asks grep to count

	$ grep -c [acgt] <your_genome.masked>

Importantly, grep is searching per line; to get an accurate count we can use fold to create single-character lines and then count lines with lowercase characters

	$ fold -w 1 <your_genome.masked> | grep -c [act]
	
Remember the '|' operator from our first lectures? Try to experiment and see what happens when you use each command separately

	$ fold -w 1 <your_genome.masked>
	
	$ fold -w 2 <your_genome.masked>

Did you produce a soft-masked sequence?
	
### 5.1.2 [EDTA](https://github.com/oushujun/EDTA)

	EDTA requires a list of coding sequences (CDS) so we will need to create protein-coding gene annotations then use these in characterizing repeats with EDTA.

### 5.2 Protein-coding gene annotation with BRAKER2

#### 5.2.1 Align RNASeq with [STAR](https://github.com/alexdobin/STAR)

```
We need to generate a STAR genome index and align our RNA-Seq with STAR.
Challenge for today: can you create your own script to do this, given the file above and the commands below?
The STAR module on FIU HPC is called star-2.7.9a and our RNA-Seq files are in /scratch/classroom. May the force be with you!

# Generate genome index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir [species_dir] --genomeSAindexNbases 12 --genomeFastaFiles [species_genome]

# Map the reads
STAR --runThreadN 12 --genomeDir [species_dir] --outSAMtype BAM Unsorted --twopassMode Basic --readFilesCommand zcat \ # if reads are zip-compressed
--readFilesIn [File_1.fastq.gz] [File_2.fastq.gz] 
```

#### 5.2.2 Run [BRAKER](http://exon.gatech.edu/genemark/braker1.html)

BRAKER is not installed on the FIU HPC and we will obtain it through miniconda. Conda (Anaconda) is a software management system that creates environments specific to individual software package needs. `miniconda' is a smaller version of conda.

First, we will need to load the FIU HPC miniconda 

	$ module load miniconda3-4.5.11-gcc-8.2.0-oqs2mbg

Activate your conda environment

	$ source activate braker2

We will need to create a directory for our species manually because of the FIU HPC environment.

	$ perl ~/.conda/envs/braker2/bin/new_species.pl --species=Odolichura
	
Copy the GENEMARK package and key from /scratch/classroom to your home directory

	$ cp /scratch/classroom/gm* ~
	
Use the script below to test your installation (note it is a bash script not a slurm job script so we can't submit it to the compute nodes).

Type vi braker.sh to create a file then 'i' for insertion mode and enter:

	#!/bin/bash
	
	braker.pl \
	--workingdir=~/braker \ # your directory with files
	--species=Odolichura \
	--cores=8 \
	--genome=GCA_022343505.1_ASM2234350v1_genomic.fna.masked \
	--bam=Aligned.out.bam \
	--GENEMARK_PATH= ~/gmes_linux_64_4 \
	--BAMTOOLS_PATH=/home/applications/spack/applications/gcc-8.2.0/bamtools-2.5.1-g5liqqe65vieinsy7unrb5lu5bflxtkx/bin \
	--softmasking \
	--useexisting
	
Hit <escape> to enter command mode and type :wq to save your file and quit vi.

Type

	$ bash braker.sh
	
And see what happens. I am actually not sure what will happen! We are still testing the FIU HPC installation. 

If it does work we run braker2 in 2 steps

First, with RNA-Seq evidence only (above)

Then with protein evidence only

	### BRAKER2
	braker.pl \
	--workingdir=[output_dir] \
	--species=[species_name] \
	--cores=8 \
	--genome=[assembly.masked] \
	--prot_seq=proteins.fa \
	--softmasking \
	--epmode
	```

If it doesn't work we move on to...

### 5.4 Functional annotation with Interproscan

Now that we have a set of protein-coding gene annotations, what can we do with them?

We still want to know what our genes and proteins may actually do and so for that we want to perform functional annotation. We will use the software Interproscan, it puts together a number of different approaches to identify protein domains with similarity to proteins with known functions. It will also associate any gene ontology terms (GO terms) and pathway information (KEGG pathways, Reactome). It is available on the FIU HPC as a module

interproscan-5.55

Your 'challenge' task for today is to generate scripts to run functional annotation on two sets of proteomes (proteins generated from protein-coding annotations) available in the /scratch/classroom directory. They are:

caenorhabditis_remanei.PRJNA577507.WBPS16.protein.fa and 356.protein.fasta

You can use your previous scripts as templates to:

1) specify your slurm (#SBATCH) configuration options
2) load the module
3) generate interproscan annotations including GO terms and pathway information. Here is a sample command to do this:

interproscan.sh -i 356.protein.fasta -f tsv -dp -goterms -pa

Here, I am asking interproscan (interproscan.sh) to perform functional annotation on my proteome (356.protein.fasta) including gene ontology annotations (-goterms) and pathway information (-pa). I am asking for a tab-separated output file (-f tsv) and -dp is disable precalculated lookup service, we will calculated everything fresh for this genome.

### 5.5 Querying our functional annotations

We completed our interproscan searches and have lists of genes with annotated protein domains. But what do these mean? We will first go through an example in class and discuss how to assign meaning to our results.

We can use Unix tools to create compiled lists, for example:

	$ cat braker.356.protein.edit.tsv | cut -f 1,12,13,14,15,16 | grep 'IPR' | sort | uniq | cut -f 2,3 | sort | uniq -c | sort -nr | tr -s '[:blank:]' | head -100

Spend some time playing with this line and figuring out what each modular command does.

Even with protein-coding annotations and functional annotations it's clear we still need more information! For that we will use BLAST and discuss today what BLAST is and how it works.

### 5.6 BLAST on the command line

The NCBI Blast web interfaces (and other web interfaces) are convenient for single sequence queries but if we have many queries we want to deal with it quickly becomes tedious. Using BLAST on the HPC allows us to scale up our work and easily format the output. 

Part of my research is studying how genes and proteins vary within populations and linking this variation to larger evolutionary patterns of change and conservation. We will take our annotated protein-coding genes and try to characterize near relatives in the model nematode *C. elegans.* * C. elegans* is one of the best-studied metazoan organisms and we know a lot about its genes and proteins.

In order to identify the nearest relatives of each protein and the sequence-level changes that have occurred between the relatives we will use BLAST to identify orthologous (i.e., descended from the same ancestral sequence) proteins between *C. remanei* and *C. elegans.* We can also use the BLAST and interproscan results to evaluate our annotation results. For example, if <50% of our predicted proteins have functional annotations or orthologous proteins it's likely our annotations are poor quality and we need to work on fine-tuning our protocol.

The first step is obtaining the set of *C. elegans* proteins. Navigate to WormBase Parasite https://parasite.wormbase.org/index.html and click on the 'Downloads' page. Right click on the link for FASTA under Proteins and click ‘Copy link location.’ Now go back to your HPC terminal and type

    $ wget ***PASTE COPIED LINK IN HERE***

Don’t actually type the part above in bold; paste the copied link in. You should have something like

    $ wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS16.protein.fa.gz

wget and curl are two handy Unix programs that can help us download data from public sites. Once you hit return you should connect to the WormBase Parasite ftp site and download the file. It is zipped, so unzip it with 

    $ gunzip *gz

Memory jog - using our wildcard \* character here. 

Take a look at the file and do some sanity checks. Use disk usage to check how large it is

    $ du -h *fa

and check how many lines and miRNA sequences you have. 

    $ wc -l *fa
    
    $ grep -c ">" *.fa

Here we are asking grep to count (-c) every instance of the .fasta header '>'. We can also get more specific and ask for lines that start with '>' with a '^' character. For regular expressions '^' denotes start of line and '$' denotes end of line. For example, try

    $ grep -c "^>" *.fa
   
and 

    $ grep -c ">$" *.fa
    
In the first command you are counting how many lines start with '>' and in the second command you are counting how many lines end with '>'. 

We can now use our *C. elegans* file to create our own BLAST database. Before we set up our script it is helpful to look at a few parameters so we know what we are setting up. Load the blast module

    $ module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k

and type

    $ module list

to make sure it loads. Type 

    $ blastp — -help

to see a list of the parameters blast requires and/or accepts. Unix BLAST has high functionality but limited help. For example, we need to make a blast database with ‘makeblastdb’ but we can not actually get UNIX BLAST to give us these kind of instructions in the Unix system. 

Why do we want to make our own database? First, the default BLAST databases don’t work well when we have very specific queries. Second, it will be much faster to run BLAST if our database is smaller (less to search) and the database is located in the same location as our working files (no input/output bottlenecks). 

In order to see the parameters required and/or accepted in making a BLAST database type

    $ makeblastdb

What happens?

Try 

    $ makeblastdb -h

What happens?

Type 

    $ makeblastdb -in ./caenorhabditis_elegans.PRJNA13758.WBPS16.protein.fa -dbtype prot -title CELEG -out CELEG

You should see a series of files with CELEG.XXX names produced. These are the index files BLAST uses to accelerate database searching. Note that we were able to do this on the head node because it is a small file and runs quickly; we need to run the BLAST search through a script. Before we make our script look at the blastp options again. The most significant are ‘db’ for which database to search and ‘query’ for which file we want to search with.  

I have created a shortened protein file so we can quickly evaluate our BLASTP output. Copy it from the /scratch/classroom to your home directory:

	$ cp /scratch/classroom/356.protein.fasta.short ~

Use vi to create a new script

    	$ vi blast_script.sh
    
Type 'i' for insertion mode and insert the following text:

    #!/bin/bash

	#SBATCH --qos=pq_bsc4934-5935
	#SBATCH --account=acc_bsc4934-5935

	# Specify number of unique nodes/servers
	#SBATCH -N 1
	# Specify the number of cores on a single node
	#SBATCH --cpus-per-task=8
	#SBATCH --threads-per-core=1

	#SBATCH --output=log

	##########################################################
	# Setup envrionmental variable. 
	##########################################################
	. $MODULESHOME/../global/profile.modules
	
	module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k
	export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
	pwd; hostname; date

	scontrol show hostname $SLURM_JOB_NODELIST > `pwd`/machines

	##########################################################
	##########################################################
	
    blastp -db CELEG -query 356.protein.fasta.short > blast_output.txt

Save your script by hitting ‘escape’ and then ‘:wq’ to exit vi. Make the script executable:

    	$ chmod +x blast_script.sh
    
and submit the script:

    	$ sbatch blast_script.sh
	
What happens?

It is sometimes easier to see what’s going on with a tabular output. Edit your script so the Blast line looks like this and run again:                         

    blastp -db CELEG -outfmt 6 -query 356.protein.fasta.short > blast_output_tabular.txt

What happens now? The columns tell you:

1) Query

2) Blast hit

3) Percent match

4) Length of aligned sequence

5) Number of polymorphisms (mismatches)

6) Gap score 

7) Query start 

8) Query end 

9) Subject start

10) Subject end

11) e-value

12) bitscore

What do you see in the results? Make sure you are running with and without the tabular format (-outfmt 6) so you can visualize the alignments as well as the matches.

One thing you can see is the multiple alignments for each query. What if you want to limit it to the first, best hit? BLAST has parameters for that but it turns out they don't work! (See: https://academic.oup.com/bioinformatics/article/35/9/1613/5106166?login=true )

We can use awk, a Unix text utility, for file manipulation. For example:

	$ awk '!_[$1]++' $FILE # here, selecting the first entry of a set of duplicates based on the value in column 1

	# calculate a column average

	$ cat $FILE | awk '{sum=sum+$3} END {print sum/NR; sum = 0}' # Here, calculating an average for column 3

	$ awk '$11 < 0.00001' $FILE # selecting all lines where column 11 is less than 0.00001 












### 6 Annotation statistics with AGAT

AGAT(https://github.com/NBISweden/AGAT#installation) is a tool for annotation editing and evaluation. We will install via conda and use it to evaluate annotation statistics. AGAT creates conflicts with some other aspects of conda and we will install/activate it into its own environment to manage the conflicts.

#### 6 Install AGAT via conda(https://anaconda.org/bioconda/agat)

	conda create -n agatenv
	conda activate agatenv
	conda install -c bioconda agat

#### 6 Count genes and other features

	conda activate agatenv
	agat_sp_statistics.pl --gff {file}.gff3


#### 6 Run Tsebra with the braker2 RNA-Seq and protein outputs

After braker2 is finished we will run Tsebra to integrate the two types of evidence

	./bin/tsebra.py -g braker1_out/augustus.hints.gtf,braker2_out/augustus.hints.gtf -c default.cfg \
    	-e braker1_out/hintsfile.gff,braker2_out/hintsfile.gff \
    	-o braker1+2_combined.gtf

```
The native output of BRAKER2/Tsebra is .gtf. We need to get it into amore usable .gff3. We will use agat to create a .gff3 and extract the mRNA and protein sequences from the .gff3 and assembled genome sequences.

	#!/bin/bash

	DIR="[BRAKER2 .gtf location]"

	conda activate agatenv

	agat_sp_extract_sequences.pl -f ${DIR}_filtered.fasta --mrna -g braker.gff3 -o braker.mRNA.fasta
	agat_sp_extract_sequences.pl -f ${DIR}_filtered.fasta -p -g braker.gff3 -o braker.protein.fasta

	tblastn -db nt -query braker.protein.fasta \
	-outfmt '6 qseqid qlen staxids bitscore std sscinames sskingdoms stitle' \
	-num_threads 4 -evalue 0.01 -max_target_seqs 2 -out blast.out
	```


