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

BRAKER is not installed on the FIU HPC and we will obtain it through miniconda. Conda is a software management system that creates environments specific to individual software package needs.

First, we will need to load the FIU HPC miniconda 

	$ module load miniconda3-4.5.11-gcc-8.2.0-oqs2mbg
	
Then create a braker2 environment for our software

	$ conda create --name braker2

Activate your conda environment

	$ conda activate braker2

BRAKER requires a lot of custom software installations. First, we will install perl modules (this script is in /scratch/classroom/)

	#!/bin/bash # Everything begins with a shebang!
	
	conda install -c anaconda perl
	conda install -c bioconda perl-app-cpanminus
	conda install -c bioconda perl-hash-merge
	conda install -c bioconda perl-parallel-forkmanager
	conda install -c bioconda perl-scalar-util-numeric
	conda install -c bioconda perl-yaml
	conda install -c bioconda perl-class-data-inheritable
	conda install -c bioconda perl-exception-class
	conda install -c bioconda perl-test-pod
	conda install -c anaconda biopython
	conda install -c bioconda perl-file-which # skip if you are not comparing to reference annotation
	conda install -c bioconda perl-mce
	conda install -c bioconda perl-threaded
	conda install -c bioconda perl-list-util
	conda install -c bioconda perl-math-utils
	conda install -c bioconda cdbtools
	
Next we will unpack the braker2 software (also available here https://github.com/Gaius-Augustus/BRAKER#supported-software-versions)

	$ cp /scratch/classroom/v2.1.6.tar.gz ~
	
	$ tar -xvf v2.1.6.tar.gz
	
Edit your .bashrc to include the BRAKER2 scripts

	$ vi ~/.bashrc
	
	Press 'i' for insert mode, navigate to the bottom of the file and add
	
	# User specific aliases and functions
	PATH=/home/[your username]/braker/BRAKER-2.1.6/scripts:$PATH
	export PATH

	Press <escape> to enter command line mode and :wq to save and quit.

We will continue with software installation on Monday. Great work!

#### End for Wednesday April 6 ####







First, with RNA-Seq evidence only

```
braker.pl \
--workingdir=[output_dir] \
--species=[species_name] \
--cores=8 \
--genome=[assembly.masked] \
--bam=Aligned.out.bam \
--softmasking

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

After braker2 is finished we will run Tsebra to integrate the two types of evidence

./bin/tsebra.py -g braker1_out/augustus.hints.gtf,braker2_out/augustus.hints.gtf -c default.cfg \
    -e braker1_out/hintsfile.gff,braker2_out/hintsfile.gff \
    -o braker1+2_combined.gtf

```
The native output of BRAKER2/Tsebra is .gtf. We need to get it into amore usable .gff3. We will use agat to create a .gff3

$

and extract the mRNA and protein sequences from the .gff3 and assembled genome sequences.

#!/bin/bash

DIR="[BRAKER2 .gtf location]"

conda activate agatenv

agat_sp_extract_sequences.pl -f ${DIR}_filtered.fasta --mrna -g braker.gff3 -o braker.mRNA.fasta
agat_sp_extract_sequences.pl -f ${DIR}_filtered.fasta -p -g braker.gff3 -o braker.protein.fasta

tblastn -db nt -query braker.protein.fasta \
-outfmt '6 qseqid qlen staxids bitscore std sscinames sskingdoms stitle' \
-num_threads 4 -evalue 0.01 -max_target_seqs 2 -out blast.out
```

#### 5.3.5 Run Tsebra with the braker2 RNA-Seq and protein outputs



### 5.4 Annotation statistics with AGAT

AGAT(https://github.com/NBISweden/AGAT#installation) is a tool for annotation editing and evaluation. We will install via conda and use it to evaluate annotation statistics. AGAT creates conflicts with some other aspects of conda and we will install/activate it into its own environment to manage the conflicts.

#### 5.4.1 Install AGAT via conda(https://anaconda.org/bioconda/agat)

	conda create -n agatenv
	conda activate agatenv
	conda install -c bioconda agat

#### 5.4.2 Count genes and other features

	conda activate agatenv
	agat_sp_statistics.pl --gff {file}.gff3

### 5.5 Create pseudo-chromosomes with Ragtag

Caenorhabditis have highly conserved chromosome structure and we expect that large-scale rearrangements have not occurred between closely related strains. This is an assumption and it may not be true but we will operate under this for now. If we assume this we can use RagTag (https://github.com/malonge/RagTag) to scaffold our fragmented assembly based on a chromosome-scale assembly for a close relative, C. remanei PX506 (https://www.ncbi.nlm.nih.gov/genome/253?genome_assembly_id=771236).

RagTag is very fast, below is a script that can align our fragmented assembly to the chromosome-scale PX506 and create a new set of coordinates for our annotated genes.

	#!/bin/bash
	
	GENOME=534

	ragtag.py scaffold GCA_010183535.1_CRPX506_genomic.fna {GENOME}.fasta -t 8

	ragtag.py updategff braker.gff3 ./ragtag_output/ragtag.scaffold.agp > ragtag.{GENOME}.gff3


## PART 6: Upload data to NCBI
