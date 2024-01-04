# GLiDe: Guide Library Designer (offline package)

## What is this?
This python script collection is the offline version of [GLiDe](https://www.thu-big.net/sgRNA_design/), a publicly avaliable web server used for genome-wide or focused sgRNA library design for CRISPRi systems. It is user-friendly for experimental biologists with no or limited programming expertise. The user only need to configure several parameters and optionally upload two standard parameters to design an sgRNA library for CRISPRi systems in a particular microorganism. The basic description of this program can be found at [BioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.25.517898v3). Please cite this paper or subsequent peer-reviewed publication if this program is useful to your work.

To use this offline version, the user only need to download the main script and several standard files (genome and annotation), edit a configure file to set several parameters needed for sgRNA design, and type in one command line (for example, cmd in Windows or terminal in MacOS) to initiate the design process.

## How to use it?
### Step 1: Install the necessary packages.

1. Install Python version 3.6.8 or above
2. Install Numpy version 1.19.1 or above
3. Install Pandas version 1.1.0 or above
4. Install RegEx version 2.5.83 or above
5. Install SeqMap. Please go to the [official cite](http://www-personal.umich.edu/~jianghui/seqmap/download/seqmap-1.0.13-src.zip) and download the seqmap-1.0.13-src.zip. Please follow the instructions of SeqMap. Briefly, unzip the file and open the command line window to use the cd command to the fold path and input "g++ -O3 -m64 -o seqmap match.cpp" or "g++ -O3 -m32 -o seqmap match.cpp" based on your computer system. After that copy one seqmap executable file to your working directory.
6. Install blast. Please go to the [official site](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and download blast version 2.15.0 or above. For Windows users, after installation, please add your_installation_path/NCBI/blast-2.15.0+/bin to PATH.

### Step 2: Prepare the necessary files. 

For one microorganism, two standard files are needed:

1. **Sequence File**

Sequence file should be single contig genome file (must in FASTA format, .fna or .fasta file eg.). In addition to the four bases (A, C, G and T), GLiDe also accepts mix-bases symbol (R, Y, M, K, S, W, H, B, V, D and N), other characters, like I (Hypoxanthine) or U (Uracil) are not accepted. Both upper and lower cases are acceptable.

The header of the sequences should include accession number and name, separated by space and leading by the ">" symbol. An example of a header structure is: ">NZ_LR881938.1 Escherichia coli str. K-12 substr. MG1655 strain K-12 chromosome MG1655, complete sequence". If the annotation file is in .ptt or .rnt format, the sequence file should only contain a single sequence contig with its corresponding header. However, if the annotation file is in .gff or .gff3 format, the sequence file can contain multiple sequence contigs, and their headers should have accession numbers that match those in the annotation file.

2. **Annotation File**
Regarding the annotation file format, GLiDe accepts the General Feature Format (GFF/GFF3 eg.), as well as two older versions: the protein table file (PTT eg.) and the RNA table file (RNT eg.). When using the PTT or RNT formats, it's important to ensure that the "Design Target" parameter aligns with the file type uploaded (choose "CDS" for PTT format and "RNA" for RNT format).

The annotation file can be customized for tailored library design, allowing users to delete specific sequences from standard files in order to design an sgRNA library for selected regions.

If a user only has the sequence file and lacks an annotation file, they can obtain one through a standard genome annotation pipeline such as the Prokaryotic Genome Annotation Pipeline (PGAP).

### Step 3: Set up the configure file (see example_configure.txt)
The configure file is used to set all the necessary parameters and tell the program where to find necessary files. This file that contains a header is in a two-column format using colon (:) as delimiter. Each line starts with one word (name of one parameter) separated with the following (setting of this parameter) by a colon delimiter. We describe each parameter as below.

**reference_file**: the gene sequence annotation file and the 

**off_threshold**: the off target penalty threshold (default=20), sgRNAs with potential off-target site carrying penalty score lower than the threshold will be eliminated. For the detailed description of the scoring method, please check our paper. Briefly, we suggest off_threshold >= 20 for library design. In situations where more sgRNAs are desired, the threshold can be decreased to 10, where the off-target effect of CRISPRi is still very slight as previously reported (Gilbert Luke et al., Cell 2014).

**GCcontent_min**: The minimal GC content of spacer region (percentage) (default=30, positive integer between 0~100 is accepted).

**GCcontent_max**: The maximum GC content of spacer region (percentage) (default=80, positive integer between 0~100 and > GCcontentMin is accepted). 

In previous reports about dCas9 based CRISPRi system, GC content of sgRNA spacer region is found to be correlated with sgRNA activity. Extreme GC content reduces sgRNA activity. Hence, we suggest the abovementioned threshold. In situations of genome with relative low or high GC content, we suggest to adjust the threshold to (10,90). 

**target**: the target for the designed sgRNA library, which can be chosen from either the coding sequence (cds) or RNA coding genes (RNA).

DNA sequence file for genes of interest (.ffn file of protein-coding genes and .frn file of RNA-coding genes for genome-scale sgRNA library design, or your customized file for focused sgRNA library design, see Step 2)

**strand**: whether the sgRNA is designed targeting (binding) to the template or nontemplate stand of a coding gene (default=nontemplate, nontemplate or template is accepted). In the case of common CRISPRi systems, the non-template strand is preferred for effective gene silencing. However, we offer the alternative choice for applications that do not have a strand preference **better activity when targeting to non-template strand in the coding region**.

**nc_number**: the number of negative control sgRNAs (sgRNA with no significant target across the genome, which is used as negative control in the following pooled screening experiment and data analysis) for the experiment(enter 0 if not needed). We strongly recommend to include the negative control sgRNAs. For the description of negative control sgRNA usage, see our paper. The number of negative control sgRNA you want to design for the experiment. If negative option is no, select 0 for this option. The default is 400. We recommand min(400, 5% of sgRNA library size) as the number of netagive control sgRNAs. 

Below is **an example configure file with default parameters**.

parameter|value
---------|-----
[configdesign]
reference_file:|example.gff,example.fna
off_threshold:|20
GCcontent_min:|30
GCcontent_max:|80
spacer_length:|20
target:|cds
strand:|nontemplate
nc_number:|400

After Step 2 and 3, check your working directory. It should looks like below:
[here](./image/files_prepared_before_library_design.png)

### Step 4：Run the pipeline
Open the command line window, cd to the working directory and run the analysis pipeline.

cd path_to_your_working_directory

python sgRNA_desgin.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 

python sgRNA_desgin_main.py example_configure.txt

For a typical laptop, the example test can be finalized within 5 minutes. The rate-limiting step is off-target site identification across the genome.

## Output description
The output files will be organized in the subdirectory whose name is specified by the 'prefix' option in configure file under the working directory (prefiex/). We term this subdirectory 'result directory' thereafter.

You can find many files under the result directory.

[your result directory after running the test](./image/resultdir_after_example_running.png)

Below is the description. For the mathematical processing, see our paper. **All .csv flat files use tab as delimiter unless mentioned**

**prefix.N20.fasta.txt**: sgRNA library sequences (N20) in .fasta format, including negative control sgRNAs if it is specified in the configure file. **This file, after adding designed flanking nucleotides, can be subjected to microarray based oligomer synthesis** to prepare the sgRNA library. 

**prefix.N20NGG.fasta.txt**: The reverse complementary of target region of each sgRNA including the PAM site in .fasta format.

**prefix.sgRNA_statistics.txt**: position information within relevant gene coding region and the GC content of each designed sgRNA in .csv format.

sgRNAID|sgRNA_position_in_gene |GCcontent
-------|-----------------------|---------
insI1b4284_32|0.028|39.13
...|...|...

**prefix.cluster.txt**: Clusters of genes with homologs (see Step 2, "Coping with multiple copy issue" section). The file has no header line and uses tab as delimiter. It is consisted of two columns: cluster name and the genes in each cluster. Genes in one cluster are separated by comma.

araFb1901|araFb1901
-------|--------------------
yjhXb4566|yjhXb4566
tufAb3339|tufAb3339,tufBb3980
...|...

**prefix.gene_statistics.txt**: file in .csv format with one header line contains the length and the designed sgRNA numbers of each gene.

gene_name|gene_length|sgRNA_number_in_gene
---------|-----------|--------------------
ydcCb1460|1137|0
insI1b4284|1152|10
pinRb1374|591|9

**prefix.fasta.txt**: target gene sequences in .fasta format. It is similar to the input gene fasta file. For genome-wide sgRNA library design, we use the combination of gene and synonym (.ptt or .rnt annotation file) to rename each gene. Hence, we give the refined .fasta file with new gene names for some convenience in following usage.

**N20_library.csv**: the sgRNA library file is at .csv formate **containing one header line**, in which there are three columns in order of id, sequence and gene respectively. **Use comma as delimiter**. If negative control (NC) sgRNAs are within this synthetic library (specified in configure file), name them NCx and assign '0' at 'gene' column of these sgRNAs. This file is used as input for the data processing subpackage after pooled screening experiment and NGS.

id|sequence|gene
--|--------|----
sgRNA1|ATCCCCCCCCCCGGGGG|recA
NC1|TGTGTGTGTGTGTGTGTGTG|0
...|...|...

**sgRNA_position.txt**: Flat file of sgRNA position (relative location of sgRNA in the coding region) information in gene **without header line**. The file contains three columns in order of gene name, sgRNAid and the relative position of sgRNA in the gene. Actually, you can also find this file as output of our library design subpackage. This file is also used as input for the data processing subpackage after pooled screening experiment and NGS.

rsmE|rsmE_9|0.012
----|------|-----
rsmE|rsmE_10|0.014
0|NC1|0
...|...|...

Two .png figures use histogram to summarize the basic information of [number of sgRNA per gene](./image/The_distrution_of_sgRNA_number_per_gene.png) and [position of sgRNAs in the coding region of relevant genes](./image/The_distrution_of_sgRNA_position_per_gene.png).
