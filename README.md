# GLiDe: Guide Library Designer (offline package)

## What is this?
This python script collection is the offline version of [GLiDe](https://www.thu-big.net/sgRNA_design/), a publicly avaliable web server used for genome-wide or focused sgRNA library design for CRISPRi systems. It is user-friendly for experimental biologists with no or limited programming expertise. Users only need to configure a few parameters and optionally upload two standard files to generate an sgRNA library for CRISPRi applications in a specific microorganism. The basic description of this program can be found at [GLiDe: a web-based genome-scale CRISPRi sgRNA design tool for prokaryotes](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-06012-0). Please cite this paper or subsequent peer-reviewed publication if this program is useful to your work.

To use this offline version, the user only need to download the main script and several standard files (genome and annotation), edit a configure file to set several parameters needed for sgRNA design, and type in one command line (for example, cmd in Windows or terminal in MacOS) to initiate the design process.

## How to use it?
### Step 1: Install the necessary packages.

1. Install Python version 3.6.8 or above
2. Install Numpy version 1.19.1 or above
3. Install Pandas version 1.1.0 or above
4. Install RegEx version 2.5.83 or above
5. Install SeqMap. Please go to the [official cite](https://jhui2014.github.io/seqmap/) and download the seqmap-1.0.13-src.zip. Please follow the instructions of SeqMap. Briefly, unzip the file and open the command line window to use the cd command to the fold path and input "g++ -O3 -m64 -o seqmap match.cpp" or "g++ -O3 -m32 -o seqmap match.cpp" based on your computer system. After that please add your_installation_path to PATH.
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

both reference files can be easily downloaded from publicly avaliable databases like NCBI, KEGG, etc.

### Step 3: Set up the configure file (see example_configure.txt)
The configure file is used to set all the necessary parameters and tell the program where to find necessary files. This file that contains a header is in a two-column format using colon (:) as delimiter. Each line starts with one word (name of one parameter) separated with the following (setting of this parameter) by a colon delimiter. We describe each parameter as below.

**reference_file**: the gene sequence file and the corresponding annotation file. GLiDe accepts multiple reference files at the same time, in this case, input each sequence file right behind its corresponding annotation file, an example is "reference_file:annotation_1.gff,sequence_1.fna,annotation_2.gff,sequence_2.fna"

**off_threshold**: The penalty score is employed for quality control of sgRNAs (default=20). GLiDe employs the seed region principle and penalty scoring metrics (see our paper for detail) to evaluate off-targets. An off-target site is identified when the penalty score is less than the user-defined threshold, considering that mismatches are generally better tolerated at the 5′ end than at the 3′ end. Three regions are categorized based on their proximity to the PAM sequence.

**GCcontent_min**: The minimal GC content of spacer region (percentage) (default=20, positive integer between 0~100 is accepted).

**GCcontent_max**: The maximum GC content of spacer region (percentage) (default=80, positive integer between 0~100 and > GCcontentMin is accepted). 

In previous reports about dCas9 based CRISPRi system, GC content of sgRNA spacer region is found to be correlated with sgRNA activity. Extreme GC content reduces sgRNA activity. Hence, we suggest the abovementioned threshold. In situations of genome with relative low or high GC content, we suggest to adjust the threshold to (20,80). 

**target**: the target for the designed sgRNA library, which can be chosen from either the coding sequence (cds) or RNA coding genes (RNA). This parameter must match the format of annotation file used in step 2, "cds" for PTT format and "RNA" for RNT format.

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

### Step 4：Run the pipeline
Open the command line window, cd to the working directory and run the analysis pipeline.

cd path_to_your_working_directory

python sgRNA_desgin.py configure.txt

We also post a toy example together with the scripts and the example_configure.txt has been edit to make it compatible. For this test, cd to the working directory, type in: 

python sgRNA_desgin_main.py example_configure.txt

For a typical laptop, the example test can be finalized within 5 minutes. The rate-limiting step is off-target site identification across the genome.

## Output description
You can find many files under the result directory, two excels are the final results: gene_info.xlsx and offtarget_log.xlsx.

Below is the description. For the mathematical processing, see our paper.

**Output 1: Guide Library List**
The main output of GLiDe is a list containing all sgRNAs (in Excel format). In the list sgRNAs are classified with their targeted genes and labeled with their start positions.

Seqid|Source|Type|Start|End|Score|Strand|Phase|Attributes|Length|Guide_seq|Guide_pos|guide_num
-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----
NC_000913.3|RefSeq|CDS|190|255|.|+|0|ID=cds-NP_414542.1;Parent=gene-b0001;Dbxref=UniProtKB/Swiss-Prot:P0AD86,Genbank:NP_414542.1,ASAP:ABE-0000006,ECOCYC:EG11277,GeneID:944742;Name=NP_414542.1;gbkey=CDS;gene=thrL;locus_tag=b0001;orig_transcript_id=gnl[b0001]mrna.NP_414542;product=thr operon leader peptide;protein_id=NP_414542.1;transl_table=11|21|['GGTGGTGCTAATGCGTTTCA', 'CGCACCGTTACCTGTGGTAA']|[210, 249]|2
...|...|...|...|...|...|...|...|...|...|...|...|...

For each gene, sgRNAs are ranked by their distance to the start codon, those closer to transcription start sites are at the top. Genes are sorted in their natural order.

Negative control sgRNAs (if designed) would be placed in the last row of the list.

This table is temporarily saved on the server, user has to download the Excel form for inspection.

**Output 2: Off-target Log List**
This list contains sgRNAs predicted to have potential off-target hits. All unselected sgRNAs with their predicted off-target sequences would be listed inside. This table is also temporarily saved on the server, user has to download the Excel form for inspection.

ID|sequence|offtarget_sequence|penalty|off_num
---|---|---|---|---
5704|AAACACCAGTTCGCCATTGC|['CGAATACAGTACGCCATTGC', 'GAAGAGCTTTCCGCCATTGC']|[17.0, 19.0]|2
...|...|...|...|...
