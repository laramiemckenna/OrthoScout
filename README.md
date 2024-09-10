# OrthoScout
OrthoScout is a pipeline for identifying candidate genes in target species via orthology analyses with genes of known or suspected phenotypes in model systems. It is built to work alongside an SQLite3 database called FloralCoreDB for use in identifying potential candidate genes in regions of reduced recombination such as Sex-Determining Regions (SDRs), apomixis loci, and Self-Incompatibility (SI) loci in species that lack a mechanism of functional verification [Figure 1].

This pipeine is a "funnel" to narrow down the possible genes in a region-of-interest to the most relevant based on mutant phenotypes in model species so that additional analyses may be done such as: Ks/Ka, multiple sequence alignments, structural inference and comparison, and more.

<p style="text-align: center;font-size: 20px;">Figure 1, A Graphical Representation of an SDR</p>

![](/assets/sex_chromosome_graphic_figure_1.png)

 **Note to Caroline:** The database described in "What is FloralCoreDB?" in the dropdown below refers to what is now called `floralcore_outdated.db` in the `databases` directory. Although it is used for the tutorial, it will not be used in our work because limiting the floral genes to the implied mutant phenotypes (IMP) evidence type yielded too little information to ultimately be useful. We now incorporate all evidence types in the updated database. You'll be using the `csv` version of this newest one, `floralcore_arabi_updated.csv`, as I continue to work on the SQLite database for it in the background.

<details>
<summary style="font-size: 24px;"><b>More about FloralCoreDB </b></summary>

Sex chromosomes contain regions of non-recombination called Sex-Determining Regions (SDRs) that contain an array of genes that orchestrate the development of a dimorphic sexual phenotype. These SDRs can range from 50 kbp to 100+ Mbp. Up until now, researchers have had to manually curate a list of candidate genes from these regions, tediously running BLAST on each gene and searching the internet for any connection it might have to their phenotype-of-interest.


Not only does this method take a lot of the researcher's time, but it lacks rigor, since BLAST alone is unable to detect remote homologs, has difficulty handling multi-domain proteins, and more. Additionally, sequence similarity is not a guarantee of functional similarity, so BLAST may return false positives.


In contrast, programs like OrthoFinder2, which this pipeline uses, overcomes many of these challenges through orthological inference where phylogenetic context is heavily weighted [1].


The focus of the **OrthoScout** project is to create the tools necessary to utilize a database I am building for my dissertation called **FloralCoreDB**.<br>
 
> :cherry_blossom: **What is FloralCoreDB?** :cherry_blossom:
>
>
> This SQLite database will include information about genes involved in plant reproductive development across Angiosperms.
>
>The purpose of **FloralCoreDB** is to build (and frequently update) a set of genes for querying the *de novo* assembled genomes and sex-determining regions of dioecious plants for aid in the discovery of genes involved in flowering and the control of sex in plants.
>
>We are starting with floral genes from implied mutant phenotypes (IMP) taken from Arabidopsis. The basic structure of the database resembles a TAIR output, but additional columns, such as GO IDs and GO SLIM IDs have been added. KEGG pathway integration is the next step for this database as well as the integration of maize and rice floral mutants.
</details><br>


## Installation
---
To install and begin running OrthoScout, follow these steps:

1. Clone the repository with `git clone`
2. Navigate into the repository directory: `cd OrthoScout`
3. Create a virtual environment: `python -m venv orthoscout_env`
4. Activate the virtual environment: `source orthoscout_env/bin/activate`
5. Install the required packages: `pip install -r requirements.txt`

If you do not have a peptide file for your target species and plan on using `run_tblastn.py`, then you will also need to install NCBI's command line tools (the NCBI C++ Toolkit). The best way to do this is through [Conda](https://anaconda.org/bioconda/blast).

#### Running OrthoScout using Singularity on HudsonAlpha's cluster: TBD

## System Requirements
---
To run OrthoScout, you need Python 3.8.8 and the following packages, which can be installed as indicated above:

- biopython==1.81 - provides tools for working with biological data
  - Bio
    - SeqIO
    - Entrez
    - Bio.SeqRecord
- pandas==2.2.2 - data manipulation and analysis

Modules already included the Python Standard Library:
- argparse - processes command-line arguments
- os - interact with the operating system
- typing - used to improve function annotations
- csv - allows user to read and write csv files
- re - operations for regular expression
- pathlib - used for filepath manipulation
- subprocess - used to run external programs and then interct with them via python

## Usage
---
<p style="text-align: center;font-size: 20px;">Figure 2, A Flowchart of the OrthoScout Pipeline</p>

![](/assets/corrected_flowchart_orthoscout.png)

Depending on the option that you choose, the following files need to be together in their own directory (just like the demonstration_data directory):
- the SQLite3 database (Option #1 & #2)
- the target species peptide file (Option #1)
- the Orthogroups.tsv file generated by OrthoFinder2 (Option #1)
- the target genome .gtf file (Option #1)
- the target genome fasta (Option #2)

This Python-only pipeline utilizes seven scripts that offers two analysis paths [Figure 2]: 

### **Option 1:** If peptide fasta files (Protein FASTA on NCBI) are available via [LiftOff](https://github.com/agshumate/Liftoff) [2] or a complete annotation:

<br>

>Note: Prior to running `orthoscout.py` you will need to run [OrthoFinder2](https://github.com/davidemms/OrthoFinder) using the peptide files from the FloralCoreDB organisms, your target organism, and a handful of other organisms (to allow for more robust phylogenetic context). In the future, a directory with the peptide files for the FloralCoreDB as well as a handful of other outgroups to improve rigor will become available.

<br>

The script `search_and_extract.py` allows you to search for a pattern in a specific column of an SQLite database and extract the matching rows to a CSV file. The usage of the script is as follows:

```
usage: python search_and_extract.py [-h] --pattern PATTERN --column-name COLUMN_NAME --database-path DATABASE_PATH

Search an SQLite database for rows that match a pattern in a specified column, and write the matching rows to a CSV file.

optional arguments:
  -h, --help            show this help message and exit
  --pattern PATTERN, -p PATTERN
                        The pattern to search for, enclosed in quotes - can be a word or exact phrase.
  --column-name COLUMN_NAME, -c COLUMN_NAME
                        The name of the column to search in.
  --database-path DATABASE_PATH, --db DATABASE_PATH
                        The path to the SQLite database file.
```

Add in the `PATTERN` you want to search for (i.e. carpel, fertility, etc.), `COLUMN_NAME` with the name of the column you want to search in, and `DATABASE_PATH` with the path to the SQLite database file.

The script will create a new CSV file in the same directory as the database file, with the name [pattern]_rows.csv, containing the rows that match the specified pattern. (At the moment, this subsetted CSV directly resembles the SQLite3 database, so if you are curious about these types of columns available, you can glance at the CSV for more information. A separate descriptive page of FloralCoreDB will be provided when it is ready for release.)


If you would prefer to use all of the loci in the database rather than a specific subset, you can use the following line to convert the database into a csv file for further processing: 

```
usage: python convert_db_to_csv.py [-h] --db DB

Convert an SQLite database to a single CSV file.

optional arguments:
  -h, --help      show this help message and exit
  --db DB, -d DB  Path to the SQLite database file
```

Next, run the script `orthoscout.py`, which allows you to process the OrthoFinder TSV, the CSV generated from the SQLite3 database, and the target protein fasta file to extract the sequences in your target species that are orthologous to the model loci you extracted using `search_and_extract.py`. The usage of the script is as follows:

```
usage: python orthoscout.py [-h] -i INPUT_FOLDER [-o OUTPUT] -t TARGET_COLUMN [-v]

Process OrthoFinder TSV, model species gene information CSV, and FASTA files.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FOLDER, --input-folder INPUT_FOLDER
                        Path to the input folder containing the OrthoFinder TSV, model species gene information CSV, and FASTA files
  -o OUTPUT, --output OUTPUT
                        Output FASTA file path
  -t TARGET_COLUMN, --target-column TARGET_COLUMN
                        Target column name in the TSV file
  -v, --verbose         Save intermediate output files in the 'verbose_output' folder
```

The `INPUT_FOLDER` argument should be the path to the folder containing the OrthoFinder TSV, subsetted CSV, and target protein fasta. The `-t `or `--target_column` option should be followed by the name of the target column in the TSV file to match with gene information. You will need to open or preview the Orthofinder.tsv to get retrieve this -- it will be in the header. The optional `-o` or `--output` option specifies the output file path (default: `target_species_orthologs.fa`). The `--verbose` option will save intermediate output files in the verbose_output folder for debugging if necessary.

The script reads the input files, matches the model locus IDs between the TSV and gene information CSV (if `--verbose`, the output of this is: `model_locus_matched_to_orthogroup_output.csv`), matches the target orthogroups (if `--verbose`, the output of this is: `target_matched_to_model_orthologs.csv`), extracts the target gene IDs, matches these with sequences in the inpu target protein fasta, and saves the final output to the specified or default output fasta file in the main directory.

A .gff and .gtf file are common outputs of most annotation programs now. The .gtf file contains information about the location of each gene and will help us narrow down on which orthologs are in our target region. To do this, run `find_orthologs_in_region.py`:

```
usage: find_orthologs_in_region.py [-h] -f FASTA -g GTF -c CHROMOSOME [-s START] [-e END] -o OUTPUT [--find_descriptions FIND_DESCRIPTIONS] [--fasta_part_index FASTA_PART_INDEX] [--fasta_separator FASTA_SEPARATOR]
                                   [--gtf_column_type GTF_COLUMN_TYPE] [--gtf_pattern GTF_PATTERN] [--postfix POSTFIX]

Find proteins in a specific range on a specific chromosome.

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Peptide fasta file.
  -g GTF, --gtf GTF     GTF file.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Chromosome.
  -s START, --start START
                        Start position.
  -e END, --end END     End position.
  -o OUTPUT, --output OUTPUT
                        Output text file.
  --find_descriptions FIND_DESCRIPTIONS
                        CSV file with gene descriptions.
  --fasta_part_index FASTA_PART_INDEX
                        Index of the part in the FASTA header (1-based index). Default is 2.
  --fasta_separator FASTA_SEPARATOR
                        Separator used in the FASTA header. Default is space (' ').
  --gtf_column_type GTF_COLUMN_TYPE
                        Type of column to extract information from in the GTF file. Default is 'transcript'.
  --gtf_pattern GTF_PATTERN
                        Regular expression pattern to extract the ID from the GTF file. Default searches anything in Column 8 (0-based index) of the gtf.
  --postfix POSTFIX     Custom post-fix to add to identifiers extracted from GTF file.
```
Where `--find_descriptions` is an optional argument that allows you to append the Model Ortholog functional description to the target genes listed in the output file. The csv input for this will be the `target_matched_to_model_orthologs.csv` located in the `verbose_output` directory. 

You can search an entire chromosome (note: you must know your chromosome of interest ahead of time), since the `-s` and `-e` arguments are optional. 

This function will return a text file summarizing the orthologs found in your area of interest, their specific coordinates, and their description.

<br>
<br>

### **Option 2:** If peptide fasta files are NOT available:

Run `search_and_extract.py`

Then, run `connect_to_ncbi.py`, which allows you to search for genes in the NCBI Nucleotide database based on the locus tags available in the output of `search_and_extract.py`, retrieve their associated protein IDs, and fetch the corresponding protein sequence information from NCBI Protein database. The usage of the script is as follows:

```
usage: python connect_to_ncbi.py [-h] --email EMAIL --input-csv INPUT_CSV --output OUTPUT

Fetch protein information from NCBI for a list of genes.

optional arguments:
  -h, --help            show this help message and exit
  --email EMAIL, -e EMAIL
                        Email address to use for Entrez queries.
  --input-csv INPUT_CSV, -i INPUT_CSV
                        Path to input folder containing CSV files.
  --output OUTPUT, -o OUTPUT
                        Path to output directory.
```

Add in your email address for Entrez queries and the path to the the CSV file generated by `search_and_extract.py` with the corresponding locus tags.

The script will generate output files with the format [locus]_output.csv and merge all retrieved sequences into a single FASTA file named `combined_output.fasta` in the specified output directory.

With this output, you can use `tblastn` via NCBI's command line tools to query your target genome [3]. This method is less robust than Option #1, so be cautious when interpreting the biological significance of hits. 

As in Option #1, I've created an option where you can specify the chromosome and chromosome range for your tblastn search. To do this, run `run_tblastn.py ` by providing the path to your target genome fasta file, the peptide file generated above, and specify your chromosome of interest and the range within the chromosome.

```
usage: python run_tblastn.py [-h] -g GENOMIC_FASTA -q QUERY_PEPTIDE_FASTA -o OUTPUT_FILE -s SUBJECT_NUM -r RANGE

Run BLAST and parse output file.

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMIC_FASTA, --genomic_fasta GENOMIC_FASTA
                        Path to the input genomic fasta file
  -q QUERY_PEPTIDE_FASTA, --query_peptide_fasta QUERY_PEPTIDE_FASTA
                        Path to the input query peptide fasta file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output file
  -s SUBJECT_NUM, --subject_num SUBJECT_NUM
                        Subject/chromosome identifier
  -r RANGE, --range RANGE
                        Sequence range in the format "start-end"
```

# Tutorial

This tutorial will utilize the peptide files from the Morella cerifera genome. *M. cerifera* is a dioecious plant with a mid-sized SDR of ~7.5 Mbp on Chr08. However, for this tutorial, we will just be looking at orthologous sequences for FloralCoreDB genes that match the pattern "stamen" on the entirety of Chr08 and show example commands for generating the output of both options.

## Option 1:

IMPORTANT: decompress the zipped files in `test_data` before proceeding. 

I want to find any Arabidopsis genes with "stamen" in their mutant phenotype description. To do this, I will run the following command:

```
python3 search_and_extract.py --pattern stamen --column-name description --database-path ./test_data/floralcore_outdated.db
```
Yay! 37 rows from the original database were found containing this pattern. Some of these rows might include the same locus tag but a different gene model -- the number of rows does not necessarily correspond with the number of unique genes.

Alternatively, you can simply use the entire database for your analyses: 

```
python convert_db_to_csv.py --db ./test_data/floralcore.db
```
Now, we can run Orthoscout to retrieve the orthologs of these genes for our target species. Remember, this script will read whatever csv file is present in the `test_data` directory, so make sure only the one you want is in there.

Since it takes a while to run this pipeline on the entire dataset, we're just going to more forward with the stamen subset: 

```
python3 orthoscout.py --input-folder ./test_data --target-column Morella_cerifera.pep --output morella_cerifera_orthologs.fa --verbose
```
With the output of this, I now have orthologs of those original Arabidopsis genes I called in the first step. Using the `gtf` file generated from the annotation, I can now search for these orthologs in a particular chromosome or in a particular region on a chromosome. The output `txt` file will return the protein name, associated orthogroup, gene location, and the corresponding Model Ortholog description.  

```
python find_orthologs_in_region.py -f morella_cerifera_orthologs.fa -g ./test_data/morella_cerifera_annotation.gtf -c chr08 -o chr_08_stamen_orthologs.txt --find_descriptions ./verbose_output/target_matched_to_model_orthologs.csv --fasta_part_index 2 --fasta_separator " " --gtf_column_type transcript --gtf_pattern "(.*)"
```

***Note:*** you can generate standard format `gtf` files from `gff` files using the program [AGAT](https://agat.readthedocs.io/en/latest/gff_to_gtf.html).


---

## References:
[1] Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20, 1-14.

[2] Shumate, A., & Salzberg, S. L. (2021). Liftoff: accurate mapping of gene annotations. Bioinformatics, 37(12), 1639-1643.

[3] Command line tools. (2023). Command line tools. NIH

[4] Jia, H. M., Jia, H. J., Cai, Q. L., Wang, Y., Zhao, H. B., Yang, W. F., ... & Gao, Z. S. (2019). The red bayberry genome and genetic basis of sex determination. Plant Biotechnology Journal, 17(2), 397-409.

[5] TAIR - Home Page. (2023). Arabidopsis.org. https://www.arabidopsis.org/
