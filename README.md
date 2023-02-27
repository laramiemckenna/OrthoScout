# OrthoScout
Project proposal for BIOL7180-D01: Scripting for Biologists

## Project Purpose: 
The focus of the **OrthoScout** project is to create the tools necessary to utilize a database I am building for my dissertation called **FloralCoreDB**. 

> :cherry_blossom: **What is FloralCoreDB?** :cherry_blossom:
> 
>
> This SQLite database will include information about functionally-verified genes involved in plant reproductive development across Angiosperms. For this semester, the focus is only on adding genes from *Arabidopsis* and *Zea mays* in carpel and stamen organ categories.
>
>The purpose of **FloralCoreDB** is to build (and frequently update) a set of genes for querying the *de novo* assembled genomes and sex-determining regions of dioecious plants for aid in the discovery of genes involved in flowering and the control of sex in plants.
>    
>Since it can take over a decade for some of the species we work with to reach reproductive maturity, functional verification of a potential sex-determining gene is unrealistic in these species. 
>
>Instead, we intended to use comparative genomic analyses, such as orthology inference, to identify the most likely gene candidates for sex-determination within the non-recombining regions of sex chromosomes.

I intend to create three sets of tools for this database using Python: 

1. **Interact & Extract:** this set of Python commands will allow users to query the database and extract gene IDs or functional descriptions for further analysis. (**Python Packages:** [pandas v1.5.3](https://pandas.pydata.org/docs/index.html), [pysqlite v2.8.3](https://pypi.org/project/pysqlite/))

    -   These commands would take regular expressions as an input and a column identifier, then output a list of rows where that column pattern matches. 

2. **Retrieve:** this tool will allow users to retrieve lists of gene transcripts and protein sequences from NCBI based on a list of gene IDs. (**Python Package:** [biopython v1.81](https://biopython.org/), specifically the [Bio.Entrez](https://biopython.org/docs/1.76/api/Bio.Entrez.html) package) 

    -  This tool would take a text file with a list of gene IDs as an input and output a fasta file with the desired sequences corresponding to those gene IDs, which will further serve as the input for the OrthoFinder2 run needed to complete "Parse" (see below).

3. **Parse:** two commands that will help users match their gene IDs with the appropriate orthogroups generated by Orthofinder. (**Python Libraries:** csv)

   -  The first command would take the text file from Retrieve and the output of OrthoFinder2, then return a text file corresponding to the Orthogroups associated with each gene ID in retrieve. This is important because OrthoFinder2 runs typically involve multiple species to correctly allocate the orthogroups, so the output is often extremely unwieldy and difficult to find the results for your species of interest. In the second command, the text file produced as the output in the previous step will be used as an input along with the species identifier for your target species, which will be used to find the corresponding gene IDs in your target species based on the Orthogroups identified through the first command using the initial gene IDs from FloralCoreDB. 
   
   -  **Note:** This tool requires the output of [OrthoFinder2](https://github.com/davidemms/OrthoFinder).
   -  **Next Steps:** After identifying the orthologs of these floral genes in their target species with Orthofinder2, the users will be able to harness alignment tools or even [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569839/) to create a local database of their target genome/chromosome and see if any orthologs appear within their regions of interest. 

**Data:** The main datatypes being used for this project are gene IDs and peptide sequences. I do not currently have an example of FloralCoreDB that I can share for this proposal, but I have uploaded a small dataset of *Arabidopsis* flowering time genes that provides a general idea of what I will be doing: identifying sets of functionally-verified genes from literature and then providing those gene IDs alongside descriptions from the literature about that gene and its involvement in the plant's reproductive cycles. Although I intend to develop scripts that will summarize other aspects of this database to build figures for a review paper, that is outside of the scope of this project proposal. 
