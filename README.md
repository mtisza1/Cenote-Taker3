# Cenote-Taker3

Discover and annotate the virome

![Logo](images/cenote-taker_3_logo.png)

Works on your laptop or HPC (compatible with MacOS and Linux)

`Cenote-Taker 3` is a virus bioinformatics tool that scales from individual genomes sequences to massive metagenome assemblies to:

1)  Identify sequences containing genes specific to viruses (virus hallmark genes)

2)  Annotate virus sequences including:

---a) adaptive ORF calling

---b) a large catalog of HMMs from virus gene families for functional annotation

---c) Hierarchical taxonomy assignment based on hallmark genes

---d) mmseqs2-based CDD database search

---e) tabular (.tsv) and interactive genome map (.gbf) outputs

**Also, `Cenote-Taker 3` is very fast, many many times faster than `Cenote-Taker 2`, and faster than comparable annotation using `pharokka` (in my hands)**

Image of example genome map:

![Map](images/genome_map1.png | width=100)

## Schematic

![Schematic](images/cenote-taker_3_schematic.png)

## Installation Instructions

*A bioconda package is forthcoming, but this is considered a beta build, so I'm holding off on that*

**This should work on MacOS and Linux**

1)  Clone this GitHub repo

2)  Using `mamba` (package manager within `conda`) and the provided yaml file, make the environment:

`mamba env create -f Cenote-Taker3/environment/ct3_env.yaml`

*Versions used in test installations*

mamba 1.5.1

conda 23.7.4

3)  Activate the conda environment.

`conda activate ct3_env`

4)  Change to repo and `pip` install command line tool.

`cd Cenote-Taker3`

`pip install .`

*You should be able to type `cenotetaker3` and 'get_ct3_dbs' in terminal to bring up help menu now*

5)  Change to a directory where you'd like to install databases and run database script, specify DB directory with `-o`.

*Total DB file size of 3.0 GB after file decompression*

`cd ..`

`get_ct3_dbs -o ct3_DBs --hmm T --mmseqs_tax T --mmseqs_cdd T --domain_list T`

<details>

  <summary>With optional hhsuite databases</summary>
  
  Warning: due to inconsistent server speed, these downloads may take over 2 hours.
  
  You may download one or more hhsuite DB.
  
  The data footprint is:
  
  | Database | Size   |
  |----------|--------|
  | CDD      | 6.1 GB |
  | pfam     | 4.6 GB |
  | pdb70    | 56 GB  |
  
  ```         
  get_ct3_dbs -o ct3_DBs --hmm T --mmseqs_tax T --mmseqs_cdd T --domain_list T --hhCDD T --hhPFAM T --hhPDB T
  ```

</details>

6)  Set the database directory as a conda environmental variable.

`conda env config vars set CENOTE_DBS=/path/to/ct3_DBs`

## Running Cenote-Taker 3

*Make sure conda environment is activated*

### Help Menu

```         
cenotetaker3 -h
```

### Test contigs

```         
cenotetaker3 -c Cenote-Taker3/test_data/testcontigs_DNA_ct2.fasta -r test_ct3 -p T
```

### Default Discover and Annotate

```         
cenotetaker3 -c my_metagenome_contigs.fna -r my_meta_ct3 -p T
```

### Discover and Annotate, Force `prodigal` (prodigal-gv is default)

```         
cenotetaker3 -c my_metagenome_contigs.fna -r my_meta_ct3pr -p T --caller prodigal
```

### Just Annotate

```         
cenotetaker3 -c my_virus_contigs.fna -r my_virs_ct3 -p F -am T
```

### Choose which HMM DBs are hallmark (virion rdrp is default)

```         
cenotetaker3 -c my_metagenome_contigs.fna -r my_meta_ct3 -p T -db virion rdrp dnarep
```

### Calculate coverage level with reads

```         
cenotetaker3 -c my_metagenome_contigs.fna -r my_meta_ct3 -p T --reads my_reads/*fastq
```

## Output Files

<pre>
{run_title}/
|   {run_title}_virus_summary.tsv                 <b><- main summary file for each virus</b>
|   {run_title}_virus_sequences.fna               <b><- all virus genome seqs</b>
|   {run_title}_virus_AA.faa                      <b><- all virus AA seqs</b>
|   {run_title}_prune_summary.tsv                 <b><- summary of pruning of each sequence</b>
|   final_genes_to_contigs_annotation_summary.tsv <b><- annotation info, all genes</b>
|   run_arguments.txt                             <b><- arguments used in this run</b>
│   {run_title}_cenotetaker.log                   <b><- main log file</b>
│
└───sequin_and_genome_maps/
│   │   {run_title}*gbf                           <b><- genome maps</b>
│   │   {run_title}*fsa                           <b><- genome sequence</b>
│   │   {run_title}*gtf                           <b><- feature table gtf format</b>
│   │   {run_title}*tbl                           <b><- feature table sequin format</b>
│   │   {run_title}*sqn                           <b><- non-human-readable sequin file for GenBank sub</b>
│   │   {run_title}*cmt                           <b><- sequin comment file</b>
│
└───ct_processing/
│   │   <b>--- many intermediate files ---</b>
</pre>

### Notes

`Cenote-Taker 3` is under active development, so please open an issue if anything seems unusual or any errors occur. It's likely that I've not tested every parameter combination, and bugs will be a simple fix.

### To-do list

-   instructions for manual curation -\> GenBank deposit of `Cenote-Taker 3` output
