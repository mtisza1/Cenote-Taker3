# Cenote-Taker3

Discover and annotate the virome

`Cenote-Taker 3` is a bioinformatics tool that scales from individual genomes sequences to massive metagenome assemblies to:

1)  Identify sequences containing genes specific to viruses (virus hallmark genes)

2)  Annotate virus sequences including:

---a) adaptive ORF calling

---b) a large catalog of HMMs from virus gene families

---c) Hierarchical taxonomy assignment based on hallmark genes

---d) mmseqs2-based CDD database search

---e) tabular (.tsv) and interactive genome map (.gbf) outputs

**Also, `Cenote-Taker 3` is very fast, many many times faster than `Cenote-Taker 2`, and faster than comparable annotation using `pharokka` (in my hands)**

## Installation Instructions

*A bioconda package is forthcoming, but this is considered a beta build, so I'm holding off on that*

**This should work on MacOS and Linux**

1)  Clone this GitHub repo

2)  Using `mamba` (package manager within `conda`) and the provided yaml file, make the environment:

`mamba env create -f Cenote-Taker3/environment/ct3_beta_env.yaml`

*Versions used in test installations*

mamba 1.5.1

conda 23.7.4


3)  Activate the conda environment.

`conda activate ct3_beta`

4)  Change to repo and `pip` install command line tool.

`cd Cenote-Taker3`

`pip install .`

*You should be able to type `cenotetaker3` in terminal to bring up help menu now*

5)  Change to a directory where you'd like to install databases and run database script, specify DB directory with `-o`.

*Total DB file size of 3.0 GB after file decompression*

`cd ..`

`python Cenote-Taker3/src/cenote/get_ct3_databases.py -o ct3_DBs --hmm T --mmseqs_tax T --mmseqs_cdd T --domain_list T`

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

### Discover and Annotate, Force `prodigal` (faster)

```
cenotetaker3 -c my_metagenome_contigs.fna -r my_meta_ct3pr -p T --caller prodigal
```

### Just Annotate

```
cenotetaker3 -c my_virus_contigs.fna -r my_virs_ct3 -p F -am T
```

### Notes

`Cenote-Taker 3` is under active development, so please open an issue if anything seems unusual or any errors occur. It's likely that I've not tested every parameter combination, and bugs will be a simple fix.


### To-do list

* Add module to use `HHsearch` for gene annotation

* Add RDRP database as option for virus discovery

* Nucleotide or kmer-based species-level taxonomy

* Incorporate `prodigal-gv` as alternative ORF caller

* Update HMM database to increase consistency of names/functions between similar models