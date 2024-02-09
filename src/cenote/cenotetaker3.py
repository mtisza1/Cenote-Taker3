#!/usr/bin/env python

import argparse
import sys, os
import subprocess
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path
from Bio import SeqIO
import time
from datetime import timedelta
import random
import string
import re
import logging
from distutils.spawn import find_executable

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

## taken directly from pharokka 
## (https://github.com/gbouras13/pharokka/blob/master/bin/input_commands.py)
def validate_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            print("FASTA checked.")
        else:
            sys.exit("Error: Input file is not in the FASTA format.\n")

## important random color ASCII-logo printer
def art_for_arts_sake():
    color_list = [
    '\033[31m', #red
    '\033[32m', #green
    '\033[33m', #orange
    '\033[34m', #blue
    '\033[35m', #purple
    '\033[36m', #cyan
    '\033[37m', #lightgrey
    '\033[91m', #lightred
    '\033[92m', #lightgreen
    '\033[93m', #yellow
    '\033[94m', #lightblue
    '\033[95m', #pink
    '\033[96m' #lightcyan
    ]
    ENDC = '\033[0m'

    rand_color = random.sample(color_list, k=3)


    print(f"{rand_color[0]}000000000000000000000000000000{ENDC}")
    print(f"{rand_color[0]}000000000000000000000000000000{ENDC}")
    print(f"{rand_color[0]}0000000000",f"{rand_color[1]}^^^^^^^^",f"{rand_color[0]}0000000000{ENDC}")
    print(f"{rand_color[0]}0000000",f"{rand_color[1]}^^^^^^^^^^^^^^",f"{rand_color[0]}0000000{ENDC}")
    print(f"{rand_color[0]}00000",f"{rand_color[1]}^^^^^",f"{rand_color[2]}CENOTE",
        f"{rand_color[1]}^^^^^",f"{rand_color[0]}00000{ENDC}")
    print(f"{rand_color[0]}00000",f"{rand_color[1]}^^^^^",f"{rand_color[2]}TAKER!",
        f"{rand_color[1]}^^^^^",f"{rand_color[0]}00000{ENDC}")
    print(f"{rand_color[0]}000000",f"{rand_color[1]}^^^^^^^^^^^^^^^^",f"{rand_color[0]}000000{ENDC}")
    print(f"{rand_color[0]}0000000",f"{rand_color[1]}^^^^^^^^^^^^^^",f"{rand_color[0]}0000000{ENDC}")
    print(f"{rand_color[0]}0000000000",f"{rand_color[1]}^^^^^^^^",f"{rand_color[0]}0000000000{ENDC}")
    print(f"{rand_color[0]}000000000000000000000000000000{ENDC}")
    print(f"{rand_color[0]}000000000000000000000000000000{ENDC}")

### entry point function for cenote-taker
def cenotetaker3():   
    ct_starttime = time.perf_counter()
    pathname = os.path.dirname(__file__)  
    cenote_script_path = os.path.abspath(pathname)      
    print(f"this script dir: {str(cenote_script_path)}")

    parentpath = Path(pathname).parents[1]

    __version__ = "3.2.1"

    Def_CPUs = os.cpu_count()

    Def_workdir = os.getcwd()

    parser = argparse.ArgumentParser(description='Cenote-Taker 3 is a pipeline for virus discovery \
                                    and thorough annotation of viral contigs and genomes. Visit \
                                    https://github.com/mtisza1/Cenote-Taker3 for help. \
                                    Version ' + str(__version__))

    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for Cenote-Taker 3 ')

    required_args.add_argument("-c", "--contigs", dest="original_contigs", type=str, required=True, 
                            help='Contig file with .fasta extension in fasta format. Each header must be unique \
                                before the first space character')
    required_args.add_argument("-r", "--run_title", dest="run_title", type=str, required=True, 
                            help='Name of this run. A directory of this name will be created. Must be unique from older \
                                runs or older run will be renamed. Must be less than 18 characters, using ONLY letters, \
                                numbers and underscores (_)')

    required_args.add_argument("-p", "--prune_prophage", dest="PROPHAGE", type=str2bool, required=True, 
                            help='True or False. Attempt to identify and remove flanking chromosomal regions from \
                                non-circular contigs with viral hallmarks (True is highly recommended for sequenced material \
                                not enriched for viruses. Virus-enriched samples probably should be False (you might check \
                                enrichment with ViromeQC). Also, please use False if --lin_minimum_hallmark_genes is set to 0)')


    optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for Cenote-Taker 3. See \
                                            https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#ModifiersPage for more \
                                            information on GenBank metadata fields')

    optional_args.add_argument("-t", "--cpu", dest="CPU", type=int, default=Def_CPUs, 
                            help=f"Default: {Def_CPUs} -- Example: 32 -- Number of CPUs available for Cenote-Taker 3. ")
    optional_args.add_argument('--version', action='version', version=str(__version__))
    optional_args.add_argument("-am", "--annotation_mode", dest="ANNOTATION_MODE", type=str2bool, default="False", 
                            help='Default: False -- Annotate sequences only (skip discovery). Only use if you believe \
                                each provided sequence is viral')
    optional_args.add_argument("-wd", "--working_directory", dest="c_workdir", type=str, default=Def_workdir, 
                            help=f"Default: {Def_workdir} -- Set working directory with absolute or relative path. \
                                run directory will be created within.")
    optional_args.add_argument("--template_file", dest="template_file", type=str, 
                            default=str(cenote_script_path) + '/dummy_template.sbt', 
                            help='Template file with some metadata. Real one required for GenBank submission. Takes a \
                                couple minutes to generate: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ ')
    optional_args.add_argument("--reads", nargs="+",
                                dest="READS", default="none", 
                                help='read file(s) in .fastq format. You can specify more than one separated by a space')
    optional_args.add_argument("--minimum_length_circular", dest="circ_length_cutoff", type=int, default='1000', 
                            help='Default: 1000 -- Minimum length of contigs to be checked for circularity. Bare minimun is \
                                1000 nts')
    optional_args.add_argument("--minimum_length_linear", dest="linear_length_cutoff", type=int, default='1000', 
                            help='Default: 1000 -- Minimum length of non-circualr contigs to be checked for viral \
                                hallmark genes.')
    optional_args.add_argument("-db", "--virus_domain_db", dest="HALL_LIST", type=str, choices=['virion', 'rdrp', 'dnarep'],
                                default=['virion', 'rdrp'], nargs="+",
                            help='default: virion rdrp -- Hits to which domain types should count as hallmark genes?\
                                \'virion\' database: genes encoding virion structural proteins, packaging \
                                proteins, or capsid maturation proteins (DNA and RNA genomes) with LOWEST false discovery \
                                rate. \'rdrp\' database: For RNA virus-derived RNA-dependent RNA polymerase. \
                                \'dnarep\' database: replication genes of DNA viruses. mostly useful for small \
                                DNA viruses, e.g. CRESS viruses')
    optional_args.add_argument("--lin_minimum_hallmark_genes", dest="LIN_MINIMUM_DOMAINS", type=int, default='1', 
                            help='Default: 1 -- Number of detected viral hallmark genes on a non-circular contig to be \
                                considered viral and recieve full annotation. \
                                \'2\' might be more suitable, yielding a false positive rate near 0. ')
    optional_args.add_argument("--circ_minimum_hallmark_genes", dest="CIRC_MINIMUM_DOMAINS", type=int, default='1', 
                            help='Default:1 -- Number of detected viral hallmark genes on a circular contig to be \
                                considered viral and recieve full annotation. For samples physically enriched for virus \
                                particles, \'0\' can be used, but please treat circular contigs without known viral \
                                domains cautiously. For unenriched samples, \'1\' might be more suitable. ')
    #optional_args.add_argument("--known_strains", dest="handle_knowns", type=str, default='do_not_check_knowns', help='Default: do_not_check_knowns -- do not check if putatively viral contigs are highly related to known sequences (via MEGABLAST). \'blast_knowns\': REQUIRES \'--blastn_db\' option to function correctly. ')
    #optional_args.add_argument("--blastn_db", dest="BLASTN_DB", type=str, default='none', help='Default: none -- Set a database if using \'--known_strains\' option. Specify BLAST-formatted nucleotide datase. Probably, use only GenBank \'nt\' database, \'nt viral\', or a subset therof, downloaded from ftp://ftp.ncbi.nlm.nih.gov/ Headers must be GenBank record format')
    #optional_args.add_argument("--enforce_start_codon", dest="ENFORCE_START_CODON", type=str2bool, default=False, 
    #                        help='Default: False -- For final genome maps, require ORFs to be initiated by a typical \
    #                            start codon? GenBank submissions containing ORFs without start codons can be rejected. \
    #                            However, if True,  important but incomplete genes could be culled from the final output. \
    #                            This is relevant mainly to contigs of incomplete genomes ')
    optional_args.add_argument("-hh", "--hhsuite_tool", dest="HHSUITE_TOOL", type=str, 
                               choices=['none', 'hhblits', 'hhsearch'], default='none', 
                            help=' default: none -- hhblits: query any of PDB, pfam, and CDD (depending on what is installed)\
                                to annotate ORFs escaping identification via upstream methods.\
                                hhsearch: a more sensitive tool, will \
                                query PDB, pfam, and CDD (depending on what is installed) to annotate ORFs. \
                                (WARNING: hhsearch takes much, much longer than hhblits and can extend the duration of the \
                                run many times over. Do not use on large input contig files). \'none\': forgoes \
                                annotation of ORFs with hhsuite. Fastest way to complete a run. ')

    optional_args.add_argument("--caller", dest="CALLER", type=str, choices=['prodigal-gv', 'prodigal', 'phanotate', 'adaptive'],
                            default='prodigal-gv', 
                            help=' ORF caller for viruses. default: prodigal-gv\
                                prodigal-gv: prodigal-gv only (prodigal with extra models for unusual viruses) (meta mode)\
                                prodigal: prodigal classic only (meta mode). phanotate: phanotate only \
                                Note: phanotate takes longer than prodigal, exponentially so for LONG input contigs.\
                                adaptive: will choose based on preliminary taxonomy call \
                                (phages = phanotate, others = prodigal-gv)')

    ## should be host prediction, instead
    #optional_args.add_argument("--crispr_file", dest="CRISPR_FILE", type=str, default='none', help='Tab-separated file with CRISPR hits in the following format: CONTIG_NAME HOST_NAME NUMBER_OF_MATCHES. You could use this tool: https://github.com/edzuf/CrisprOpenDB. Then reformat for Cenote-Taker 3')

    optional_args.add_argument("--isolation_source", dest="isolation_source", type=str, default='unknown', 
                            help='Default: unknown -- Describes the local geographical source of the organism from \
                                which the sequence was derived')
    optional_args.add_argument("--collection_date", dest="collection_date", type=str, default='unknown', 
                            help='Default: unknown -- Date of collection. this format: 01-Jan-2019, i.e. DD-Mmm-YYYY')
    optional_args.add_argument("--metagenome_type", dest="metagenome_type", type=str, default='unknown', 
                            help='Default: unknown -- a.k.a. metagenome_source')
    optional_args.add_argument("--srr_number", dest="srr_number", type=str, default='unknown', 
                            help='Default: unknown -- For read data on SRA, run number, usually beginning with \'SRR\' \
                                or \'ERR\' ')
    optional_args.add_argument("--srx_number", dest="srx_number", type=str, default='unknown', 
                            help='Default: unknown -- For read data on SRA, experiment number, usually beginning \
                                with \'SRX\' or \'ERX\' ')
    optional_args.add_argument("--biosample", dest="biosample", type=str, default='unknown', 
                            help='Default: unknown -- For read data on SRA, sample number, usually beginning \
                                with \'SAMN\' or \'SAMEA\' or \'SRS\' ')
    optional_args.add_argument("--bioproject", dest="bioproject", type=str, default='unknown', 
                            help='Default: unknown -- For read data on SRA, project number, usually beginning \
                                with \'PRJNA\' or \'PRJEB\' ')
    optional_args.add_argument("--assembler", dest="ASSEMBLER", type=str, default='unknown_assembler', 
                            help='Default: unknown_assembler -- Assembler used to generate contigs, if applicable. \
                                Specify version of assembler software, if possible. ')
    optional_args.add_argument("--molecule_type", dest="MOLECULE_TYPE", type=str, default='DNA', 
                            help='Default: DNA -- viable options are DNA - OR - RNA ')
    optional_args.add_argument("--data_source", dest="DATA_SOURCE", type=str, default='original', 
                            help=' default: original -- original data is not taken from other researchers\' public or \
                                private database. \'tpa_assembly\': data is taken from other researchers\' public or \
                                private database. Please be sure to specify SRA metadata.  ')

    ### I need to account for this or remove it:
    optional_args.add_argument("--filter_out_plasmids", dest="FILTER_PLASMIDS", type=str2bool, default=True, 
                            help=argparse.SUPPRESS)
    ### I need to account for this or remove it:
    #optional_args.add_argument("--scratch_directory", dest="SCRATCH_DIR", type=str, default="none", 
    #                        help='Default: none -- When running many instances of Cenote-Taker 3, it seems to run more \
    #                            quickly if you copy the hhsuite databases to a scratch space temporarily. Use this argument \
    #                            to set a scratch directory that the databases will be copied to (at least 100GB of scratch \
    #                            space are required for copying the databases)')
    optional_args.add_argument("--cenote-dbs", dest="C_DBS", type=str, default="default", 
                            help='DB path. If not set here, Cenote-Taker looks for environmental variable CENOTE_DBS. \
                                Then, if this variable is unset, DB path is assumed to be ' + str(parentpath))
    optional_args.add_argument("--hmmscan_dbs", dest="HMM_DBS", type=str, default="v3.1.1", 
                            help='HMMscan DB version. looks in cenote_db_path/hmmscan_DBs/')
    optional_args.add_argument("--wrap", dest="WRAP", type=str2bool, default="True", 
                            help='Default: True -- Wrap/rotate DTR/circular contigs so the start codon of an ORF is \
                                the first nucleotide in the contig/genome')
    optional_args.add_argument("--genbank", dest="GENBANK", type=str2bool, default="True", 
                            help='Default: True -- Make GenBank files (.gbf, .sqn, .fsa, .tbl, .cmt, etc)?')

    args = parser.parse_args()

    ## annotation mode overrides other arguments
    if str(args.ANNOTATION_MODE) == "True":
        args.LIN_MINIMUM_DOMAINS = 0
        args.CIRC_MINIMUM_DOMAINS = 0
        args.circ_length_cutoff = 1
        args.linear_length_cutoff = 1
        args.PROPHAGE = "False"

    HALL_TYPE = ' '.join(map(str,args.HALL_LIST))

    ## make out directory (rename any existing directory)

    out_directory = os.path.join(str(args.c_workdir), str(args.run_title))

    if not os.path.isdir(out_directory):
        os.makedirs(out_directory)
    else:
        randID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
        movdir = f"{out_directory}_old_{randID}"
        os.rename(out_directory, movdir)
        os.makedirs(out_directory)

    #### define logger #####
    logger = logging.getLogger("cenote_logger")
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(os.path.join(out_directory, f"{str(args.run_title)}_cenotetaker.log"))
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    #########################

    art_for_arts_sake()

    validate_fasta(args.original_contigs)

    ## joins read files together when provided with spaces in between
    if not args.READS == "none":
        READS = ' '.join(map(str,args.READS))
    else:
        READS = str(args.READS)

    ## DB path check/change
    if args.C_DBS == "default" and os.getenv('CENOTE_DBS') != None:
        args.C_DBS = os.getenv('CENOTE_DBS')
    elif args.C_DBS == "default":
        args.C_DBS = parentpath

    ## check that all DBs exist
    def check_ct3_dbs():
        ## checking db path
        if not os.path.isdir(str(args.C_DBS)):
            logger.warning(f"database directory is not found at {str(args.C_DBS)}. Exiting.")
            logger.warning("Check instructions at https://github.com/mtisza1/Cenote-Taker3 for installing databases\
                           and setting CENOTE_DBS environmental variable")
            sys.exit()
        ## checking hmm dir
        if not os.path.isdir(os.path.join(str(args.C_DBS), 'hmmscan_DBs', str(args.HMM_DBS))):
            logger.warning(f"hmm db directory version {str(args.HMM_DBS)} is not found at")
            logger.warning(f"{os.path.join(str(args.C_DBS), 'hmmscan_DBs', str(args.HMM_DBS))}")
            logger.warning(f"looking for others in {os.path.join(str(args.C_DBS), 'hmmscan_DBs')}")
            for path, subdirs, files in os.walk(os.path.join(str(args.C_DBS), 'hmmscan_DBs')):
                for name in files:
                    if name == 'Virion_HMMs.h3m':
                        bottom_dir = path.split('/')[-1:][0]
                        sec_dir = '/'.join(path.split('/')[0:-1])
                        if sec_dir == os.path.join(str(args.C_DBS), 'hmmscan_DBs'):
                            args.HMM_DBS = bottom_dir
                        else:
                            logger.warning("Check instructions at https://github.com/mtisza1/Cenote-Taker3\
                                           for installing databases\
                                           and setting CENOTE_DBS environmental variable")
                            logger.warning("Exiting.")
                            sys.exit()
        ## checking hmm files
        hmm_flist = ['Virion_HMMs.h3m', 'DNA_rep_HMMs.h3m', 'RDRP_HMMs.h3m', 
                     'Useful_Annotation_HMMs.h3m', 'phrogs_for_ct.h3m']
        for hf in hmm_flist:
            if not os.path.isfile(os.path.join(str(args.C_DBS), 'hmmscan_DBs', str(args.HMM_DBS), hf)):
                logger.warning(f"hmm db file is not found at")
                logger.warning(f"{os.path.join(str(args.C_DBS), 'hmmscan_DBs', str(args.HMM_DBS), hf)}")
                logger.warning("Check instructions at https://github.com/mtisza1/Cenote-Taker3 for installing databases\
                            and setting CENOTE_DBS environmental variable")
                logger.warning("Exiting.")
                sys.exit()
        ## checking mmseqs tax db
        if not os.path.isfile(os.path.join(str(args.C_DBS), 'mmseqs_DBs', 'refseq_virus_prot_taxDB')):
            logger.warning(f"mmseqs tax db file is not found at")
            logger.warning(f"{os.path.join(str(args.C_DBS), 'mmseqs_DBs', 'refseq_virus_prot_taxDB')}")
            logger.warning("Check instructions at https://github.com/mtisza1/Cenote-Taker3 for installing databases\
                           and setting CENOTE_DBS environmental variable")
            logger.warning("Exiting.")
            sys.exit()
        ## checking mmseqs cdd db
        if not os.path.isfile(os.path.join(str(args.C_DBS), 'mmseqs_DBs', 'CDD')):
            logger.warning(f"mmseqs CDD db file is not found at")
            logger.warning(f"{os.path.join(str(args.C_DBS), 'mmseqs_DBs', 'CDD')}")
            logger.warning("Check instructions at https://github.com/mtisza1/Cenote-Taker3 for installing databases\
                           and setting CENOTE_DBS environmental variable")
            logger.warning("Exiting.")
            sys.exit()
         ## checking virus domain list file
        if not os.path.isfile(os.path.join(str(args.C_DBS), 'viral_cdds_and_pfams_191028.txt')):
            logger.warning(f"virus domain list file is not found at")
            logger.warning(f"{os.path.join(str(args.C_DBS), 'viral_cdds_and_pfams_191028.txt')}")
            logger.warning("Check instructions at https://github.com/mtisza1/Cenote-Taker3 for installing databases\
                           and setting CENOTE_DBS environmental variable")
            logger.warning("Exiting.")
            sys.exit()

    check_ct3_dbs()

    ### check dependencies
    def is_tool(name):
        """Check whether `name` is on PATH."""
        return find_executable(name) is not None

    tool_dep_list = ['samtools', 'minimap2', 'tRNAscan-SE', 'seqkit', 'hhblits', 
                     'bedtools', 'phanotate.py', 'mmseqs']
    
    for tool in tool_dep_list:
        if not is_tool(tool):
            logger.warning(f"{tool} is not found. Exiting.")
            sys.exit()   


    reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
    installed_packages = [r.decode().split('==')[0] for r in reqs.split()]

    python_dep_list = ['pyhmmer', 'numpy', 'pandas', 'biopython', 'pyrodigal-gv']

    for pydep in python_dep_list:
        if pydep not in installed_packages:
            logger.warning(f"{pydep} not found in installed python packages. Exiting.")
            sys.exit() 

    ## check run_title suitability
    if re.search(r'^[a-zA-Z0-9_]+$', str(args.run_title)) and \
        len(str(args.run_title)) <= 18:
        logger.info(str(args.run_title))
    else:
        logger.warning(f"{str(args.run_title)} is not a valid name for the run title ( -r argument)")
        logger.warning( "the run title needs to be only letters, numbers and underscores (_) and \
              18 characters or less. Exiting.")
        sys.exit()


    #### define logging of subprocess (cenote_main.sh) ####
    def log_subprocess_output(pipe):
        for line in iter(pipe.readline, b''): # b'\n'-separated lines
            logger.info(line.decode("utf-8").rstrip('\n'))

    ### run the main script
    process = Popen(['bash', str(cenote_script_path) + '/cenote_main.sh', str(cenote_script_path), 
                    str(args.original_contigs), str(args.run_title), str(args.PROPHAGE), str(args.CPU),  
                    str(__version__), str(args.ANNOTATION_MODE), str(out_directory), 
                    str(args.template_file), str(READS), str(args.circ_length_cutoff), 
                    str(args.linear_length_cutoff), str(args.CIRC_MINIMUM_DOMAINS), 
                    str(args.LIN_MINIMUM_DOMAINS), str(HALL_TYPE), str(args.C_DBS), str(args.HMM_DBS), 
                    str(args.WRAP), str(args.CALLER), str(args.HHSUITE_TOOL), 
                    str(args.isolation_source), str(args.collection_date), str(args.metagenome_type), 
                    str(args.srr_number), str(args.srx_number), str(args.biosample), 
                    str(args.bioproject), str(args.ASSEMBLER), str(args.MOLECULE_TYPE), 
                    str(args.DATA_SOURCE), str(args.GENBANK)],
                    stdout=PIPE, stderr=STDOUT)

    with process.stdout:
        log_subprocess_output(process.stdout)
    exitcode = process.wait()

    ct_endtime = time.perf_counter()

    time_taken = ct_endtime - ct_starttime

    time_taken = round(time_taken, 2) 

    logger.info("This Cenote-Taker run finished in " + str(timedelta(seconds=time_taken)))

if __name__ == "__main__":
    cenotetaker3()