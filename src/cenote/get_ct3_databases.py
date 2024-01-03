#!/usr/bin/env python

import argparse
import sys, os
import subprocess

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_ct3_dbs(): 
    pathname = os.path.dirname(sys.argv[0])
    cenote_script_path = os.getcwd()   
    #print(cenote_script_path) 

    parser = argparse.ArgumentParser(description='Update and/or download databases associated with \
                                    Cenote-Taker 3. HMM (hmmer) databases: updated October 8th 2023.\
                                    RefSeq Virus taxonomy DB compiled July 31, 2023.')

    required_args = parser.add_argument_group(' REQUIRED ARGUMENTS')

    required_args.add_argument("-o", dest="C_DBS", type=str, required=True, 
                            help='output directory when database will be downloaded')

    optional_args = parser.add_argument_group('Use options to pick databases to update.')

    optional_args.add_argument("--hmm", dest="HMM_DB", type=str2bool, default=False, 
                            help=' Default: False -- choose: True -or- False')
    optional_args.add_argument("--mmseqs_tax", dest="MMSEQS_TAX", type=str2bool, default=False, 
                            help=' Default: False -- choose: True -or- False')
    optional_args.add_argument("--mmseqs_cdd", dest="MMSEQS_CDD", type=str2bool, default=False, 
                            help=' Default: False -- choose: True -or- False')
    optional_args.add_argument("--domain_list", dest="DOM_LIST", type=str2bool, default=False, 
                            help=' Default: False -- choose: True -or- False')
    optional_args.add_argument("--hhCDD", dest="HHCDD", type=str2bool, default=False, 
                               help=' Default: False -- choose: True -or- False')
    optional_args.add_argument("--hhPFAM", dest="HHPFAM", type=str2bool, default=False, 
                               help=' Default: False -- choose: True -or- False')
    optional_args.add_argument("--hhPDB", dest="HHPDB", type=str2bool, default=False, 
                               help=' Default: False -- choose: True -or- False')


    args = parser.parse_args()

    if not os.path.isdir(str(args.C_DBS)):
        os.makedirs(str(args.C_DBS))

    def is_tool(name):
        """Check whether `name` is on PATH."""
        from distutils.spawn import find_executable
        return find_executable(name) is not None

    if not is_tool("mmseqs") :
        print("mmseqs is not found. Exiting.")
        quit()

    if str(args.HMM_DB) == "True":
        # https://zenodo.org/records/8436361/files/hmmscan_DBs.tgz
        print ("running HMM database update/install")
        HMM_outdir = os.path.join(str(args.C_DBS), "hmmscan_DBs")
        if not os.path.isdir(HMM_outdir):
            os.makedirs(HMM_outdir, exist_ok=True)
        subprocess.call(['wget', '--directory-prefix=' + str(HMM_outdir), 
                        'https://zenodo.org/records/8436361/files/hmmscan_DBs.tgz'])
        subprocess.call(['tar', '-xvf', os.path.join(HMM_outdir, 'hmmscan_DBs.tgz'),
                        '-C',  str(HMM_outdir)])
        subprocess.call(['rm', '-f', os.path.join(HMM_outdir, 'hmmscan_DBs.tgz')])

    if str(args.MMSEQS_TAX) == "True":
        # https://zenodo.org/records/8436361/files/refseq_virus_prot.fasta.gz
        print ("running mmseqs taxdb database update/install")
        if not is_tool("mmseqs") :
            print("mmseqs is not found. Exiting. Is conda environment activated?")
            quit()
        MMSEQS_outdir = os.path.join(str(args.C_DBS), "mmseqs_DBs")
        if not os.path.isdir(MMSEQS_outdir):
            os.makedirs(MMSEQS_outdir, exist_ok=True)
        subprocess.call(['wget', '--directory-prefix=' + str(MMSEQS_outdir), 
                        'https://zenodo.org/records/8436361/files/refseq_virus_prot.fasta.gz'])
        subprocess.call(['gunzip', '-d', os.path.join(MMSEQS_outdir, 'refseq_virus_prot.fasta.gz')])
        #subprocess.call(['rm', '-f', os.path.join(MMSEQS_outdir, 'refseq_virus_prot.fasta.gz')])
        subprocess.call(['wget', '--directory-prefix=' + str(MMSEQS_outdir), 
                        'https://zenodo.org/records/8436361/files/refseq_virus_prot_taxids.mmseqs_fmt.tsv'])
        subprocess.call(['mmseqs', 'createdb', os.path.join(MMSEQS_outdir, 'refseq_virus_prot.fasta'), 
                        os.path.join(MMSEQS_outdir, 'refseq_virus_prot_taxDB')])
        subprocess.call(['mmseqs', 'createtaxdb', os.path.join(MMSEQS_outdir, 'refseq_virus_prot_taxDB'), 
                        os.path.join(MMSEQS_outdir, 'tmp'), '--tax-mapping-file', 
                        os.path.join(MMSEQS_outdir, 'refseq_virus_prot_taxids.mmseqs_fmt.tsv')])

    if str(args.MMSEQS_CDD) == "True":
        #https://zenodo.org/records/8436361/files/cddid_all.tbl
        print ("running mmseqs CDD database update/install")
        subprocess.call(['wget', '--directory-prefix=' + str(MMSEQS_outdir), 
                        'https://zenodo.org/records/8436361/files/cddid_all.tbl'])
        if not is_tool("mmseqs") :
            print("mmseqs is not found. Exiting. Is conda environment activated?")
            quit()
        MMSEQS_outdir = os.path.join(str(args.C_DBS), "mmseqs_DBs")
        if not os.path.isdir(MMSEQS_outdir):
            os.makedirs(MMSEQS_outdir, exist_ok=True)

        subprocess.call(['mmseqs', 'databases', 'CDD', os.path.join(MMSEQS_outdir, 'CDD'), 
                        os.path.join(MMSEQS_outdir, 'tmp')])


    #https://zenodo.org/records/8436361/files/viral_cdds_and_pfams_191028.txt
    if str(args.DOM_LIST) == "True":
        subprocess.call(['wget', '--directory-prefix=' + str(args.C_DBS), 
                        'https://zenodo.org/records/8436361/files/viral_cdds_and_pfams_191028.txt'])

    if str(args.HHCDD) == "True":
        print ("running hhsuite CDD database update/install")
        subprocess.call(['rm', '-r', '-f', str(args.C_DBS) + '/NCBI_CD/'])
        isExist = os.path.exists(str(args.C_DBS) + '/NCBI_CD')
        if not isExist:
            os.makedirs(str(args.C_DBS) + '/NCBI_CD/', exist_ok=True)
        subprocess.call(['wget', '--directory-prefix=' + str(args.C_DBS), 
                         'https://zenodo.org/records/3660537/files/NCBI_CD_hhsuite.tgz'])
        subprocess.call(['tar', '-xvf', 'NCBI_CD_hhsuite.tgz'])
        subprocess.call(['rm', '-f', 'NCBI_CD_hhsuite.tgz'])

    if str(args.HHPFAM) == "True":
        print ("running PFAM database update/install")
        subprocess.call(['rm', '-r', '-f', str(args.C_DBS) + '/pfam_32_db/'])
        isExist = os.path.exists(str(args.C_DBS) + '/pfam_32_db')
        if not isExist:
            os.makedirs(str(args.C_DBS) + '/pfam_32_db/', exist_ok=True)
        subprocess.call(['wget', '--directory-prefix=' + str(args.C_DBS), 
                         'http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_32.0.tar.gz'])
        subprocess.call(['tar', '-xvf', 'pfamA_32.0.tar.gz', '--directory', str(args.C_DBS) + '/pfam_32_db'])
        subprocess.call(['rm', '-f', 'pfamA_32.0.tar.gz'])

    if str(args.HHPDB) == "True":
        print ("running PDB database update/install. This could take around 2 hours.")
        subprocess.call(['rm', '-r', '-f', str(args.C_DBS) + '/pdb70/'])
        isExist = os.path.exists(str(args.C_DBS) + '/pdb70')
        if not isExist:
            os.makedirs(str(args.C_DBS) + '/pdb70/', exist_ok=True)
        subprocess.call(['wget', '--directory-prefix=' + str(args.C_DBS), 
                         'http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz'])
        subprocess.call(['tar', '-xvf', 'pdb70_from_mmcif_latest.tar.gz', '--directory', str(args.C_DBS) + '/pdb70'])
        subprocess.call(['rm', '-f', 'pdb70_from_mmcif_latest.tar.gz'])

if __name__ == "__main__":
    get_ct3_dbs()