# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/7 15:46
@Auth ： zhenyaliu
@File ：signal_rlk_rlp.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import argparse
import re
import subprocess
import os 


def is_file_empty(file_path):
    return os.stat(file_path).st_size == 0

def create_file_empty(file_path):
    with open(file_path,"w") as f:
        pass

def parse_args():
    parser = argparse.ArgumentParser(description='Predicting signal peptide module of Rpredictor')
    parser.add_argument('--fasta', type=str, default='./data', help='path to the fasta file')
    return parser.parse_args()

def ProteinToDict(path):
    #Return a dictionary with id as key and sequence as value
    protein_seq = ''
    new_protein_seq = ''
    protein_id = ''
    protein_id_seq = {}
    fin_protein = open(path,'r')
    for i in fin_protein.readlines():
        i = i.strip('\n')
        if re.search('>',i) != None and protein_seq != '':
            for j in protein_seq:
                    new_protein_seq += j
            protein_id_seq[protein_id] = new_protein_seq
            new_protein_seq = ''
            protein_seq = ''
            protein_id = ''
        if re.search('>',i) != None:
            protein_id = i.split(" ")[0]
            continue
        protein_seq += i
    for j in protein_seq:
        new_protein_seq += j
        protein_id_seq[protein_id] = new_protein_seq
        fin_protein.close()
    return protein_id_seq

def generate_protein(protein,dic,outpath):
    with open(outpath,"w") as w:
        for key in dic.keys():
            for kkey in protein.keys():
                if re.search(key,kkey) != None:
                    w.write(kkey+"\n")
                    w.write(protein[kkey]+"\n")
                    break
                else:
                    continue

def generate_protein_nopknb(protein,dic,outpath):
    with open(outpath,"w") as w:
        for key in protein.keys():
            if key not in list(dic.keys()):
                w.write(key+"\n")
                w.write(protein[key]+"\n")
            else:
                continue

def signalp(path):
    target = {}
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if "#" in i:
                continue
            else:
                if "signal_peptide" in i:
                    s = re.findall("\\d+\\.?\\d*", i.split("signal_peptide")[1])
                    if i.split(" ")[0] not in target.keys():
                        target[i.split("\t")[0]] = []
                        target[i.split("\t")[0]].append(s[0]+"-"+s[1])
                    else:
                        target[i.split("\t")[0]].append(s[0]+"-"+s[1])
                else:
                    continue
    return target

def tmhmm(path):
    target = {}
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if "#" in i:
                continue
            else:
                if "TMhelix" in i:
                    tm = re.findall("\\d+\\.?\\d*", i.split("TMhelix")[1])
                    if i.split("\t")[0] not in target.keys():
                        target[i.split("\t")[0]] = []
                        target[i.split("\t")[0]].append(tm[0]+"-"+tm[1])
                    else:
                        target[i.split("\t")[0]].append(tm[0] + "-" + tm[1])
                else:
                    continue
    return target

def main(args):
    if is_file_empty(args.fasta.split(".")[0]+"_pk_tm.fasta") != True:
        rlk_signal = "signalp6 --fastafile "+args.fasta.split(".")[0]+"_pk_tm.fasta"+" --format txt --organism eukarya --output_dir "+args.fasta.split(".")[0]+"_pk_tm_signal"+" --mode fast"
        subprocess.run(rlk_signal, shell=True, check=True)
        rlk_s = signalp(args.fasta.split(".")[0]+"_pk_tm_signal/output.gff3")
        protein_pk_tm = ProteinToDict(args.fasta.split(".")[0]+"_pk_tm.fasta")
        generate_protein(protein_pk_tm,rlk_s,args.fasta.split(".")[0]+"_pk_tm_s.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_pk_tm_s.fasta")

    if is_file_empty(args.fasta.split(".")[0]+"_pk_notm.fasta") != True:
        pk_notm_signal = "signalp6 --fastafile "+args.fasta.split(".")[0]+"_pk_notm.fasta"+" --format txt --organism eukarya --output_dir "+args.fasta.split(".")[0]+"_pk_notm_signal"+" --mode fast"
        subprocess.run(pk_notm_signal, shell=True, check=True)
        pk_notm_s = signalp(args.fasta.split(".")[0]+"_pk_notm_signal/output.gff3")
        protein_pk_notm = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm.fasta")
        generate_protein(protein_pk_notm,pk_notm_s,args.fasta.split(".")[0]+"_pk_notm_s.fasta")
        protein_pk_notm_s = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm_s.fasta")
        generate_protein_nopknb(protein_pk_notm,protein_pk_notm_s,args.fasta.split(".")[0]+"_pk_notm_nos.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_pk_notm_s.fasta")
        protein_pk_notm = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm.fasta")
        protein_pk_notm_s = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm_s.fasta")
        generate_protein_nopknb(protein_pk_notm,protein_pk_notm_s,args.fasta.split(".")[0]+"_pk_notm_nos.fasta")

    if is_file_empty(args.fasta.split(".")[0]+"_nopknb.fasta") != True:
        rlp_signal = "signalp6 --fastafile "+args.fasta.split(".")[0]+"_nopknb.fasta"+" --format txt --organism eukarya --output_dir "+args.fasta.split(".")[0]+"_nopknb_signal"+" --mode fast"
        subprocess.run(rlp_signal, shell=True, check=True)
        rlp_s = signalp(args.fasta.split(".")[0]+"_nopknb_signal/output.gff3")
        protein_no_pknb = ProteinToDict(args.fasta.split(".")[0]+"_nopknb.fasta")
        generate_protein(protein_no_pknb,rlp_s,args.fasta.split(".")[0]+"_nopknb_s.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_nopknb_s.fasta")


    #rlk_s = signalp(args.fasta.split(".")[0]+"_pk_tm_signal/output.gff3")
    #pk_notm_s = signalp(args.fasta.split(".")[0]+"_pk_notm_signal/output.gff3")
    #rlp_s = signalp(args.fasta.split(".")[0]+"_nopknb_signal/output.gff3")
    #protein_pk_tm = ProteinToDict(args.fasta.split(".")[0]+"_pk_tm.fasta")
    #protein_pk_notm = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm.fasta")
    #protein_no_pknb = ProteinToDict(args.fasta.split(".")[0]+"_nopknb.fasta")
    #generate_protein(protein_pk_tm,rlk_s,args.fasta.split(".")[0]+"_pk_tm_s.fasta")
    #generate_protein(protein_pk_notm,pk_notm_s,args.fasta.split(".")[0]+"_pk_notm_s.fasta")
    #protein_pk_notm_s = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm_s.fasta")
    #generate_protein_nopknb(protein_pk_notm,protein_pk_notm_s,args.fasta.split(".")[0]+"_pk_notm_nos.fasta")
    #generate_protein(protein_no_pknb,rlp_s,args.fasta.split(".")[0]+"_nopknb_s.fasta")

    rlp_s_tm_env = "perl /root/tool/tmhmm-2.0c/bin/tmhmm "+args.fasta.split(".")[0]+"_nopknb_s.fasta"+" > "+args.fasta.split(".")[0]+"_nopknb_s_tm.txt"
    subprocess.run(rlp_s_tm_env,shell=True,check=True)
    rlp_s_tm = tmhmm(args.fasta.split(".")[0]+"_nopknb_s_tm.txt")
    protein_rlp_s = ProteinToDict(args.fasta.split(".")[0]+"_nopknb_s.fasta")
    generate_protein(protein_rlp_s,rlp_s_tm,args.fasta.split(".")[0]+"_nopknb_s_tm.fasta")
    protein_rlp_s_tm = ProteinToDict(args.fasta.split(".")[0]+"_nopknb_s_tm.fasta")
    generate_protein_nopknb(protein_rlp_s,protein_rlp_s_tm,args.fasta.split(".")[0]+"_nopknb_s_notm.fasta")

if __name__ == '__main__':
    args = parse_args()
    main(args)
