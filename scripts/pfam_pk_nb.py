# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/6 11:01
@Auth ： zhenyaliu
@File ：pipeline.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import subprocess
import argparse
import os
import re
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='Predicting Pkinase and NB-ARC module of Rpredictor')
    parser.add_argument('--fasta', type=str, default='./data', help='path to the fasta file')
    parser.add_argument('--type', '-t', type=str,
                        choices=['RLK', 'RLP', 'NLR', 'all'], default='all',
                        help='Types of disease resistance gene (default: all)')
    parser.add_argument('--pfam_envir', type=str, default='Pfam_scan', help='virtual environment name of Pfam_scan')
    parser.add_argument('--esmlrr_envir', type=str, default='ESM-LRR', help='virtual environment name of ESM-LRR')
    parser.add_argument('--tm_envir', type=str, default='TM', help='virtual environment name of TMHMM-2.0')
    parser.add_argument('--signal_envir', type=str, default='Signal', help='virtual environment name of SignalP6.0')
    parser.add_argument('--cc_envir', type=str, default='CC', help='virtual environment name of CoCoPRED')
    parser.add_argument('--dir', type=str, default='./hmm', help='path to the hmm')
    return parser.parse_args()

def clean_fasta(fasta_path):
    records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        record.seq = record.seq.replace("*", "")
        records.append(record)
    with open(fasta_path, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")

def create_file_empty(file_path):
    with open(file_path,"w") as f:
        pass

def process_pfam(path,first):
    if first == True:
        pk = {}
        nbarc = {}
        with open(path,"r") as f:
            for i in f.readlines():
                i = i.strip()
                if i != "" and "#" not in i:
                    if re.search("CL0016",i) != None:
                        target = re.findall(r"\d+\.?\d*", i.split("Domain")[0])
                        if i.split(" ")[0] not in pk.keys():
                            pk[i.split(" ")[0]] = []
                            pk[i.split(" ")[0]].append(target[1]+"-"+target[2])
                        else:
                            pk[i.split(" ")[0]].append(target[1] + "-" + target[2])
                    if re.search("CL0023",i) != None:
                        target = re.findall(r"\d+\.?\d*", i.split("Domain")[0])
                        if i.split(" ")[0] not in nbarc.keys():
                            nbarc[i.split(" ")[0]] = []
                            nbarc[i.split(" ")[0]].append(target[1]+"-"+target[2])
                        else:
                            nbarc[i.split(" ")[0]].append(target[1] + "-" + target[2])
                else:
                    continue
        return pk,nbarc
    else:
        dic = {}
        with open(path,"r") as f:
            for i in f.readlines():
                i = i.strip()
                if i != "" and "#" not in i:
                    target = re.findall(r"\d+\.?\d*", i.split("Domain")[0])
                    dic[i.split(" ")[0]] = []
                    dic[i.split(" ")[0]].append(target[1]+"-"+target[2])
        return dic

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

def filterprotein(protein):
    for k in protein.keys():
        if protein[k][-1] == "*":
            protein[k] = protein[k][:-1]
        else:
            continue
    return protein

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

def generate_tkp(pk,protein,outpath):
    id_list = []
    for k in pk.keys():
        if 2 <= len(pk[k]) <= 5:
            id_list.append(k)
        else:
            continue
    if len(id_list) == 0:
        create_file_empty(outpath)
    else:
        fw = open(outpath,"w")
        for i in id_list:
            fw.write(">"+i+"\n")
            fw.write(protein[">"+i]+"\n")
        fw.close()

def tmhmm(path):
    target = {}
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if "#" in i:
                continue
            else:
                if "TMhelix" in i:
                    tm = re.findall(r"\d+\.?\d*", i.split("TMhelix")[1])
                    if i.split("\t")[0] not in target.keys():
                        target[i.split("\t")[0]] = []
                        target[i.split("\t")[0]].append(tm[0]+"-"+tm[1])
                    else:
                        target[i.split("\t")[0]].append(tm[0] + "-" + tm[1])
                else:
                    continue
    return target

def difference(dic1,dic2):
    #the difference set of two dictionaries
    dic = {}
    for k in dic1.keys():
        if k not in dic2.keys():
            dic[k]=dic1[k]
        else:
            continue
    return dic

def main(args):
    os.makedirs(args.dir+"/outcome",exist_ok=True)
    clean_fasta(args.fasta)
    new_dir = args.dir+"/tmp"
    os.makedirs(new_dir,exist_ok=True)
    pkinase = "pfam_scan.pl -fasta "+args.fasta+" -dir "+args.dir+"/hmm/PK_NB_HMM -outfile "+new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk_nb.txt"
    subprocess.run(pkinase,shell=True,check=True)
    protein = filterprotein(ProteinToDict(args.fasta))
    pk,nb = process_pfam(new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk_nb.txt",True)
    generate_tkp(pk,protein,args.dir+"/outcome"+"/"+args.fasta.split("/")[-1].split(".")[0]+"_tkp.fasta")
    generate_protein(protein,pk,new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk.fasta")
    generate_protein(protein,nb,new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_nb.fasta")
    pk_protein = ProteinToDict(new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk.fasta")
    nb_protein = ProteinToDict(new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_nb.fasta")
    pknb = {**pk_protein,**nb_protein}
    generate_protein_nopknb(protein,pknb,new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_nopknb.fasta")
    
    #pk_tm
    pk_tm_env = "perl /public/home/liuzhenya001/project/r-predictor/tool/tmhmm-2.0c/bin/tmhmm "+new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk.fasta"+" > "+new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk_tm.txt"
    subprocess.run(pk_tm_env,shell=True,check=True)
    pk_tm = tmhmm(new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk_tm.txt")
    generate_protein(pk_protein,pk_tm,new_dir+"/"+args.fasta.split(".")[0].split("/")[-1]+"_pk_tm.fasta")

if __name__ == '__main__':
    args = parse_args()
    main(args)
