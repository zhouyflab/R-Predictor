# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/9 19:21
@Auth ： zhenyaliu
@File ：pfam_lysm.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import subprocess
import argparse
import os
import re
import os


def is_file_empty(file_path):
    return os.stat(file_path).st_size == 0

def create_file_empty(file_path):
    with open(file_path,"w") as f:
        pass

def parse_args():
    parser = argparse.ArgumentParser(description='Predicting Lysm domain module of Rpredictor')
    parser.add_argument('--fasta', type=str, default='./data', help='path to the fasta file')
    return parser.parse_args()

def process_pfam(path,first):
    if first == True:
        pk = {}
        nbarc = {}
        with open(path,"r") as f:
            for i in f.readlines():
                i = i.strip()
                if i != "" and "#" not in i:
                    if re.search("CL0016",i) != None:
                        target = re.findall("\d+\.?\d*", i.split("Domain")[0])
                        pk[i.split(" ")[0]] = []
                        pk[i.split(" ")[0]].append(target[1]+"-"+target[2])
                        continue
                    if re.search("CL0023",i) != None:
                        target = re.findall("\d+\.?\d*", i.split("Domain")[0])
                        nbarc[i.split(" ")[0]] = []
                        nbarc[i.split(" ")[0]].append(target[1]+"-"+target[2])
                        continue
                else:
                    continue
        return pk,nbarc
    else:
        dic = {}
        with open(path,"r") as f:
            for i in f.readlines():
                i = i.strip()
                if i != "" and "#" not in i:
                    target = re.findall("\d+\.?\d*", i.split("Domain")[0])
                    dic[i.split(" ")[0]] = []
                    dic[i.split(" ")[0]].append(target[1]+"-"+target[2])
        return dic

def pfam_other(path):
    dic = {}
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if i != "" and "#" not in i:
                if i.split(" ")[0] not in dic.keys():
                    dic[i.split(" ")[0]] = []
                    dic[i.split(" ")[0]].append(i.split(" ")[-1])
                else:
                    dic[i.split(" ")[0]].append(i.split(" ")[-1])
    pk_only = {}
    for k in dic.keys():
        if len(dic[k]) == 1 and dic[k][0] == "CL0016":
            pk_only[k] = dic[k][0]
    return pk_only

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

def main(args):
    #lysm-rlk and s+tm+pk
    if is_file_empty(args.fasta.split(".")[0]+"_pk_tm_s_nolrr.fasta") != True:
        rlk_lysm_path = "pfam_scan.pl -fasta "+args.fasta.split(".")[0]+"_pk_tm_s_nolrr.fasta"+" -dir /root/autodl-tmp/pfam/LysM_HMM -outfile "+args.fasta.split(".")[0]+"_pk_tm_s_nolrr_lysm.txt"
        subprocess.run(rlk_lysm_path,shell=True,check=True)
        protein_rlk_s_nolrr = ProteinToDict(args.fasta.split(".")[0]+"_pk_tm_s_nolrr.fasta")
        rlk_lysm = process_pfam(args.fasta.split(".")[0]+"_pk_tm_s_nolrr_lysm.txt",False)
        generate_protein(protein_rlk_s_nolrr,rlk_lysm,args.fasta.split(".")[0]+"_pk_tm_s_nolrr_lysm.fasta")
        pk_tm_s_nolrr_path = "pfam_scan.pl -fasta "+args.fasta.split(".")[0]+"_pk_tm_s_nolrr.fasta"+" -dir /root/autodl-tmp/pfam/PfamA -outfile "+args.fasta.split(".")[0]+"_pk_tm_s_nolrr_other.txt"
        subprocess.run(pk_tm_s_nolrr_path,shell=True,check=True)
        pk_tm_s_only = pfam_other(args.fasta.split(".")[0]+"_pk_tm_s_nolrr_other.txt")
        generate_protein(protein_rlk_s_nolrr,pk_tm_s_only,args.fasta.split(".")[0]+"_pk_tm_s_only.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_pk_tm_s_nolrr_lysm.fasta")
        create_file_empty(args.fasta.split(".")[0]+"_pk_tm_s_only.fasta")

    #pk
    if is_file_empty(args.fasta.split(".")[0]+"_pk_notm_nos_nolrr.fasta") != True:
        pk_notm_nos_nolrr_path = "pfam_scan.pl -fasta "+args.fasta.split(".")[0]+"_pk_notm_nos_nolrr.fasta"+" -dir /root/autodl-tmp/pfam/PfamA -outfile "+args.fasta.split(".")[0]+"_pk_notm_nos_nolrr_other.txt"
        subprocess.run(pk_notm_nos_nolrr_path,shell=True,check=True)
        protein_pk_notm_nos_nolrr = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm_nos_nolrr.fasta")
        pk_only = pfam_other(args.fasta.split(".")[0]+"_pk_notm_nos_nolrr_other.txt")
        generate_protein(protein_pk_notm_nos_nolrr,pk_only,args.fasta.split(".")[0]+"_pk_only.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_pk_only.fasta")

    #lysm-rlp
    if is_file_empty(args.fasta.split(".")[0]+"_nopknb_s_tm_nolrr.fasta") != True:
        rlp_lysm_path = "pfam_scan.pl -fasta "+args.fasta.split(".")[0]+"_nopknb_s_tm_nolrr.fasta"+" -dir /root/autodl-tmp/pfam/LysM_HMM -outfile "+args.fasta.split(".")[0]+"_nopknb_s_tm_nolrr_lysm.txt"
        subprocess.run(rlp_lysm_path, shell=True, check=True)
        protein_rlp_s_tm_nolrr = ProteinToDict(args.fasta.split(".")[0] + "_nopknb_s_tm_nolrr.fasta")
        rlp_lysm = process_pfam(args.fasta.split(".")[0] + "_nopknb_s_tm_nolrr_lysm.txt", False)
        generate_protein(protein_rlp_s_tm_nolrr, rlp_lysm, args.fasta.split(".")[0] + "_nopknb_s_tm_nolrr_lysm.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0] + "_nopknb_s_tm_nolrr_lysm.fasta")

    #s+lysm
    if is_file_empty(args.fasta.split(".")[0]+"_nopknb_s_notm.fasta") != True:
        lysm_path = "pfam_scan.pl -fasta "+args.fasta.split(".")[0]+"_nopknb_s_notm.fasta"+" -dir /root/autodl-tmp/pfam/LysM_HMM -outfile "+args.fasta.split(".")[0]+"_nopknb_s_notm_lysm.txt"
        subprocess.run(lysm_path, shell=True, check=True)
        protein_rlp_s_notm = ProteinToDict(args.fasta.split(".")[0] + "_nopknb_s_notm.fasta")
        lysm = process_pfam(args.fasta.split(".")[0] + "_nopknb_s_notm_lysm.txt", False)
        generate_protein(protein_rlp_s_notm, lysm, args.fasta.split(".")[0] + "_nopknb_s_notm_lysm.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0] + "_nopknb_s_notm_lysm.fasta")

if __name__ == '__main__':
    args = parse_args()
    main(args)
