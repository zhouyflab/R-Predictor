# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/11 9:21
@Auth ： zhenyaliu
@File ：pfam_tir_rpw8.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import subprocess
import argparse
import os
import re


def parse_args():
    parser = argparse.ArgumentParser(description='Predicting coiled-coil by Paircoil2')
    parser.add_argument('--fasta1', type=str, default='./data', help='path to the fasta file')
    parser.add_argument('--fasta2', type=str, default='./data', help='path to the fasta file')
    parser.add_argument('--dir', type=str, default='./hmm', help='path to the hmm')
    return parser.parse_args()

def is_file_empty(file_path):
    return os.stat(file_path).st_size == 0

def create_file_empty(file_path):
    with open(file_path,"w") as f:
        pass

def process_pfam(path,first):
    if first == True:
        tir = {}
        rpw8 = {}
        with open(path,"r") as f:
            for i in f.readlines():
                i = i.strip()
                if i != "" and "#" not in i:
                    if re.search("CL0173",i) != None:
                        if "Family" in i:
                            target = re.findall("\d+\.?\d*", i.split("Family")[0])
                        else:
                            target = re.findall("\d+\.?\d*", i.split("Domain")[0])
                        tir[i.split(" ")[0]] = []
                        tir[i.split(" ")[0]].append(target[1]+"-"+target[2])
                        continue
                    else:
                        target = re.findall("\d+\.?\d*", i.split("Family")[0])
                        rpw8[i.split(" ")[0]] = []
                        rpw8[i.split(" ")[0]].append(target[1]+"-"+target[2])
                        continue
                else:
                    continue
        return tir,rpw8
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

def prosite(path):
    tir = {}
    rpw8 = {}
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if i.split("\t")[3] == "PS50104":
                if i.split("\t")[0] not in tir.keys():
                    tir[i.split("\t")[0]] = []
                    tir[i.split("\t")[0]].append(i.split("\t")[1]+"-"+i.split("\t")[2])
                else:
                    tir[i.split("\t")[0]].append(i.split("\t")[1]+"-"+i.split("\t")[2])
            else:
                if i.split("\t")[0] not in rpw8.keys():
                    rpw8[i.split("\t")[0]] = []
                    rpw8[i.split("\t")[0]].append(i.split("\t")[1]+"-"+i.split("\t")[2])
                else:
                    rpw8[i.split("\t")[0]].append(i.split("\t")[1]+"-"+i.split("\t")[2])
    return tir,rpw8

def paircoil2(path):
    target = {}
    tmp = []
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if "Sequence Code" in i:
                if tmp != []:
                    target[list(target.keys())[-1]].append(tmp[0] + "-" + tmp[-1])
                    tmp = []
                target[">"+i.split(": ")[1]] = []

            else:
                if "#" not in i:
                    cc = re.findall("\d+\.?\d*", i)
                    if float(cc[1]) <= 0.1:
                        tmp.append(cc[0])
                    else:
                        if tmp != [] and float(cc[1]) > 0.1:
                            target[list(target.keys())[-1]].append(tmp[0]+"-"+tmp[-1])
                            tmp = []
                else:
                    continue
        if tmp != []:
            target[list(target.keys())[-1]].append(tmp[0] + "-" + tmp[-1])
    target = {key: value for key, value in target.items() if target[key] != []}
    return target

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

def writeprotein(protein,path):
    with open(path,"w") as f:
        for k in protein.keys():
            f.write(k+"\n")
            f.write(protein[k]+"\n")

def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
            
def main(args):
    create_folder(args.dir+"/outcome")
    #nl cnl
    if is_file_empty(args.fasta1) != True:
        protein_nb_lrr_notir_norpw8 = ProteinToDict(args.fasta1)
        #here to run paircoil2
        cnl_path = " /public/home/liuzhenya/tool/paircoil2/./paircoil2 -win 21 "+args.fasta1+" "+args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8_cc.txt"
        subprocess.run(cnl_path,shell=True,check=True)
        nb_lrr_notir_norpw8_cc = paircoil2(args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8_cc.txt")
        generate_protein(protein_nb_lrr_notir_norpw8,nb_lrr_notir_norpw8_cc,args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8_cc.fasta")
        protein_nb_lrr_notir_norpw8_cc = ProteinToDict(args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8_cc.fasta")
        writeprotein(protein_nb_lrr_notir_norpw8_cc,args.dir+"/outcome/"+args.fasta1.split("_")[0]+"_cnl.fasta")
        generate_protein_nopknb(protein_nb_lrr_notir_norpw8,nb_lrr_notir_norpw8_cc,args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8_nocc.fasta")
        protein_nb_lrr_notir_norpw8_nocc = ProteinToDict(args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8_nocc.fasta")
        writeprotein(protein_nb_lrr_notir_norpw8_nocc,args.dir+"/outcome/"+args.fasta1.split("_")[0]+"_nl.fasta")
    else:
        create_file_empty(args.fasta1.split("_")[0]+"_nb_lrr_notir_norpw8.fasta")
        create_file_empty(args.dir+"/outcome/"+args.fasta1.split("_")[0]+"_nl.fasta")
        create_file_empty(args.fasta1.split("_")[0]+"_nb_lrr_cc.fasta")
        create_file_empty(args.dir+"/outcome/"+args.fasta1.split("_")[0]+"_cnl.fasta")

    #n cn
    if is_file_empty(args.fasta2) != True:
        protein_nb_nolrr_notir_norpw8 = ProteinToDict(args.fasta2)
        #here to run paircoil2
        cn_path = "/public/home/liuzhenya/tool/paircoil2/./paircoil2 -win 21 "+args.fasta2+" "+args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_cc.txt"
        subprocess.run(cn_path,shell=True,check=True)
        nb_nolrr_notir_norpw8_cc = paircoil2(args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_cc.txt")
        generate_protein(protein_nb_nolrr_notir_norpw8,nb_nolrr_notir_norpw8_cc,args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_cc.fasta")
        protein_nb_nolrr_notir_norpw8_cc = ProteinToDict(args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_cc.fasta")
        writeprotein(protein_nb_nolrr_notir_norpw8_cc,args.dir+"/outcome/"+args.fasta2.split("_")[0]+"_cn.fasta")
        generate_protein_nopknb(protein_nb_nolrr_notir_norpw8,nb_nolrr_notir_norpw8_nocc,args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_nocc.fasta")
        protein_nb_nolrr_notir_norpw8_nocc = ProteinToDict(args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_nocc.fasta")
        writeprotein(protein_nb_nolrr_notir_norpw8_nocc,args.dir+"/outcome/"+args.fasta2.split("_")[0]+"_n.fasta")
    else:
        create_file_empty(args.fasta2.split("_")[0]+"_nb_nolrr_cc.fasta")
        create_file_empty(args.dir+"/outcome/"+args.fasta2.split("_")[0]+"_cn.fasta")
        create_file_empty(args.fasta2.split("_")[0]+"_nb_nolrr_notir_norpw8_nocc.fasta")
        create_file_empty(args.dir+"/outcome/"+args.fasta2.split("_")[0]+"_n.fasta")




if __name__ == '__main__':
    args = parse_args()
    main(args)
