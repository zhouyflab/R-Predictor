# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/7 20:04
@Auth ： zhenyaliu
@File ：esm-lrr.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import random
from collections import Counter
from tqdm import tqdm

import torch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import esm
import scipy
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.svm import SVC, SVR
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression, SGDRegressor
from sklearn.pipeline import Pipeline
import pickle
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
    parser = argparse.ArgumentParser(description='ESM-LRR module of Rpredictor')
    parser.add_argument('--fasta', type=str, default='./data', help='path to the fasta file')
    parser.add_argument('--dir', type=str, default='./models', help='path to models')
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

def generateseqsegment(sequence,start,end):
    #Sequence represents the input sequence
    #Start represents the start position
    #End represents the end position
    pattern = ''
    len_seq = len(sequence)
    all_pattern = len_seq - end +1
    segment = []
    for i in range(0,all_pattern):
        for i in range(start,end):
            pattern += sequence[i]
        segment.append(pattern)
        start += 1
        end += 1
        pattern = ''
    return segment

def generatewindow(seq,start,end):
    #return[[],[],[],...]
    Window = []
    for i in seq:
        window = generateseqsegment(i,start,end)
        Window.append(window)
    return Window

def generatesegment(protein,outpath):
    fw = open(outpath, "w")
    count = 0
    for key in protein.keys():
        loc_dic = {}
        for n in range(23, 24, 1):
            b = generatewindow([protein[key]], 0, n)
            for i in range(len(b[0])):
                loc_dic[i + 1] = b[0][i]
        for ki in loc_dic.keys():
            fw.write(">" + str(count) + "|" + key.split(">")[1].split(" ")[0] + "@" + str(ki) + "|" + "\n")
            fw.write(loc_dic[ki] + "\n")
            count += 1

def predict(fasta_path,emb_path,outpath,esm_lrr_path):
    Xs = []
    Xs_id = []
    EMB_LAYER = 33
    for header, _seq in esm.data.read_fasta(fasta_path):
        fn = f'{emb_path}/{header}.pt'
        Xs_id.append(fn.split('|')[1])
        embs = torch.load(fn)
        Xs.append(embs['mean_representations'][EMB_LAYER])
    Xs = torch.stack(Xs, dim=0).numpy()

    with open(esm_lrr_path, 'rb') as file:
        classifier = pickle.load(file)

    predictions = classifier.predict(Xs)
    fw = open(outpath, "w")
    for i in range(len(predictions)):
        fw.write(Xs_id[i] + "\t" + str(predictions[i]) + "\n")

def read_predict(path,thr):
    predict = {}
    with open(path,"r") as f:
        for i in f.readlines():
            i = i.strip()
            if i.split("\t")[0].split("@")[0] not in predict.keys():
                predict[i.split("\t")[0].split("@")[0]] = {}
                if float(i.split("\t")[1]) >= thr:
                    predict[i.split("\t")[0].split("@")[0]][i.split("\t")[0].split("@")[1]] = float(i.split("\t")[1])
                else:
                    continue
            else:
                if float(i.split("\t")[1]) >= thr:
                    predict[i.split("\t")[0].split("@")[0]][i.split("\t")[0].split("@")[1]] = float(i.split("\t")[1])
                else:
                    continue
    return predict

def max_peak(keys_to_check,max_dic):
    max_value = float('-inf')
    max_key = None
    for key in keys_to_check:
        if key in max_dic and max_dic[key] > max_value:
            max_value = max_dic[key]
            max_key = key
    return max_key,max_value

def only_peak(predict):
    new_predict = {}
    for key in predict.keys():
        new_predict[key] = {}
        if predict[key] != {} and len(list(predict[key].keys())) >= 2:
            keys_to_check = []
            for n in range(len(list(predict[key].keys()))-1):
                if float(list(predict[key].keys())[n+1]) - float(list(predict[key].keys())[n]) <= 5:
                    keys_to_check.append(list(predict[key].keys())[n])
                    if n+1 == len(list(predict[key].keys()))-1:
                        keys_to_check.append(list(predict[key].keys())[n+1])
                        max_key, max_value = max_peak(keys_to_check, predict[key])
                        keys_to_check = []
                        new_predict[key][max_key] = max_value
                else:
                    if keys_to_check != []:
                        keys_to_check.append(list(predict[key].keys())[n])
                        #print(keys_to_check)
                        max_key,max_value = max_peak(keys_to_check,predict[key])
                        keys_to_check = []
                        new_predict[key][max_key] = max_value
                        if n + 1 == len(list(predict[key].keys())) - 1:
                            new_predict[key][list(predict[key].keys())[n+1]] = predict[key][list(predict[key].keys())[n+1]]
                    else:
                        new_predict[key][list(predict[key].keys())[n]] = predict[key][list(predict[key].keys())[n]]
                        if n + 1 == len(list(predict[key].keys())) - 1:
                            new_predict[key][list(predict[key].keys())[n+1]] = predict[key][list(predict[key].keys())[n+1]]
        else:
            if len(list(predict[key].keys())) == 1:
                new_predict[key][list(predict[key].keys())[0]] = predict[key][list(predict[key].keys())[0]]
            else:
                continue
    return new_predict

def DeleteOther(onecluster):          #make sure the best score infos >= 20aa
    #print(onecluster)
    new_onecluster = []
    Blacklist = []
    for i in range(len(onecluster)):
        if i <= len(onecluster) - 2:
            if int(onecluster[i].split('\t')[0]) + 20 > int(onecluster[i+1].split('\t')[0]):
                #distance < 20aa
                if float(onecluster[i].split('\t')[1]) > float(onecluster[i+1].split('\t')[1]):
                    #onecluster[i]'s score >= onecluster[i+1]'s score
                    if onecluster[i] not in new_onecluster and onecluster[i] not in Blacklist:
                        if new_onecluster != []:
                            if int(new_onecluster[-1].split('\t')[0]) + 20 > int(onecluster[i].split('\t')[0]):
                                #the distance between new_onecluster[-1] and onecluster[i] < 20aa
                                if float(new_onecluster[-1].split('\t')[1]) < float(onecluster[i].split('\t')[1]):
                                    #compare the score of both
                                    del new_onecluster[-1]
                                    new_onecluster.append(onecluster[i])
                                    #onecluster[i] replace new_onecluster[-1] successfully
                                else:
                                    continue
                                    #onecluster[i] faid to replace new_onecluster[-1]
                            else:
                                new_onecluster.append(onecluster[i])
                                #the distance between new_onecluster[-1] and onecluster[i] >= 20aa
                        else:
                            new_onecluster.append(onecluster[i])
                            #new_onecluster is [] so onecluster[i] become the new_onecluster[-1] directly
                    if onecluster[i+1] not in Blacklist:
                        Blacklist.append(onecluster[i+1])
                        #take the onecluster[i+1] to blacklist
                else:
                    if onecluster[i+1] not in new_onecluster and onecluster[i+1] not in Blacklist:
                        if new_onecluster != []:
                            if int(new_onecluster[-1].split('\t')[0]) + 20 > int(onecluster[i+1].split('\t')[0]):
                                if float(new_onecluster[-1].split('\t')[1]) < float(onecluster[i+1].split('\t')[1]):
                                    del new_onecluster[-1]
                                    new_onecluster.append(onecluster[i+1])
                                else:
                                    continue
                            else:
                                new_onecluster.append(onecluster[i+1])
                        else:
                            new_onecluster.append(onecluster[i+1])
                    if onecluster[i] not in Blacklist:
                        Blacklist.append(onecluster[i])
            else:
                if onecluster[i] not in new_onecluster and onecluster[i] not in Blacklist:
                    new_onecluster.append(onecluster[i])
                if i == len(onecluster) - 2:
                    new_onecluster.append(onecluster[i+1])
        else:
            break
    return new_onecluster

def get_len(li):
    new_li = []
    if len(li) >= 2:
        for n in range(len(li)-1):
            if int(li[n+1].split("\t")[0])-int(li[n].split("\t")[0]) <= 30:
                new_li.append(li[n].split("\t")[0]+"-"+str(int(li[n+1].split("\t")[0])-1))
            else:
                new_li.append(li[n].split("\t")[0]+"-"+str(int(li[n].split("\t")[0])+22))
        new_li.append(li[-1].split("\t")[0]+"-"+str(int(li[-1].split("\t")[0])+22))
    elif len(li) == 1:
        new_li.append(li[0].split("\t")[0]+"-"+str(int(li[0].split("\t")[0])+22))
    else:
        new_li = []
    return new_li

def write_lrr(target,protein,outpath):
    loc_target = {}
    for key in target.keys():
        onecluster = []
        for key_i in target[key].keys():
            onecluster.append(key_i + "\t" + str(target[key][key_i]))
        new_oncluster = DeleteOther(onecluster)
        new_oncluster = get_len(new_oncluster)
        loc_target[key] = new_oncluster
    with open(outpath,"w") as w:
        for key in loc_target.keys():
            if loc_target[key] != []:
                pattern = re.compile(rf".*{re.escape(key)}.*")
                kkey = [s for s in list(protein.keys()) if pattern.match(s)]
                w.write(kkey[0] + "\n")
                for i in loc_target[key]:
                    w.write(protein[kkey[0]][int(i.split("-")[0])-1:int(i.split("-")[1])]+"\t"+i+"\t"+str(target[key][i.split("-")[0]])+"\n")
            else:
                continue

def esmlrr(path):
    target = {}
    with open(path, "r") as f:
        for i in f.readlines():
            i = i.strip()
            if ">" in i:
                target[i.split(" ")[0]] = []
            else:
                target[list(target.keys())[-1]].append(i.split("\t")[1])
    return target

def writeprotein(protein,path):
    with open(path,"w") as f:
        for k in protein.keys():
            f.write(k+"\n")
            f.write(protein[k]+"\n")

def main(args):
    plm_path = args.dir+"/models/esm1v_t33_650M_UR90S_1.pt"
    esm_lrr_path = args.dir+"/models/esm_lrr.pickle"

    #nb
    if is_file_empty(args.fasta.split(".")[0]+"_nb.fasta") != True:
        protein_nb = ProteinToDict(args.fasta.split(".")[0]+"_nb.fasta")
        generatesegment(protein_nb,args.fasta.split(".")[0]+"_nb_segment.fasta")
        nb_lrr_extract_env = f"python extract.py {plm_path} "+args.fasta.split(".")[0]+"_nb_segment.fasta "+args.fasta.split(".")[0]+"_nb_segment_emb --repr_layers 33 --include mean"
        subprocess.run(nb_lrr_extract_env, shell=True, check=True)
        predict(args.fasta.split(".")[0]+"_nb_segment.fasta",args.fasta.split(".")[0]+"_nb_segment_emb",args.fasta.split(".")[0]+"_nb_lrr_score.txt",esm_lrr_path)
        nb_lrr_info = read_predict(args.fasta.split(".")[0]+"_nb_lrr_score.txt",1.2)
        nb_lrr_target = only_peak(nb_lrr_info)
        write_lrr(nb_lrr_target,protein_nb,args.fasta.split(".")[0]+"_nb_lrr.txt")
        nb_lrr = esmlrr(args.fasta.split(".")[0]+"_nb_lrr.txt")
        generate_protein(protein_nb,nb_lrr,args.fasta.split(".")[0]+"_nb_lrr.fasta")
        protein_nb_lrr = ProteinToDict(args.fasta.split(".")[0]+"_nb_lrr.fasta")
        generate_protein_nopknb(protein_nb,protein_nb_lrr,args.fasta.split(".")[0]+"_nb_nolrr.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_nb_lrr.fasta")
        create_file_empty(args.fasta.split(".")[0]+"_nb_nolrr.fasta")


    #rlk
    if is_file_empty(args.fasta.split(".")[0]+"_pk_tm_s.fasta") != True:
        protein_rlk_s = ProteinToDict(args.fasta.split(".")[0]+"_pk_tm_s.fasta")
        generatesegment(protein_rlk_s,args.fasta.split(".")[0]+"_pk_tm_s_segment.fasta")
        rlk_s_lrr_extract_env = f"python extract.py {plm_path} "+args.fasta.split(".")[0]+"_pk_tm_s_segment.fasta "+args.fasta.split(".")[0]+"_pk_tm_s_segment_emb --repr_layers 33 --include mean"
        subprocess.run(rlk_s_lrr_extract_env, shell=True, check=True)
        predict(args.fasta.split(".")[0]+"_pk_tm_s_segment.fasta",args.fasta.split(".")[0]+"_pk_tm_s_segment_emb",args.fasta.split(".")[0]+"_pk_tm_s_lrr_score.txt",esm_lrr_path)
        rlk_s_lrr_info = read_predict(args.fasta.split(".")[0]+"_pk_tm_s_lrr_score.txt",1.2)
        rlk_s_lrr_target = only_peak(rlk_s_lrr_info)
        write_lrr(rlk_s_lrr_target,protein_rlk_s,args.fasta.split(".")[0]+"_pk_tm_s_lrr.txt")
        rlk_s_lrr = esmlrr(args.fasta.split(".")[0]+"_pk_tm_s_lrr.txt")
        generate_protein(protein_rlk_s,rlk_s_lrr,args.fasta.split(".")[0]+"_pk_tm_s_lrr.fasta")
        protein_rlk_s_lrr = ProteinToDict(args.fasta.split(".")[0]+"_pk_tm_s_lrr.fasta")
        writeprotein(protein_rlk_s_lrr,args.dir+"/outcome/"+args.fasta.split(".")[0].split("/")[-1]+"_lrr_rlk.fasta")
        generate_protein_nopknb(protein_rlk_s,protein_rlk_s_lrr,args.fasta.split(".")[0]+"_pk_tm_s_nolrr.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_pk_tm_s_lrr.fasta")
        create_file_empty(args.dir+"/outcome/"+args.fasta.split(".")[0].split("/")[-1]+"_lrr_rlk.fasta")
        create_file_empty(args.fasta.split(".")[0]+"_pk_tm_s_nolrr.fasta")

    #pk_notm_nos
    if is_file_empty(args.fasta.split(".")[0]+"_pk_notm_nos.fasta") != True:
        protein_pk_notm_nos = ProteinToDict(args.fasta.split(".")[0]+"_pk_notm_nos.fasta")
        generatesegment(protein_pk_notm_nos,args.fasta.split(".")[0]+"_pk_notm_nos_segment.fasta")
        pk_notm_nos_lrr_extract_env = f"python extract.py {plm_path} "+args.fasta.split(".")[0]+"_pk_notm_nos_segment.fasta "+args.fasta.split(".")[0]+"_pk_notm_nos_segment_emb --repr_layers 33 --include mean"
        subprocess.run(pk_notm_nos_lrr_extract_env, shell=True, check=True)
        predict(args.fasta.split(".")[0]+"_pk_notm_nos_segment.fasta",args.fasta.split(".")[0]+"_pk_notm_nos_segment_emb",args.fasta.split(".")[0]+"_pk_notm_nos_lrr_score.txt",esm_lrr_path)
        pk_notm_nos_lrr_info = read_predict(args.fasta.split(".")[0]+"_pk_notm_nos_lrr_score.txt",1.2)
        pk_notm_nos_lrr_target = only_peak(pk_notm_nos_lrr_info)
        write_lrr(pk_notm_nos_lrr_target,protein_pk_notm_nos,args.fasta.split(".")[0]+"_pk_notm_nos_lrr.txt")
        pk_notm_nos_lrr = esmlrr(args.fasta.split(".")[0]+"_pk_notm_nos_lrr.txt")
        generate_protein_nopknb(protein_pk_notm_nos,pk_notm_nos_lrr,args.fasta.split(".")[0]+"_pk_notm_nos_nolrr.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_pk_notm_nos_nolrr.fasta")

    #rlp
    if is_file_empty(args.fasta.split(".")[0]+"_nopknb_s_tm.fasta") != True:
        protein_rlp_s_tm = ProteinToDict(args.fasta.split(".")[0]+"_nopknb_s_tm.fasta")
        generatesegment(protein_rlp_s_tm,args.fasta.split(".")[0]+"_nopknb_s_tm_segment.fasta")
        rlp_s_tm_lrr_extract_env = f"python extract.py {plm_path} "+args.fasta.split(".")[0]+"_nopknb_s_tm_segment.fasta "+args.fasta.split(".")[0]+"_nopknb_s_tm_segment_emb --repr_layers 33 --include mean"
        subprocess.run(rlp_s_tm_lrr_extract_env, shell=True, check=True)
        predict(args.fasta.split(".")[0]+"_nopknb_s_tm_segment.fasta",args.fasta.split(".")[0]+"_nopknb_s_tm_segment_emb",args.fasta.split(".")[0]+"_nopknb_s_tm_lrr_score.txt",esm_lrr_path)
        rlp_s_tm_lrr_info = read_predict(args.fasta.split(".")[0]+"_nopknb_s_tm_lrr_score.txt",1.2)
        rlp_s_tm_lrr_target = only_peak(rlp_s_tm_lrr_info)
        write_lrr(rlp_s_tm_lrr_target,protein_rlp_s_tm,args.fasta.split(".")[0]+"_nopknb_s_tm_lrr.txt")
        rlp_s_tm_lrr = esmlrr(args.fasta.split(".")[0]+"_nopknb_s_tm_lrr.txt")
        generate_protein(protein_rlp_s_tm,rlp_s_tm_lrr,args.fasta.split(".")[0]+"_nopknb_s_tm_lrr.fasta")
        protein_rlp_s_tm_lrr = ProteinToDict(args.fasta.split(".")[0]+"_nopknb_s_tm_lrr.fasta")
        writeprotein(protein_rlp_s_tm_lrr,args.dir+"/outcome/"+args.fasta.split(".")[0].split("/")[-1]+"_lrr_rlp.fasta")
        generate_protein_nopknb(protein_rlp_s_tm,protein_rlp_s_tm_lrr,args.fasta.split(".")[0]+"_nopknb_s_tm_nolrr.fasta")
    else:
        create_file_empty(args.fasta.split(".")[0]+"_nopknb_s_tm_lrr.fasta")
        create_file_empty(args.dir+"/outcome/"+args.fasta.split(".")[0].split("/")[-1]+"_lrr_rlp.fasta")
        create_file_empty(args.fasta.split(".")[0]+"_nopknb_s_tm_nolrr.fasta")

if __name__ == '__main__':
    args = parse_args()
    main(args)
