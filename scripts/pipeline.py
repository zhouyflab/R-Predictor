# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/7 14:33
@Auth ： zhenyaliu
@File ：pipeline.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import subprocess
import os


# a single protein fasta file
protein_path = "/public/home/liuzhenya001/project/r-predictor/test.fasta"
work_path = "/public/home/liuzhenya001/project/r-predictor/work"
subprocess.run(f"conda run -n pfam_scan python pfam_pk_nb.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n signalp python signal_rlk_rlp.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n esm-lrr python esm-lrr.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n pfam_scan python pfam_lysm.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n pfam_scan python pfam_tir_rpw8.py --fasta {protein_path} --dir {work_path}",shell=True)

# multiple protein fasta files
#def list_all_file_paths(folder_path):
#     return [
#         os.path.abspath(os.path.join(folder_path, f))
#         for f in os.listdir(folder_path)
#         if os.path.isfile(os.path.join(folder_path, f))
#     ]

#protein_paths = list_all_file_paths("/public/home/zhangwei1/R-Predictor/data_banana/pep_noclean")
#work_path = "/public/home/zhangwei1/R-Predictor/work"
#for protein_path in protein_paths:
#    subprocess.run(f"conda run -n pfam_scan python pfam_pk_nbzw.py --fasta {protein_path} --dir {work_path}",shell=True)
#   subprocess.run(f"conda run -n signalp python signal_rlk_rlp.py --fasta {protein_path}",shell=True)
#   subprocess.run(f"conda run -n esm-lrr python esm-lrr.py --fasta {protein_path} --dir {work_path}",shell=True)
#   subprocess.run(f"conda run -n pfam_scan python pfam_lysm.py --fasta {protein_path} --dir {work_path}",shell=True)
#   subprocess.run(f"conda run -n pfam_scan python pfam_tir_rpw8.py --fasta {protein_path} --dir {work_path}",shell=True)
