# -*- coding: utf-8 -*-
"""
@Time ： 2024/3/7 14:33
@Auth ： zhenyaliu
@File ：pipeline.py
@IDE ：PyCharm
@Mail：zhenyaliu77@gmail.com

"""
import subprocess

#protein_path = "/root/autodl-tmp/Rpredictor/data_grape/protein.fasta"
protein_path = "/root/autodl-tmp/Rpredictor/test/test.fasta"
work_path = "/root/autodl-tmp/Rpredictor"
subprocess.run(f"conda run -n pfam_scan python pfam_pk_nb.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n signalp python signal_rlk_rlp.py --fasta {protein_path}",shell=True)
subprocess.run(f"conda run -n esm-lrr python esm-lrr.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n pfam_scan python pfam_lysm.py --fasta {protein_path} --dir {work_path}",shell=True)
subprocess.run(f"conda run -n pfam_scan python pfam_tir_rpw8.py --fasta {protein_path} --dir {work_path}",shell=True)
