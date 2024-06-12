import os
import csv
import pandas as pd

def get_vcf_names(vcf_path):
    with open(vcf_path, "r") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names

def ana(file_name):
    names = get_vcf_names(folder_dir + file_name)
    vcf = pd.read_csv(folder_dir + file_name, comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names)
    l = []
    for a in vcf:
        for i in range(len(a)):
            l.append(a['#CHROM'][i])
    s = set(l)
    result = dict()
    for i in s:
        result[i] = l.count(i)
    tem_name = file_name.split('.')
    tem = tem_name[0]
    print(tem)
    print(result)
    print()

folder_dir = "../example/sars/"
for images in os.listdir(folder_dir):
    if (images.endswith(".vcf")):
        ana(images)