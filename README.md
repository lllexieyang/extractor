# node-extract
根据blast结果提取fasta文件中含指定基因的contig

#!/usr/bin/env python

import re
import os
import sys, getopt


def Inc_result(isolate, gene):
    file = blast_input + isolate + '.blast'
    f = open(file, "r")
    line = f.readline()
    node_num = 'apple'
    feedback = output + 'node_info.txt'
    while line:
        find_gene = line.find(gene)
        if find_gene >= 0:
            r = open(feedback, "a")
            node_num = re.findall(r"NODE_(.+?)_", line)[0]
            r.write(isolate + "\t" + line + "\n")
            r.close()  
        line = f.readline()
    f.close()
    return node_num

def node_extract(isolate, node, gene):
    input_name = fasta_input + isolate + '.fasta' 
    start_point = "NODE_" + node
    stop_point = ">NODE_" + str(int(node) + 1)
    output_name = output + isolate + '_' + gene + '_' + start_point + '.fasta'
    i = open(input_name, "r")
    line = i.readline()
    while line:
        if line.find(start_point) >= 0:
            o = open(output_name, "a")
            while line:
                if line.find(stop_point) >= 0:
                    o.close()
                    break
                o.write(line)
                line = i.readline()
        line = i.readline()
    i.close()
    return 0

#options
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hf:b:o:")
except getopt.GetoptError:
    print 'Error: please try node_extract3.0.py -h for help'   
    sys.exit(2)

fasta_input = ""
blast_input = ""
output = ""


for op, value in opts:
    if op == "-f":
        fasta_input = value
    elif op == "-b":
        blast_input = value
    elif op == "-o":
        output = value
    elif op == "-h":
        print '\n python node_extract3.0.py -b <blastpath> -f <fastapath> -o <outputpath>'
        print ' for example:'
        print ' python node_extract3.0.py -b /disk1/cau/cvmylx/blast/ -f /disk1/cau/cvmylx/contigs/ -o /disk1/cau/
cvmylx/nodes/' 
        print '\n'
        sys.exit()

#get file id
path = fasta_input
file_name = []
files = os.listdir(path)
for file in files:
    parts = file.split('.')
    id = parts[0]
    file_name.append(id)

#get genes
input_genes =raw_input("please input your genes:")
genes = input_genes.split()

print file_name
print genes
for name in file_name:
    for gene in genes:
        node_num = Inc_result(name, gene)
        if node_num != 'apple':
            node_extract(name,node_num, gene)
