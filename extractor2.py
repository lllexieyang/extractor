# Copy right reserved : Lu Yang (yanglu2016@cau.edu.cn)
# Last change: Nov 8 2021
# Version 2.1

import os
from optparse import OptionParser
import collections
from Bio import SeqIO


def get_arguments():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input", action="store", dest="input", help="input fasta file", default=[])
    parser.add_option("-e", "--extract", action="store", dest="extract", help="extract contigs from file",
                      default=[])
    parser.add_option("-r", "--remove", action="store", dest="remove", help="remove contigs from file",
                      default=[])
    parser.add_option("-o", "--output", action="store", dest="output", help="output path, default is ./", default="")
    parser.add_option("-b", "--blast", action="store", help="database for blast to find target contigs",
                      default="")
    parser.add_option("-p", "--plasmid", action="store", dest="plasmid", help="provide mlplasmids result for division",
                      default="")
    parser.add_option("-x", "--prefix", action="store", dest="prefix", help="filename output prefix", default='extract')
    parser.add_option("-l", "--locustag", action="store", dest="locustag", help="locustag prefix", default='')

    return parser.parse_args()


def extract(f, ids, out, locus):
    n = 1
    ids_found_in_fasta = []
    print ('Open file ', f)
    print ('Extract contigs: ', ids)
    for record in SeqIO.parse(f, "fasta"):
        if record.id in ids:
            out_name = out[:-6] + '_' + str(n) + '.fasta'
            o = open(out_name, "w")
            ids_found_in_fasta.append(record.id)
            record.id = f[:-6] + '_' + locus + '_' + str(n)
            o.write(record.format("fasta"))
            n += 1
    o.close()
    save = 0
    for i in ids:
        if i not in ids_found_in_fasta:
            print("id \033[1;31m%s\033[0m not found in the fasta file. Please provide the FULL contig name." % i)
        else:
            save = 1
    if save:
        print('Saved in ', out)
    return 0


def remove(f, ids, out):
    o = open(out, "w")
    print('Open file ', f)
    print('Remove contigs: ', ids)
    save = 0
    for record in SeqIO.parse(f, "fasta"):
        if record.id not in ids:
            o.write(record.format("fasta"))
        else:
            save = 1
    o.close()
    if save:
        print('Contigs have been removed. New file was saved as ', out)
    else:
        os.remove(out)
        print("id \033[1;31m%s\033[0m not found in the fasta file." % ids)
    return 0


def divide(f, ids, o):
    ids_found_in_fasta = []
    out1 = o + f[:-6] + '_plasmid.fasta'
    out2 = o + f[:-6] + '_chromosome.fasta'
    o1 = open(out1, "w")
    o2 = open(out2, "w")
    print('\nOpen file ', f)
    for record in SeqIO.parse(f, "fasta"):
        if record.id in ids:
            o1.write(record.format("fasta"))
            ids_found_in_fasta.append(record.id)
        else:
            o2.write(record.format("fasta"))
    o1.close()
    o2.close()
    save1 = 0
    for i in ids:
        if i not in ids_found_in_fasta:
            print("id \033[1;31m%s\033[0m not found in the fasta file. Please provide the correct input file." % i)
        else:
            save1 = 1
    if save1:
        print(f, 'division succeed.')
        print('Chromosome is saved in ', out2)
        print('Plasmids is saved in ', out1)
    return 0


def get_blast(f,db,type):
    if not os.path.exists(db + '.nin'):
        if os.path.exists(db):
            os.system("makeblastdb -in %s -dbtype nucl" % db)
        else:
            print(db, 'not exists in current path')
    blast = f[:-6] + ".blast"
    if os.path.exists(blast):
        print('Blast results are already exist in current path, and will be used directly.')
    else:
        os.system("blastn -query " + f + " -db " + db +
                  " -outfmt 6 -out " + blast)
#                  " -outfmt '6 qseqid sacc pident qlen length sstart send' -out " + blast)

    b = open(blast, "r")
    for line in b:
        if type == 2:
            if line.split("\t")[0] not in lists[f]:
                lists[f].append(line.split("\t")[0])
        elif type == 4:
            seqs = []
            seqs.append(line.split("\t")[0])
            seqs.append(line.split("\t")[6])
            seqs.append(line.split("\t")[7])
            lists[f].append(seqs)
    return lists


def main():
    (options, args) = get_arguments()
    locus = options.locustag
    if options.input:
        if options.extract:
            lists[options.input].append(options.extract)
        elif options.remove:
            lists[options.input].append(options.remove)
        elif options.blast:
            get_blast(options.input, options.blast, 2)
#            print(get_blast(options.input, options.blast, 4))
        elif options.plasmid:
            p = open(options.plasmid, "r")
            header = 1
            for line in p:
                if header == 0:
                    i = (line.split()[3]).replace('"', '')
                    lists[options.input].append(i)
                header = 0
            for file_name in lists:
                ids = lists[file_name]
                divide(file_name, ids, options.output)
            exit()
    elif options.extract or options.remove:
        if options.extract:
            f = open(options.extract, "r")
        elif options.remove:
            f = open(options.remove, "r")
        for line in f:
            line = line.replace("\n", "\t\n")
            file_name = line.split("\t")[0]
            id_list = line.split("\t")[1]
            lists[file_name].append(id_list)
            print(lists)
    else:
        print("No input sequence file provided (-i)\n -h for help")

    for file_name in lists:
        ids = lists[file_name]
        if options.remove:
            outname = options.output + file_name[:-6] + '_clean.fasta'
            remove(file_name, ids, outname)
        else:
            outname = options.output + file_name[:-6] + '_' + options.prefix + '.fasta'
            extract(file_name, ids, outname, locus)

    return 0


lists = collections.defaultdict(list)


if __name__ == '__main__':
    main()

