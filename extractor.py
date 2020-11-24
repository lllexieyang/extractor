# 1. extract/remove single contig from one fasta file (-i example.fasta -e/-r 1)
# 2. extract/remove contigs from multiple fasta files according to a list (-e/-r list.txt)
#    File of list contains fasta file names and contig names to extract or remove, separated by tab(\t).
#    Column 1: file name, Column 2: contig name
# 3. run blast to find the best hit and extract contigs. (-i example.fasta -b PlasmidFinder.fasta)
# 4. use mlplasmids output to divide genome into chromosome and plasmids. (-i example.fasta -p example.tab)

# for batch running: 1. use "-e" or "-r" with a list file, WITHOUT "-i";
#                    2. use "-i" and run it in a loop.

# Copy right reserved : Lu Yang (yanglu2016@cau.edu.cn)
# Last change: Nov 24 2020


import os
from optparse import OptionParser
import collections
from Bio import SeqIO


def get_arguments():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-l", "--list", action="store", dest="list", default="",
                      help="")
    parser.add_option("-i", "--input", action="store", dest="input", help="input fasta file", default=[])
    parser.add_option("-e", "--extract", action="store", dest="extract", help="extract single contig from file",
                      default=[])
    parser.add_option("-r", "--remove", action="store", dest="remove", help="remove contigs from file",
                      default=[])
    parser.add_option("-o", "--output", action="store", dest="output", help="output path, default is ./", default="")
    parser.add_option("-b", "--blast", action="store", help="database for blast to find target contigs",
                      default="")
    parser.add_option("-p", "--plasmid", action="store", dest="plasmid", help="provide mlplasmids result for division",
                      default="")

    return parser.parse_args()


def extract(f, ids, out):
    ids_found_in_fasta = []
    o = open(out, "w")
    print ('Open file ', f)
    print ('Extract contigs: ', ids)
    for record in SeqIO.parse(f, "fasta"):
        if record.id in ids:
            o.write(record.format("fasta"))
            ids_found_in_fasta.append(record.id)
    o.close()
    save = 0
    for i in ids:
        if i not in ids_found_in_fasta:
            print("id \033[1;31m%s\033[0m not found in the fasta file. Please provide the correct contig name." % i)
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


def main():
    (options, args) = get_arguments()
    lists = collections.defaultdict(list)
    if options.input:
        if options.extract:
            lists[options.input].append(options.extract)
        elif options.remove:
            lists[options.input].append(options.remove)
        elif options.blast:
            db = options.blast
            if not os.path.exists(db + '.nin'):
                if os.path.exists(db):
                    os.system("makeblastdb -in %s -dbtype nucl" % db)
                else:
                    print(db, 'not exists in current path')
            blast = options.input[:-6] + ".blast"
            if not os.path.exists(blast):
                os.system("blastn -query " + options.input + " -db " + db +
                          " -outfmt '6 qseqid sacc pident qlen length sstart send' -out " + blast)
            b = open(blast, "r")
            for line in b:
                if line.split("\t")[0] not in lists[options.input]:
                    lists[options.input].append(line.split("\t")[0])
            b.close()
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
            outname = options.output + file_name[:-6] + '_extract.fasta'
            extract(file_name, ids, outname)

    return 0


if __name__ == '__main__':
    main()
