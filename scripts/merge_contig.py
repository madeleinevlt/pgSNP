###
# Créer un fichier d'alignement pour chaque contig
###

#!/usr/bin/python3

import argparse
import glob
from Bio import AlignIO
import subprocess
import os

def get_argv() :
    parser = argparse.ArgumentParser(description='Take a multifasta from vcfs2superalignement to make the alignment of all genome ')
    parser.add_argument("file_in", metavar='input', type=str, help='multifasta file')
    parser.add_argument("-name_outfile", default="result_contig.fasta", type=str, help='name of output')
    return parser.parse_args()

def read_genome(genome_file):
## arbo : nom_contig{nom_genome : ['ATATATAGA']}
    dico_contigs={}
    dico_genome={}
    old_name_contig=""
    flag=0
    with open(genome_file) as fillin:
        for line in fillin :
            if line.startswith(">"):
                name_contig=line.strip(">").strip("\n").split(" ")[-1]
                name_genome=line.strip(">").strip("\n").split(" ")[0]
                if old_name_contig=="" :
                    old_name_contig=name_contig
                if name_contig!=old_name_contig :
                    if old_name_contig in dico_contigs :
                        dico_contigs[old_name_contig].append(dico_genome)
                    else :
                        dico_contigs[old_name_contig]=[dico_genome]
                    dico_genome={}
                    old_name_contig=name_contig

            else :
                dico_genome[name_genome]=line.strip("\n")
        dico_contigs[name_contig].append(dico_genome)
    #print(dico_contigs["17Q002737NODE_19_length_89388_cov_41.48267375667"])
    return(dico_contigs)
def write_genome(dico, folder):

    for contig in dico :
        #print(contig["17Q002737NODE_19_length_89388_cov_41.48267375667"])
        #print(contig)
         #minimum pour que iqtree fasse un arbre avec boostrap
        with open(folder + contig +"_aligned.fasta","w") as fillout :
            for genome in dico[contig] :
                #print(contig)
                #print("GENOME")
                #print(genome)
                for elem in genome :
                    #print("ELEMENT")
                    #print(elem)
                    #print(genome[elem])
                    #break
                    #break
                    #seq=dico[contig][genome]
                    #list_key=sorted(list(dico[list(dico.keys())[0]].keys()))
                    if len(genome[elem].replace("-","")) !=0:
                        fillout.write(">" + elem + "\n")
                        fillout.write(genome[elem] +"\n")

                    #if len(genome[elem].replace("-","")) > (len(genome[elem]) /100*10):
                    #    fillout.write(">" + elem + "\n")
                    #    fillout.write(genome[elem] +"\n")

def veref_aligment (files) :
    for file in files :
        if os.stat(file).st_size != 0 :
            alignment = AlignIO.read(file, "fasta")
            if len(alignment) < 4 : #taile mini pour iqtree
                subprocess.call("rm " + file, shell=True)
            #else :
            #    subprocess.call("python3 /global/scratch/m.desousaviolante/mbandaka/pgSNP_mbandaka/prune_aln_cols.py --all-gap -i" + file + "> " + file + "_gap_free", shell=True)
        else :
            subprocess.call("rm " + file, shell=True)


if __name__ == "__main__":
    args = get_argv()
    dico=read_genome(args.file_in)

    write_genome(dico,args.name_outfile)
    veref_aligment(glob.glob("*aligned.fasta"))
