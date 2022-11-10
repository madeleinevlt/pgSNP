#!/usr/bin/python3

import argparse
import glob
import subprocess


def get_argv() :
    parser = argparse.ArgumentParser(description='Change iqtree mldist file to input for Erable. Mldist and alignement file MUST BE in the same order (ex : -i treeA.mldist treeB.mldist -a treeA.alignement treeB.alignement')
    parser.add_argument("-i", metavar='input', nargs="+", type=str, help='iqtree mldist files')
    parser.add_argument("-a", metavar='alignment', nargs="+", type=str, help='alignement file')
    parser.add_argument("-o", metavar='output', type=str, help='output file')
    return parser.parse_args()

def get_all_matrix_with_alignement_length(list_file, lengths, output_file) :
    print("get_all_matrix")
    with open(output_file,"w") as fillout : 
        fillout.write("\n" + str(len(list_file)) + "\n\n")
        i=0
        for file in list_file :
            with open(file, "r") as fillin :
                nb_ech=fillin.readline().strip()
                #fillout.write(nb_ech + " " + "1" + "\n")
                fillout.write(nb_ech + " " + str(lengths[i]) + "\n")
                for line in fillin :
                    fillout.write(line)
            
            fillout.write("\n")
            i=i+1    
            
    
def verif_matrix(list_file, lengths) :
    res=[]
    res_l=[]
    j=0
    for file in list_file :
        pass_file = True
        #print(file)
        with open(file, "r") as fillin :
            fillin.readline()
            for line in fillin :
                l=line.strip().split()
                #print(l)
                l.pop(0)
                for i in l :
                    if float(i) > 1.5 :
                        pass_file=False
                        print(file, " FAILED")
        if pass_file :
            res.append(file)
            res_l.append(lengths[j])
        j=j+1    
    return res, res_l
                        
  
def get_all_length(list_alignement) :
    lengths = []
    for file in list_alignement :
        with open(file, "r") as fillin :
            fillin.readline()
            lengths.append(len(fillin.readline().strip()))
    return lengths        
    
    
if __name__ == "__main__": 
    args = get_argv()
    l = get_all_length(args.a)
    print(l)
    files, lengths = verif_matrix(args.i, l)
    print(len(files))
    print(len(lengths))
    get_all_matrix_with_alignement_length(files, lengths, args.o)
    
