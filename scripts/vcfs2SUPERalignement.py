#!/usr/bin/python3

## AU NIVEAU DE L'ARBORECENCE DE FICHIER :
#  chaque VCF est dans un dossier de résultat SNIPPY
#  tous les dossiers de résultats sont dans un même dossier, spécifié en input
## 

import argparse
import glob
from Bio import Align
import subprocess


def get_argv() :
    parser = argparse.ArgumentParser(description='Merge VCF to create a fasta of all_contigs')
    parser.add_argument("file_in", metavar='input', type=str, help='folder which contains ALL SNIPPY folders')
    parser.add_argument("-o", default="result.fasta", type=str, help='name of output')
    parser.add_argument("-nom_a_suppr", default="_snip", type=str, help='si vos noms de dossiers ont des choses ajoutés')
    parser.add_argument("--only_snp", action='store_true',help='Si oui, alors ne gère que les SNP (et non pas les insertions/deletions)')
    
    return parser.parse_args()

    

def read_VCF(VCF_files) :
## arbo : nom_vcf{nom_contig{position : ['listes des informations']}
## contient QUE les insertions
    dictionnaire_contigs={} ## stock les noms de contigs
    
    with open(VCF_files) as fillin:
        for line in fillin :
            if not line.startswith("#"):
                line_splitted = line.strip("\n").split("\t")
                if len(line_splitted[3])<len(line_splitted[4]):
                    if len(line_splitted[3])> 1 : ##si c'est un complexe en faite les SNP sont déjà changé mais il manque l'insertion du dernier
                        line_splitted[1]=int(line_splitted[1])+len(line_splitted[3])-1
                        line_splitted[3]=line_splitted[3][-1]
                        line_splitted[4]=line_splitted[4][-2:]
                        
                
                    if not line_splitted[0] in dictionnaire_contigs:
                        dictionnaire_contigs[line_splitted[0]]={line_splitted[1]:line_splitted[3:5]}
                    else:
                        dictionnaire_contigs[line_splitted[0]][line_splitted[1]]=line_splitted[3:5]
    #print(dictionnaire_contigs)                
    return(dictionnaire_contigs)

def read_genome(genome_file):
## arbo : nom_genome{nom_contig : ['ATATATAGA']}

    dico_contigs={}
    flag=0
    with open(genome_file) as fillin:
        for line in fillin :
            if line.startswith(">"):
                if flag==1:
                    dico_contigs[name_contig]=sequence
                name_contig=line.strip(">").strip("\n")
                sequence=""
                flag=1
            else :
                sequence=sequence+line.strip("\n")
        dico_contigs[name_contig]=sequence
    return(dico_contigs)
  

                
def constru_seq_with_INDEL(all_VCF,genome_files, list_en_DOUBLE, output_name):
##arbo : nom_genome{contig : ['ATATA']}
    with open("taille_contigs","w") as fillout_len_seq :
        with open(output_name,"w") as fillout : 
            taille_seq=[]
            for genome in genome_files :  # Pour chaque genome dans les génomes
                print("GENOME " +genome)
                le_premier=""
                for contig in genome_files[genome]: # Pour chaque contig dans le génome
                    print("Contig " +contig)                   
                    next=""
                    new_position=0
                    seq=genome_files[genome][contig].replace("N","-")  

                    print(all_VCF[contig])
                    if contig in all_VCF:  # Si un contig a des variants 
                        flag=""
                        i=0
                        list_double_done=[]
                        ## format de ins : position, ref, alt, genome qui a cette variation, ref aligné, alt aligné
                        for ins in all_VCF[contig] : # Pour tous les variants dans le contig
                            print("ins :")
                            print(ins)
                            ins[0]=int(ins[0])
                            index=int(ins[0])-1 + new_position 
                            ref=ins[1]
                            alt=ins[2]
                            
                            if contig in list_en_DOUBLE : # si le contig est dans la liste en double 
                                if ins[0] in list_en_DOUBLE[contig] : 
                                    if i==0 or ins[0] != all_VCF[contig][i-1][0]: 
                                        flag="premier"
                                    if flag=="premier":
                                        j=i
                                        plus_grand=len(ins[5])
                                        while int(ins[0]) == int(all_VCF[contig][j+1][0]) :
                                            if len(all_VCF[contig][j+1][5]) > plus_grand :
                                                plus_grand=len(all_VCF[contig][j+1][5])
                                            flag="yaunplusgrand"    
                                            j=j+1
                                            if len(all_VCF[contig])==j+1 :
                                                break
                                        if flag=="yaunplusgrand" :
                                            j=i
                                            ins[4]=ins[4] + (plus_grand - len(ins[4])) *"-"
                                            ins[5]=ins[5] + (plus_grand - len(ins[5])) *"-"
                                            while int(ins[0]) == int(all_VCF[contig][j+1][0]): #plutot un alignement
                                                all_VCF[contig][j+1][4] = all_VCF[contig][j+1][4] + (plus_grand - len(all_VCF[contig][j+1][4])) * "-"
                                                all_VCF[contig][j+1][5] = all_VCF[contig][j+1][5] + (plus_grand - len(all_VCF[contig][j+1][5])) * "-"
                                                j=j+1
                                               
                                                if len(all_VCF[contig])==j+1 :
                                                    break
                                        flag=""        
                                    
                                    if genome in ins[3]:
                                        list_double_done.append(int(ins[0]))
                                    
                                    else :
                                        if int(ins[0]) in list_double_done :
                                            flag="DOUBLEPLUSTARD"
                                        elif len(all_VCF[contig])!=i+1 :
                                            if int(ins[0]) == int(all_VCF[contig][i+1][0]) :
                                                flag="DOUBLEPLUSTARD"

                                                                           
                                
                            if flag!="DOUBLEPLUSTARD" : # si on peut traiter le variant mtn

                                if genome in ins[3]:  # si le genome a l'insertion, ou si c'est juste une VICTIME d'une insertion dans un autre 
                                    a_remplacer=ins[5]
                                else :
                                    a_remplacer=ins[4]
                                
                                if len(ref) > 1 : 
                                    to_search=seq[index:index+len(ref)]
                                else :
                                    to_search=seq[index]
                                
                                if to_search == ref : #si j'ai bien trouvé la position
                                    new_seq=seq[:int(index)]+a_remplacer+seq[int(index)+1:]
                                    seq=new_seq
                                    next=next+a_remplacer
                 
                                else :
                                    if to_search == alt[0] :
                                        a_remplacer = ins[5]
                                        new_seq=seq[:int(index)]+a_remplacer+seq[int(index)+1:]
                                        seq=new_seq
                                        next=next+a_remplacer

                                    elif to_search == "-" or to_search== "n" :
                                        new_seq=seq[:int(index)]+len(a_remplacer)*to_search+seq[int(index)+1:]
                                        seq=new_seq
                                        next=next+len(a_remplacer)*"X"
                                             
                                    else :
                                        
                                        aligned_ref, pattern, aligned_alt=make_alignment(to_search, a_remplacer)

                                        if a_remplacer==ins[4] : #si normalement c'est la ref, faut changer la ref

                                            new_seq=seq[:int(index)]+aligned_ref+seq[int(index)+1:]
                                            a_remplacer=aligned_ref
                                            
                                        else :
                                            new_seq=seq[:int(index)]+aligned_alt+seq[int(index)+1:]
                                            a_remplacer=aligned_alt
                                        seq=new_seq
                                        next=next+a_remplacer

                                new_position=new_position+len(a_remplacer)-1

                                taille_seq.append(len(seq))    

                            i=i+1
                            flag=""

                        fillout_len_seq.write(genome + " " + contig + " " + str(len(seq)) + "\n") 


                    fillout.write(">" + genome + " " + contig +"\n")
                    seq=seq.replace("X","-").replace("N","-").replace(".","-")
                    fillout.write(seq.replace("n","N") + "\n")

                   
def get_seq_aligned(alignment):
##pour récupérer le target et la query avec les gaps
##renvoit seq1, le pattern et seq2 alignés de façon à pouvoir récupérer les elem 1 par 1
    query = alignment.query 
    target = alignment.target 
     
    seq1 = str(target) 
    seq2 = str(query) 
    n1 = len(seq1) 
    n2 = len(seq2) 
    aligned_seq1 = "" 
    aligned_seq2 = "" 
    pattern = "" 
    path = alignment.path 
    end1, end2 = path[0] 
    if end1 > 0 or end2 > 0: 
        end = max(end1, end2) 
        aligned_seq1 += "." * (end - end1) + seq1[:end1] 
        aligned_seq2 += "." * (end - end2) + seq2[:end2] 
        pattern += '.' * end 
    start1 = end1 
    start2 = end2 
    for end1, end2 in path[1:]: 
        gap = 0 
        if end1 == start1: 
            gap = end2 - start2 
            aligned_seq1 += '-' * gap 
            aligned_seq2 += seq2[start2:end2] 
            pattern += '-' * gap 
        elif end2 == start2: 
            gap = end1 - start1 
            aligned_seq1 += seq1[start1:end1] 
            aligned_seq2 += '-' * gap 
            pattern += '-' * gap 
        else: 
            s1 = seq1[start1:end1] 
            s2 = seq2[start2:end2] 
            aligned_seq1 += s1 
            aligned_seq2 += s2 
            for c1, c2 in zip(s1, s2): 
                if c1 == c2: 
                    pattern += '|' 
                else: 
                    pattern += 'X' 
        start1 = end1 
        start2 = end2 
    n1 -= end1 
    n2 -= end2 
    n = max(n1, n2) 
    aligned_seq1 += seq1[end1:] + '.' * (n - n1) 
    aligned_seq2 += seq2[end2:] + '.' * (n - n2) 
    pattern += '.' * n 
    return(aligned_seq1, pattern, aligned_seq2)

def transition(nucleotide):
    if nucleotide == "A":
        return "G"
    elif nucleotide == "G":
        return "A"
    elif nucleotide == "C":
        return "T"
    else :
        return "C"
        
def make_alignment(seq1,seq2):
##prend en entrée 2 séquence et renvoit son alignment
## remarque : je ne traite QUE des insertions. Donc j'ai forcement un nucléotide en commun dans REF et dans ALT, sauf si j'ai un SNP a cette position dans la séquence
## donc je ne dois pas oublier de le vérifier quand je refais l'alignement
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -10
    
    alignments = aligner.align(seq1,seq2)
    
    if alignments : ##ya bien un alignement
        alignment = sorted(alignments)[0]
        return get_seq_aligned(alignment)
    
    else : ##yen a pas donc je dois le faire à la main
        transi = transition(seq1)
        if transi in seq2:
            index=seq2.find(transi)
            aligned_seq1="-" * index + seq1 + "-" *  (len(seq2)-index-1)
            return(aligned_seq1, "", seq2)
            
        else : #ya meme pas de transition, donc c'est un alignement qui va être naze donc je met en 1ere position et basta
            aligned_seq1=seq1 + "-" * (len(seq2)-1)
            return(aligned_seq1, "", seq2)
            
            
def treat_INDEL(all_VCF_files):
#permet de faire une belle liste avec tous les INDELS pour chaque contigs, avec la liste des génomes concernés
    list_INDEL={}
    list_endouble={}
    for genome in all_VCF_files:
        for contig in all_VCF_files[genome]:
            #print(contig)
            for pos in all_VCF_files[genome][contig]:
                indels  = all_VCF_files[genome][contig][pos]
                if not list_INDEL: ## si la liste est vide
                    #make alignment
                    aligned_seq1,pattern,aligned_seq2=make_alignment(indels[0],indels[1])
                    list_INDEL[contig]=[[pos, indels[0], indels[1], [genome], aligned_seq1, aligned_seq2]]
                
                elif contig not in list_INDEL: ##si le contig n'a pas encore été ajouté
                    #make alignment
                    aligned_seq1,pattern,aligned_seq2=make_alignment(indels[0],indels[1])
                    list_INDEL[contig]=[[pos, indels[0], indels[1], [genome], aligned_seq1, aligned_seq2]]
               
                else :
                    exist=False
                    for elem in list_INDEL[contig]:
                        if elem[0]==pos: 
                            if elem[2]==indels[1] :
                                elem[3].append(genome)
                                exist=True
                            else : #position en double
                                if not list_endouble :
                                    list_endouble[contig]=[pos]
                                elif contig not in list_endouble :
                                    list_endouble[contig]=[pos]
                                else :
                                    list_endouble[contig].append(pos)
                    if exist==False :
                        #make alignment
                        aligned_seq1,pattern,aligned_seq2=make_alignment(indels[0],indels[1])
                        list_INDEL[contig].append([pos, indels[0], indels[1], [genome], aligned_seq1, aligned_seq2])
    return(list_INDEL, list_endouble)
        
def sort_list_INDEL(list_INDEL):
    new_list={}
    for contig in list_INDEL :
        list_position=[]
        for ins in list_INDEL[contig]:
            list_position.append(int(ins[0]))
        position_sorted=sorted(list_position)
        
        for position in position_sorted:
            for ins in list_INDEL[contig] :
                if int(ins[0])==position :
                    if not new_list : 
                        new_list[contig]=[ins]
                    elif contig not in new_list:
                        new_list[contig]=[ins]
                    else :
                        new_list[contig].append(ins)
                    del list_INDEL[contig][list_INDEL[contig].index(ins)]
                    break
    return(new_list)
                
def constru_seq_only_SNP(genome_files, output_name):
    with open(output_name,"w") as fillout :
        for genome in genome_files :  # Pour chaque genome dans les génomes
            for contig in genome_files[genome]:
            
                seq=genome_files[genome][contig].replace("N","-")
                fillout.write(">" + genome + " " + contig +"\n")
                seq=seq.replace("X","-").replace("N","-").replace(".","-")
                fillout.write(seq.replace("n","N") + "\n")
        
    
if __name__ == "__main__": 
    args = get_argv()
    
    all_VCF_files={}
    all_genome={}
    
    print(len(glob.glob(args.file_in + "*/snps.vcf")))
    for VCF_files in glob.glob(args.file_in + "*/snps.vcf"):
        all_VCF_files[VCF_files.split("/")[-2].split(args.nom_a_suppr)[0]]=read_VCF(VCF_files)

    #print(all_VCF_files)
    print(len(glob.glob(args.file_in + "*/snps.aligned.fa")))
    for genome_file in glob.glob(args.file_in + "*/snps.aligned.fa"):
        #subprocess.call("cp " + genome_file + " " + genome_file + "_withINDEL",shell=True)
        all_genome[genome_file.split("/")[-2].split(args.nom_a_suppr)[0]]=read_genome(genome_file)
    
    print(args.only_snp)
    ##INDELs gestion -> DO NOT USE
    if not args.only_snp : 
        list_INDEL, list_en_DOUBLE=treat_INDEL(all_VCF_files)
        list_INDEL_sorted=sort_list_INDEL(list_INDEL)
        constru_seq_with_INDEL(list_INDEL_sorted, all_genome, list_en_DOUBLE, args.o)
    else :
        ##SNP gestion
        constru_seq_only_SNP(all_genome, args.o)

    
    
    
