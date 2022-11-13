#!/usr/bin/python3

## 1 ETAPE : Prendre une séquence et la considérer comme notre référence
## 2 ETAPE : Pour chaque fichier :
##              -lancer un BLAST
##              - récupérer les fragments de contig > à une taille qui appartient pas à un hsp
##              - récupérer les contigs qui n'ont pas de hit et les fragments de séquence > à une taille
##              - les ajouter à la référence
#############################

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO
import argparse
import subprocess
import glob
import os
from Bio import Seq
from Bio import SeqRecord

def get_argv() :
    parser = argparse.ArgumentParser(description='Build a reference from contigs files')
    parser.add_argument("-i", metavar='input', type=str, help='folder which contains your contigs or file which contains all path for contig location')
    parser.add_argument("-nb_contigs", default=300, type=int, help='minimum of contigs in files, default=300')
    parser.add_argument("-hit_min", default=500, type=int, help='minimum sequence length for a non hit contig, default=500')
    parser.add_argument("-identity", default=95, type=float, help='identity, default=95')
    parser.add_argument("-o", metavar='output', type=str, help='folder which contains all outputs. Must be created')
    parser.add_argument("-name_ref", type=str, default="ref.fasta", help='name of the final reference')
    parser.add_argument("-name_log", type=str, default="result", help="name of log output")
    parser.add_argument("-output_tmp", type=str, default="output_contigs_good_order", help="name of a tmp output")
    parser.add_argument("-default_path", type=str, default="/global/bio/data/GAMeR_DB/SALMONELLA/", help="default path of folder if not specifed in file_in")
    return parser.parse_args()


def write_output(output_folder, name_log, contribution_contigs, contribution_nb_bp, keys):
    with open(name_log + ".log", "w") as fillresult :
        for i in range(0,len(keys)) :
            fillresult.write(keys[i] +"," + str(contribution_contigs[i]) +"," + str(contribution_nb_bp[i]) + "\n")


def read_contigs(all_files,nb_contigs, hit_mini) :
    dict_name_seq = {}
    list_name = []
    list_length = []
    num=0
    for file in all_files :
        if len(list(SeqIO.parse(file, "fasta"))) <= nb_contigs :
            for record in SeqIO.parse(file, "fasta"):
                name=file.split("/")[-1].split("_contigs")[0] + str(record.id) + str(num)
                num=num+1
                dict_name_seq[name]=str(record.seq)
                list_name.append(name)
                list_length.append(len(record.seq))

    return dict_name_seq, list_name, list_length


def make_contig_file(dict_name_seq , list_name , list_length , len_hit_mini, output_tmp) :
    keys=[]
    contribution_contigs=[] #contribution pour chaque fichier en nombre de contigs dans la référence
    contribution_nb_bp=[] # contribution pour chaque fichier en nombre de bp dans la référence


    tri_sorted=sorted(list_length)
    with open(output_tmp,"w") as fillout :
        while tri_sorted :
            number=tri_sorted.pop()
            recup_id=list_length.index(number)
            recup_name=list_name[recup_id]
            list_length.pop(recup_id)
            list_name.pop(recup_id)
            if len(dict_name_seq[recup_name]) >= len_hit_mini :
                fillout.write(">" + recup_name + "\n")
                fillout.write(dict_name_seq[recup_name] + "\n")


def gestion_fragment(ref_in,potentiel_inser) :
    print("début de la gestion des fragments : \n\n")
    p=[]
    print(potentiel_inser.keys())
    for seq_record in SeqIO.parse(ref_in, "fasta"):
        print(seq_record.id)
        if seq_record.id in potentiel_inser.keys() :
            print("id trouve !!!!")
            list_pos=list(potentiel_inser[seq_record.id].keys())
            list_pos.sort(reverse=True)
            seq=str(seq_record.seq)
            print("before :" + str(len(seq)))
            for position in list_pos :
                seq=seq[:position]+potentiel_inser[seq_record.id][position]+seq[position:]
            seq_record.seq=Seq.Seq(seq)
            print("after :" + str(len(seq)))
            print(potentiel_inser)
        p.append(seq_record)
    if p :
        SeqIO.write(p,ref_in,"fasta-2line")


def make_reference(output_tmp, nb_contigs, len_hit_mini, IDENTITY, output_folder, ref_name, name_log):

    keys=[]
    contribution_contigs=[] #contribution pour chaque fichier en nombre de contigs dans la référence
    contribution_nb_bp=[] # contribution pour chaque fichier en nombre de bp dans la référence
    nb_contig_all=0
    ref_in = output_folder + "/" + ref_name


    #creation reference
    with open(output_tmp,"r") as fillin :
        with open(ref_in,"w") as fillout :
            fillout.write(fillin.readline())
            seq=""
            for line in fillin :
                if line[0]!=">" :
                    seq=seq+line
                else :
                    break
            fillout.write(seq)

    tmp_fasta=output_folder + "tmp.fasta"
    tmp_xml=output_folder + "test_sortie.xml"
    for record in SeqIO.parse(output_tmp, "fasta"):


        SeqIO.write(record,tmp_fasta,"fasta-2line")

        blastn_cline = NcbiblastnCommandline(query=tmp_fasta, subject=ref_in, out=tmp_xml, perc_identity=IDENTITY, outfmt=5)

        blastn_cline() #lancer un blast
        print("record numéro :" + record.id)
        #print(timetime)
        result_handle = open(tmp_xml)
        blast_record = NCBIXML.read(result_handle)
        q_dict = SeqIO.index(tmp_fasta, "fasta") # contient le contig entier
        hits = [] # liste qui va récupérer les nom de contigs qui n'ont pas de hits

        nb_potentiel=0 # afin d'avoir le nb de bout de sequences
        nb_pb_potentiel=0 # taille total de tous les nb de bouts de sequences que j'ajoute

        dico_bornes={} #dictionnaire qui contient tous les débuts et fin de chaque alignement -> permet de repérer les morceaux qui ne sont pas hsp
        dico_bornes_hit={} #dictionnaire là où c'est la fin du hit sur la référence + nom du contig
        if blast_record.alignments :
            hits.append(blast_record.query.split()[0])
            for alignment in blast_record.alignments:   #pour chaque alignement d'un contig
                for hsp in alignment.hsps:      # pour chaque hsp
                    print("identite :" + str(hsp.identities))
                    if hsp.query_start in dico_bornes :    # on ajoute le début et la fin de l'alignement dans un dictionnaire
                        max_end_query = max(dico_bornes[hsp.query_start][0], hsp.query_end)
                        if max_end_query == hsp.query_end :
                            dico_bornes[hsp.query_start]=[max_end_query, hsp.sbjct_end, alignment.hit_id]
                    else :
                        dico_bornes[hsp.query_start]= [hsp.query_end, hsp.sbjct_end, alignment.hit_id]

            list_debut=list(dico_bornes.keys())
            print(dico_bornes)
            list_debut.sort() # le début de chaque alignement est stocké dans une liste pour pouvoir la sort (car pas possible dans dictionnaire).
            print(list_debut)

            if list_debut :
                potentiel={}    # liste qui contient les contigs qui peuvent potentiellement etre unique et donc à ajouter à la reference
                potentiel_inser={} # liste qui contient les insertions qui peuvent potentiellement etre unique et donc à ajouter à la reference
                if list_debut[0] > 1 : #regarder au début si on peut prendre des fragments entre [:premierdebut]
                    potentiel[0] = [list_debut[0],0, dico_bornes[list_debut[0]][2]]

                while len(list_debut)!=0 :
                    debut1=list_debut[0]
                    if len(list_debut)>1 : #si ya + ou = que 2 elem
                        debut2=list_debut[1]
                        if dico_bornes[debut1][0] > debut2 : # si fin seq 1 > debut seq 2
                            if dico_bornes[debut1][0] < dico_bornes[debut2][0] : # si la fin de la seq 1 est plus petite que la fin de la seq 2
                                dico_bornes[debut1][0] = dico_bornes[debut2][0]
                            del list_debut[1]   # on delete seq 2 (car on a toutes les infos necessaires dans seq 1)
                        else : #sinon : (donc si on a un espace entre seq 1 et seq 2)
                            potentiel[dico_bornes[debut1][0]]= [debut2,dico_bornes[debut1][1],dico_bornes[debut1][2]]
                            del list_debut[0]
                    else : # s'il nous reste 1 dernier element faut regarder si ya un espace entre le dernier et la fin de la sequence totale
                        potentiel[dico_bornes[debut1][0]]= [len(q_dict[blast_record.query.split()[0]]),dico_bornes[debut1][1],dico_bornes[debut1][2]]
                        del list_debut[0]


                print(potentiel)
                if potentiel : #si on a trouve des fragments qui ont une taille suffisante pour être ajouté à la ref
                    file_exist=False
                    for key in potentiel :
                        to_add=str(q_dict[blast_record.query.split()[0]].seq)[key:potentiel[key][0]]
                        if len(to_add) < len_hit_mini :
                            print("taille du contig à ajouter :")
                            print(len(to_add))
                            who_add=potentiel[key][2] # sequence sur laquelle on va ajouter l'insertion
                            where_add=potentiel[key][1] # a quel position il faut ajouter l'insertion
                            if who_add not in potentiel_inser :
                                potentiel_inser[who_add] = {where_add : to_add}
                            else :
                                potentiel_inser[who_add][where_add] = to_add

                        else :
                            print("contig trop grand pour etre une insertion ")
                            file_exist=True
                            with open(output_folder + "potentiel_contigs_ref","a") as fillout :
                                fillout.write(">" + blast_record.query.split()[0] + str(key))
                                fillout.write("\n")
                                fillout.write(str(q_dict[blast_record.query.split()[0]].seq)[key:potentiel[key][0]])
                                fillout.write("\n")
                                nb_contig_all=nb_contig_all+1
                        nb_potentiel=nb_potentiel+1
                        nb_pb_potentiel=nb_pb_potentiel+len(str(q_dict[blast_record.query.split()[0]].seq)[key:potentiel[key][0]])

                    if file_exist :
                        subprocess.call("cat " + output_folder + "/potentiel_contigs_ref>> " + ref_in, shell=True)
                        subprocess.call("rm " + output_folder + "/potentiel_contigs_ref", shell=True)
                #pour les contigs qui n'ont pas de hit sur la sequence :
            misses = set(q_dict.keys()) - set(hits)
            print(misses)
            orphan_records = [q_dict[name] for name in misses]

            orphan_length_ok = []
            nb_bp=0
            for seq in orphan_records :
                nb_contig_all=nb_contig_all+1
                seq.id = record.id + seq.id
                if len(seq) > len_hit_mini : # vérification de la taille des contigs qui n'ont pas hit
                    orphan_length_ok.append(seq)
                    nb_bp=nb_bp+len(seq)
            tmp_record=output_folder + "record_mauvais"
            SeqIO.write(orphan_length_ok, tmp_record,"fasta")
            subprocess.call("cat " + tmp_record + " >> " + ref_in, shell=True)


            if potentiel_inser : #s'il y a potentiellement des insertions
                gestion_fragment(ref_in,potentiel_inser)
            keys.append(record.id)
            contribution_contigs.append(len(orphan_length_ok) + nb_potentiel)
            contribution_nb_bp.append(nb_bp + nb_pb_potentiel)
        #pas de hit, donc ajout du contig entier
        else :
            subprocess.call("cat " + tmp_fasta + " >> " + ref_in, shell=True)
            print("adding all the contig")

        #break
    p = subprocess.Popen('grep -c ">" '+ ref_in ,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #nombre de contigs dans la référence au départ
    out, err = p.communicate()
    print("number of contigs in reference :" + str(int(out)))
    print("numero contig cense avoir :" + str(nb_contig_all))
    write_output(output_folder, name_log, contribution_contigs, contribution_nb_bp, keys)


def read_snippy(path_reads) :
    format_r=""
    ids=[]
    for r in path_reads :
        ids.append(os.path.basename(os.path.normpath(r)))
        if r[-2:]=="gz" :
            if format_r :
                if format_r=="fastq" :
                    print("Format problem, please choose between fastq and fastq.gz format for your reads in the folder")
                    exit()
            format_r="gz"
        else :
            format_r="fastq"
            if format_r :
                if format_r=="gz" :
                    print("Format problem, please choose between fastq and fastq.gz format for your reads in the folder")
                    exit()
    return(ids,format)


def read_path_file(file, default_path) :
    with open(file, "r") as fillin :
        l=[]
        for line in fillin :
            fold = os.path.basename(os.path.normpath(line.strip()))
            if not os.path.exists(line.strip()) :
                l.append(default_path + fold + "/" + fold + "_contigs.fasta")
            else :
                if "_contigs.fasta" not in line :
                    l.append(line.strip() + "/" + fold + "_contigs.fasta")
                else :
                    l.append(line.strip())
    return l


if __name__ == "__main__":
    args = get_argv()

    if os.path.isfile(args.i) == True :
        paths=read_path_file(args.i, args.default_path)
        dict_name_seq, list_name, list_length = read_contigs(paths, args.nb_contigs, args.hit_min)

    elif os.path.isdir(args.i) :
        dict_name_seq, list_name, list_length = read_contigs(glob.glob(args.i +"*_contigs.fasta"),args.nb_contigs, args.hit_min)
    else :
        print("Format incorrect. Please give a folder with all of your contig file or a file which contains all paths to your contigs files")
        exit()

    make_contig_file(dict_name_seq , list_name , list_length , args.hit_min, args.output_tmp)
    make_reference(args.output_tmp , args.nb_contigs, args.hit_min,args.identity, args.o, args.name_ref, args.name_log)
