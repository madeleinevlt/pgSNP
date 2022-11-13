#!/usr/bin/python3


#pgSNP pipeline
import argparse
import subprocess
import glob
import os
import scripts.make_pangenome as mp
import scripts.vcfs2alignement as va
import scripts.merge_contig as mc
import multiprocessing
import scripts.matrix_mldist_to_erable_input as me

def get_argv() :
    parser = argparse.ArgumentParser(description='Build a reference from contigs files')
    parser.add_argument("-contig", type=str, help='folder which contains your contigs or file which contains all path for contig location',required=True)
    parser.add_argument("-reads", type=str, help='folder which contains your reads',required=True)
    parser.add_argument("-o", metavar='output', type=str, help="folder which contains all outputs. Must be created",required=True)
    parser.add_argument("-name_res", type=str, default="pangenome", help='name of the finals outputs')


    parser.add_argument("-tmp", default=False , type=bool, help="Keep all tempory files, default=false")
    parser.add_argument("-cpus", default=12 , type=int, help="Cpus available")


##Pangenome arguments
    parser.add_argument("-nb_contigs", default=300, type=int, help='minimum of contigs in files, default=300') ##this number can be changed if the contig file contains a higher number of contigs
    parser.add_argument("-hit_min", default=500, type=int, help='minimum sequence length for a non hit contig, default=500')
    parser.add_argument("-identity", default=95, type=float, help='identity, default=95')
    parser.add_argument("-name_log", type=str, default="result", help="name of log output")
    parser.add_argument("-output_tmp", type=str, default="output_contigs_good_order", help="name of a tmp output")
##snippy arguments
    parser.add_argument("-name_snip_folder", type=str, default="_snip", help="output_snippy_name")

##alignment arguments
    parser.add_argument("-output_alignment_name", default="result.fasta", type=str, help='name of output')
    parser.add_argument("-tmp_outfile", default="result_contig.fasta", type=str, help='name of output')
##tree arguments
    parser.add_argument("-model", default="TEST", type=str, help='tree model of IQTREE. Default = TEST, tested also on GTR+I+R+F on Salmonella')

    return parser.parse_args()


if __name__ == "__main__":
    args = get_argv()
    if not os.path.isdir(args.reads) :
        print("Reads forlder incorrect. Please give a folder with all of your reads file")
        exit()

    if args.cpus > multiprocessing.cpu_count() :
        print("ERROR: You have specified more threads than CPU cores available, it will fail at IQ-TREE step. CPUs available : " + str(multiprocessing.cpu_count()))
        exit()
##Pangenome pipeline
    if os.path.isfile(args.contig) == True :
        paths=read_path_file(args.contig, args.default_path)
        dict_name_seq, list_name, list_length = mp.read_contigs(paths, args.nb_contigs, args.hit_min)

    elif os.path.isdir(args.contig) :
        dict_name_seq, list_name, list_length = mp.read_contigs(glob.glob(args.contig +"*_contigs.fasta"),args.nb_contigs, args.hit_min)
    else :
        print("Incorrect format. Please give a folder with all of your contig file or a file which contains all paths to your contigs files")
        exit()

    args.o=args.o+"/"

    #mp.make_contig_file(dict_name_seq , list_name , list_length , args.hit_min, args.o + args.output_tmp)
    #mp.make_reference(args.o + args.output_tmp , args.nb_contigs, args.hit_min,args.identity, args.o, args.name_ref + ".fasta", args.o + args.name_log)


    tmp_list=args.o + "tmp.fasta" + " " + args.o + "record_mauvais" + " " + args.o + "test_sortie.xml" + " " + args.o + "log_reference_all_data.log.log" + " " + args.o + "result.log" + args.o + args.output_tmp
    if not args.tmp :
        subprocess.call("rm " + tmp_list, shell="True")

##Snippy pipeline
    if os.path.isdir(args.reads) :
        p_reads=glob.glob(args.reads + "*.fastq.gz")
        if not p_reads :
            p_reads=glob.glob(args.reads + "*.fastq")
            if p_reads :
                print("Read format incorrect. Reads in fastq format or fastq.gz format. Read folder : " + args.reads)
                exit()

        ids,format=mp.read_snippy(p_reads)
        ids_sorted=sorted(ids)
        if len(ids_sorted) % 2 ==1 :
            print("Incorrect reads pairs. Please check which read doesn't have its pairs (eg. R1 or R2 lacking)")
            exit()
        while len(ids_sorted) > 1 :
            R1=ids_sorted[0]
            R2=ids_sorted[1]
            name_folder=R1.split("_R1")[0]+args.name_snip_folder
            #subprocess.call("bin/snippy/bin/snippy --cpus "+ str(args.cpus) + " --outdir "+ args.o + "/" + name_folder + " --ref " + args.name_ref + ".fasta --R1 " + args.reads + R1 + " --R2 " +  args.reads + R2, shell=True)
            del ids_sorted[0]
            del ids_sorted[0]
##Alignment marker
    all_VCF_files={}
    all_genome={}
    for VCF_files in glob.glob(args.o + "*/snps.vcf"):
        all_VCF_files[VCF_files.split("/")[-2].split(args.name_snip_folder)[0]]=va.read_VCF(VCF_files)
    for genome_file in glob.glob(args.o + "*/snps.aligned.fa"):
        all_genome[genome_file.split("/")[-2].split(args.name_snip_folder)[0]]=va.read_genome(genome_file)
    va.constru_seq_only_SNP(all_genome, args.o + args.output_alignment_name)
    #dico=mc.read_genome(args.o + args.output_alignment_name)

    #mc.write_genome(dico,args.o)
    #mc.veref_aligment(glob.glob(args.o+"*_aligned.fasta"))
    if not args.tmp :
        subprocess.call("rm -rf "+ args.o + "*" + args.name_snip_folder, shell="True")

## trees

    #for file in glob.glob(args.o+"*_aligned.fasta") :
        #subprocess.call("iqtree -s " + file +" -m " + args.model + " -bb 1000 -nt " + str(args.cpus),shell=True)
    subprocess.call("cat " + args.o + "*.treefile > " + args.o + "all_iqtree_files && bin/FastRFS-linux/FastRFS -i " + args.o + "all_iqtree_files -o " + args.o + "file_fastRFS",shell=True)

    if not args.tmp :
        subprocess.call("rm -rf "+ args.o + args.output_alignment_name, shell="True")

## branches
    print("Rscript scripts/calculate_mldist.R " + args.o + "file_fastRFS.single " + args.o)
    subprocess.call("Rscript scripts/calculate_mldist.R " + args.o + "file_fastRFS.single " + args.o,shell=True)
    files=glob.glob(args.o + "*.treefilematrix_ok")
    samples=[x.replace(args.o + ".treefilematrix_ok","") for x in files]
    l = me.get_all_length(samples)
    files, lengths = me.verif_matrix(files, l)
    print(len(lengths))
    me.get_all_matrix_with_alignement_length(files, lengths, args.o + args.name_res)

    subprocess.call("bin/erable -i " + args.o + args.name_res + " -t fastRFS_clean", shell=True)
    subprocess.call("sed -e 's,:-[0-9\.]\+,:0.0,g' " + args.o + args.name_res + ".lengths.nwk > " + args.o + args.name_res + ".tree",shell=True) #Suppr longueur de branche négative

    if not args.tmp :
        tmp_outfile = args.o + "*aligned.fasta*" + " " + args.o + "file_fastRFS"
        subprocess.call("rm " + tmp_outfile, shell="True")
