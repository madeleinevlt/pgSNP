#!/bin/bash
#SBATCH --cpus-per-task=48
#SBATCH --job-name=only_typhimurium   # Job name
#SBATCH -o %x.%N.%j.out            # fichier où sera écrit la sortie standart STDOUT
#SBATCH -e %x.%N.%j.err            # fichier où sera écrit la sortie d'erreur STDERR
#SBATCH -p Research
source /global/conda/bin/activate
conda activate
conda activate superwgtree

id_path="/global/scratch/m.desousaviolante/paper_outbreak_typhimurium/pipeline_paper_sans_conta/id_sans_conta.txt"
paths=$(pwd)

python3 scripts/make_pangenome.py $id_path -name_ref reference.fasta -name_log log_reference_all_data.log -hit_min 500 -frag_min 500 -identity 95 ./

cat $id_path | while read name ; do
        echo "source /global/conda/bin/activate && conda activate && conda activate superwgtree && snippy --cpus 16 --outdir $name'_snip' --ref $paths'/reference.fasta' --R1 /global/bio/data/GAMeR_DB/SALMONELLA/$name/$name'_R1.fastq.gz' --R2 /global/bio/data/GAMeR_DB/SALMONELLA/$name/$name'_R2.fastq.gz'"  |  sed s/\'//g >> all_snippy_to_go
done


conda activate base

job=$(sarray -J all_snippy_to_go -o %x.%N.%j.out --dependency=$(squeue --noheader --format %i --name job_snippy2 | sed -r 's/_.*//' ) -p Research --cpus-per-task=16 all_snippy_to_go)


IFS=' ' read -ra ADDR <<< "$job"
num_job=${ADDR[3]}
#je lance tous les jobs avec sarray
a=$(squeue --job $num_job)
python3 /global/scratch/m.desousaviolante/prog/vcfs2SUPERalignement.py ./ -o result_with_non_codant.fasta --only_snp

python3 /global/scratch/m.desousaviolante/prog/merge_contig.py result_with_non_codant.fasta

#for file in *_aligned.fasta;
#do
#        echo "source /global/conda/bin/activate && conda activate && conda activate phylogeny && iqtree -s $file -m GTR+I+R+F -bb 1000 -nt 38 ">> all_iqtree
#done

#conda activate base
#job2=$(sarray -J all_iqtree -o %x.%N.%j.out --dependency=$(squeue --noheader --format %i --name job_iqtree | sed -r 's/_.*//' ) -p Research --cpus-per-task=16 all_iqtree)


#IFS=' ' read -ra ADDR <<< "$job2"
#num_job=${ADDR[3]}
#echo $num_job
#je lance tous les jobs avec sarray
#a=$(squeue --job $num_job)

#while [ ${#a} != 84 ]; do a=$(squeue --job $num_job); done


#cat *.treefile > all_iqtree_files

#conda activate superwgtree

#java -jar /global/scratch/m.desousaviolante/bin/astral/ASTRAL/astral.5.6.3.jar -i all_iqtree_files -o result_ASTRAL > tree_astral.log
#/global/scratch/m.desousaviolante/bin/ASTRID-linux -i all_iqtree_files -o file_astrid
#/global/scratch/m.desousaviolante/bin/FastRFS-linux/FastRFS -i all_iqtree_files -o file_fastRFS


files=$(ls *.treefilematrix_ok)
SAMPLES=`for i in $files;do printf $i | sed 's/.treefilematrix_ok/ /' ;done`
python3 /global/scratch/m.desousaviolante/prog/matrix_mldist_to_erable_input.py -i $files -a $SAMPLES -o distances

/global/scratch/m.desousaviolante/bin/erable/erable -i distances -t fastRFS_clean
