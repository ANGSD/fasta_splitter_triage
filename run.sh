#/bin/bash

#cut -f6,25,39  /maps/projects/caeg/data/db/aeDNA-refs/resources/20230719/ncbi/metadata/wgs.tsv |bgzip -c -@8 >wgs_f6_25_39.tsv

#touch seqs.fai.gz
#for i in /maps/projects/caeg/data/db/aeDNA-refs/temp/filter/nt/v5/generic/nt.fas /maps/projects/caeg/data/db/aeDNA-refs/temp/filter/custom/generic/artic.fas /maps/projects/caeg/data/db/aeDNA-refs/temp/filter/custom/generic/artic_other.fas /maps/projects/caeg/data/db/aeDNA-refs/temp/filter/custom/generic/norPlantCom.fas
#do
#    cat ${i}.fai |bgzip -c -@ 8 >>seqs.fai.gz    
#done

for i in  /maps/projects/caeg/data/db/aeDNA-refs/temp/filter/wgs/*/*.fas
do
    NAM=$(basename "${i}" .fas)
    HIT=$(grep "${NAM}" wgs_f6_25_39.tsv)
    if [ -z "$HIT" ]
    then
	echo "Problem finding $NAM in wgs_f6_25_39.tsv"
	exit 1;
    fi
    echo ${HIT}
done |tr " " "\t" >wgs.list

cut -f1-2 wgs.list |sort -k1 -n -S 10G|bgzip -c  >wgs_taxid_bp.txt.gz

time ./a.out makefai 1 -acc2taxid_flist seqs.acc2taxid.list  -meta_file seqs_old.fai -outname seqs_old.extended.fai.gz

#uncompress, extract size and taxid, sort by taxid, then take the sum of bp within each tqxid
gunzip -c seqs_old.extended.fai.gz |cut -f3,7|sort -k2 -n -S 10G |datamash -g 2 sum 1 |bgzip -c >seqs_taxid_bp.txt.gz
