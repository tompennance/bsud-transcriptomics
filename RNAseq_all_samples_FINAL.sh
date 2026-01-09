#!/usr/bin/env bash

sample_folder=/nfs3/IB/Blouin_Lab/steinauer_lab/data/mrnaseq/usftp21.novogene.com/raw_data
work_folder=/nfs3/IB/Blouin_Lab/steinauer_lab/analysis/mrnaseq_analysis
genome_file=/nfs3/IB/Blouin_Lab/steinauer_lab/data/Bsud111_ccs_assembly.fasta
annotation_gff=/nfs3/IB/Blouin_Lab/steinauer_lab/data/bs_isoseq/PolyA_Annotation/Bsudanica_annotation_polyA_v2024.gtf

# change sample_files=XXX below to the files you want to run, or leave as the diretory (i.e. sample folder) if wanting to run all.

sample_files=$(ls -d $sample_folder"/"*/ | sed "s/\/$//" | sed "s/.*\///" | grep ^K)
#sample_files=B8
Illumina_adapters=$work_folder"/"TruSeq3-PE-2.fa
threads=30

TRIM=TRUE
HISAT_Host_With_Parasite=TRUE
HISAT_Parasite=TRUE
Parasite_FILTER=TRUE
HISAT_Host_Without_Parasite=FALSE
COUNT=TRUE

hisat_align () {
  echo "Running Hisat2 with:"
  echo "function_descriptor: "$function_descriptor
  echo "function_genome:" $function_genome
  echo "function_gff: "$function_gff
  echo "function_read_descriptor: "$function_read_descriptor
  echo "function_read_location: "$function_read_location
  echo "function_with_annotation: "$function_with_annotation
  echo ""
  echo ""

  rm -r $work_folder"/HISAT2_"$function_descriptor
  mkdir $work_folder"/HISAT2_"$function_descriptor
  mkdir $work_folder"/HISAT2_"$function_descriptor"/Reference"
  mkdir $work_folder"/HISAT2_"$function_descriptor"/Stats"

  if [ "$function_with_annotation" == "TRUE" ]
  then
    echo "HISAT Reference"
    #gffread -T -o $work_folder"/HISAT2_"$function_descriptor"/Reference/anot.gtf" $function_gff
    hisat2_extract_splice_sites.py $function_gff > $work_folder"/HISAT2_"$function_descriptor"/Reference/Splice_sites.ss"
    hisat2_extract_exons.py  $function_gff > $work_folder"/HISAT2_"$function_descriptor"/Reference/Exon_sites.exon"

    echo "HISAT Index"
    hisat2-build -q -p $threads $function_genome --ss $work_folder"/HISAT2_"$function_descriptor"/Reference/Splice_sites.ss" --exon $work_folder"/HISAT2_"$function_descriptor"/Reference/Exon_sites.exon" $work_folder"/HISAT2_"$function_descriptor"/Reference/Genome.ref"
  else
    echo "HISAT Index"
    hisat2-build -q -p $threads $function_genome $work_folder"/HISAT2_"$function_descriptor"/Reference/Genome.ref"
  fi

  for sample in $sample_files
  do
    echo "HISAT Alignment.. $sample"
    hisat2 -p $threads -1 $function_read_location"/"$sample"_"$function_read_descriptor"_1P.gz" -2 $function_read_location"/"$sample"_"$function_read_descriptor"_2P.gz" -x $work_folder"/HISAT2_"$function_descriptor"/Reference/Genome.ref" -S $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.sam" 2> $work_folder"/HISAT2_"$function_descriptor"/Stats/"$sample".Stats.txt"

    echo "Sorting bam file.. $sample"
    samtools sort -O BAM -o $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.bam" $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.sam"
    echo "Indexing bam file.. $sample"
    samtools index $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.bam"

    echo "$sample BAM file created and sorted"
    rm $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.sam"
  done
}


if [ "$TRIM" == TRUE ]
then
  rm -r $work_folder"/trimmed_reads"
  mkdir $work_folder"/trimmed_reads"

  for sample in $sample_files
  do
    echo $sample_folder"/"$sample"/"$sample"_1.fq.gz"
    echo $sample_folder"/"$sample"/"$sample"_2.fq.gz"
    trimmomatic.sh PE $sample_folder"/"$sample"/"$sample"_1.fq.gz" $sample_folder"/"$sample"/"$sample"_2.fq.gz" -baseout $work_folder"/trimmed_reads/"$sample'_trimmed_read.gz' -threads $threads ILLUMINACLIP:$Illumina_adapters":2:30:10" SLIDINGWINDOW:5:20 MINLEN:25
    echo "Done Trimming: "$sample
  done
fi


if [ "$HISAT_Host_With_Parasite" == TRUE ]
then
  function_descriptor=With_Parasite
  function_genome=/nfs3/IB/Blouin_Lab/steinauer_lab/data/Bsud111_ccs_assembly.fasta
  function_gff=/nfs3/IB/Blouin_Lab/steinauer_lab/data/bs_isoseq/PolyA_Annotation/Bsudanica_annotation_polyA_v2024.gtf
  function_read_descriptor=trimmed_read
  function_read_location=$work_folder"/trimmed_reads"
  function_with_annotation=TRUE

  hisat_align
fi


if [ "$HISAT_Parasite" == TRUE ]
then
  function_descriptor=Parasite
  function_genome=/nfs3/IB/Blouin_Lab/steinauer_lab/data/Smansoni_genome/schistosoma_mansoni.PRJEA36577.WBPS17.genomic.fa
  function_gff=/nfs3/IB/Blouin_Lab/steinauer_lab/data/Smansoni_genome/schistosoma_mansoni.PRJEA36577.WBPS17.annotations.gff3
  function_read_descriptor=trimmed_read
  function_read_location=$work_folder"/trimmed_reads"
  function_with_annotation=TRUE

  hisat_align

fi


if [ "$Parasite_FILTER" == TRUE ]
then
  rm -r $work_folder"/FILTERED_READS"
  mkdir $work_folder"/FILTERED_READS"
  function_descriptor=Parasite

  for sample in $sample_files
  do
    samtools view -b -f 4 $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.bam" > $work_folder"/FILTERED_READS/"$sample"_filtered_trimmed_read.bam"
    bedtools bamtofastq -i $work_folder"/FILTERED_READS/"$sample"_filtered_trimmed_read.bam" -fq $work_folder"/FILTERED_READS/"$sample"_filtered_trimmed_read"'_1P' -fq2 $work_folder"/FILTERED_READS/"$sample"_filtered_trimmed_read"'_2P'

    gzip $work_folder"/FILTERED_READS/"$sample"_filtered_trimmed_read"'_1P'
    gzip $work_folder"/FILTERED_READS/"$sample"_filtered_trimmed_read"'_2P'
    echo "Parasite reads removed. bam and fastq.gz generated"
  done
fi

if [ "$HISAT_Host_Without_Parasite" == TRUE ]
then
  function_descriptor=Without_Parasite
  function_genome=/nfs3/IB/Blouin_Lab/steinauer_lab/data/Bsud111_ccs_assembly.fasta
  function_gff=/nfs3/IB/Blouin_Lab/steinauer_lab/data/bs_isoseq/PolyA_Annotation/Bsudanica_annotation_polyA.gtf
  function_read_descriptor=filtered_trimmed_read
  function_read_location=$work_folder"/FILTERED_READS"
  function_with_annotation=TRUE

  for sample in $sample_files
  do
    echo "HISAT Alignment.. $sample"
    hisat2 -p $threads -1 $function_read_location"/"$sample"_"$function_read_descriptor"_1P.gz" -2 $function_read_location"/"$sample"_"$function_read_descriptor"_2P.gz" -x $work_folder"/HISAT2_"$function_descriptor"/Reference/Genome.ref" -S $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.sam" 2> $work_folder"/HISAT2_"$function_descriptor"/Stats/"$sample".Stats.txt"

    echo "Sorting bam file.. $sample"
    samtools sort -O BAM -o $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.bam" $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.sam"
    echo "Indexing bam file.. $sample"
    samtools index $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.bam"

    echo "$sample BAM file created and sorted"
    rm $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.sam"
  done

  echo "Host without parasite alignment complete"

fi


if [ "$COUNT" == TRUE ]
then
  #rm -r $work_folder"/Counts"
  #mkdir $work_folder"/Counts"
  function_descriptor=With_Parasite
  for sample in $sample_files
  do
    if [ ! -f $work_folder"/Counts/"$sample"_reads.count" ]
    then
      echo "Counting: "$sample
      htseq-count -r pos --stranded=yes --idattr=gene_id --format=bam --type=exon -q $work_folder"/HISAT2_"$function_descriptor"/"$sample"_alignment.bam" $annotation_gff > $work_folder"/Counts/"$sample"_reads.count"
      echo "DONE WITH THE SAMPLE!!"
    else
      echo Sample Repeated. Moving on.
    fi
  done

  count_da_genes=$(awk -F "\t" '{print $1}' $work_folder"/Counts/"$sample"_reads.count" | grep -v "__" )

  echo "GenID "$sample_files | tr " " "\t" >>  $work_folder"/Counts/All_reads.count"

  for gen in $count_da_genes
  do
    echo $gen
    record=$gen
    for sample in $sample_files
    do
      counting=$(grep -w -F $gen $work_folder"/Counts/"$sample"_reads.count" | awk -F "\t" '{print $2}')
      record=$(echo $record" "$counting)
    done
    echo $record | tr " " "\t" >>  $work_folder"/Counts/All_reads.count"
  done
fi
echo "Done with everything!!!"
