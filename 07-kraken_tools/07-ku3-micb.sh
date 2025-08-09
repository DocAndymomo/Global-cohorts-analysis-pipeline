#!/bin/bash

# Define your variables here
db_ku="/mnt/DATA4/DB_3/k2_database/01-Microbial2023"
folder_rcr="/mnt/DATA1/GEO-1/99-test/08-rcr"
folder_r2="/mnt/DATA1/GEO-1/99-test/09-ku3"
kt_code_path="/root/1-Biosoft/KrakenTools-1.2"

# Create output directory if it doesn't exist
mkdir -p $folder_r2

# Continue with the rest of your script using sorted input
ls "${folder_rcr}"/*_rcr_1.fastq | sort | while IFS= read -r filepath; do
  filename=$(basename "$filepath")
  # Rest of your script...
  sample_id="${filename%_rcr_1.fastq.gz}" 
  kraken_report="${folder_r2}/${sample_id}.kraken.report.txt" 
  kraken_out="${folder_r2}/${sample_id}.kraken.output.txt" 
  faq_rcr1="${folder_mr}/${sample_id}_rcr_1.fastq" 
  faq_rcr2="${folder_mr}/${sample_id}_rcr_2.fastq" 

  if [ -f "${folder_r2}/${sample_id}.kraken.report.txt" ]; then
    echo "Skipping sample ${sample_id} as ku2 output already exists."
    continue
  fi



  krakenuniq \
        --db $db_ku \
        --threads 40 \
        --preload \
        --report-file $kraken_report \
        --paired $faq_rcr1 $faq_rcr2 &

  mpa_output="${folder_r2}/${sample_id}.kraken.mpa.txt" 
  python "${kt_code_path}/kreport2mpa.py" -r $kraken_report -o $mpa_output 
  krona_output="${folder_r2}/${sample_id}.krona.txt" 
  python "${kt_code_path}/kreport2krona.py" -r $kraken_report -o $krona_output 

done
###########################################################


#!/bin/bash

# Define your variables here
db_ku="/home/rstudio/SAHMI/k2_pluspf_20230605"
folder_rcr="/mnt/DATA1/TCGA/07-rcr"
folder_r2="/mnt/DATA1/TCGA/08-ku2"
kt_code_path="/root/1-Biosoft/KrakenTools-1.2"

# Create output directory if it doesn't exist
mkdir -p $folder_r2

# Continue with the rest of your script using sorted input
ls "${folder_rcr}"/*_rcr_1.fastq.gz | sort -r | while IFS= read -r filepath; do
  filename=$(basename "$filepath")
  # Rest of your script...
  sample_id="${filename%_rcr_1.fastq.gz}" 
  kraken_report="${folder_r2}/${sample_id}.kraken.report.txt" 
  kraken_out="${folder_r2}/${sample_id}.kraken.output.txt" 
  faq_rcr1="${folder_mr}/${sample_id}_rcr_1.fastq.gz" 
  faq_rcr2="${folder_mr}/${sample_id}_rcr_2.fastq.gz" 

  if [ -f "${folder_r2}/${sample_id}.kraken.report.txt" ]; then
    echo "Skipping sample ${sample_id} as ku2 output already exists."
    continue
  fi

  /root/1-Biosoft/kraken2-2.1.3/kraken2 --db "$db_ku" --threads 32 --paired --report-minimizer-data --report $kraken_report --use-names --output $kraken_out --paired $faq_rcr1 $faq_rcr2 
  mpa_output="${folder_r2}/${sample_id}.kraken.mpa.txt" 
  python "${kt_code_path}/kreport2mpa.py" -r $kraken_report -o $mpa_output 
  krona_output="${folder_r2}/${sample_id}.krona.txt" 
  python "${kt_code_path}/kreport2krona.py" -r $kraken_report -o $krona_output 

done



# 合并数据到mpa
############################################################
bash
conda activate rna
folder3="/mnt/DATA1/TCGA/07-ku2"
folder4="/mnt/DATA1/TCGA/08-otu-mpa"
kt_code_path="/root/1-Biosoft/KrakenTools-1.2"

python ${kt_code_path}/combine_mpa.py \
      -i $folder3/*.mpa.txt \
      -o $folder4/0-combine_mpa.txt



for file in "${folder3}"/*.mpa.txt; do
    sample_id=$(echo $file | sed 's/\.kraken\.mpa\.txt//')
    echo $sample_id >> "${folder4}"/0-sample_id.txt
done

sort -u "${folder4}"/0-sample_id.txt -o "${folder4}"/0-sample_id.txt

############################################################


# 合并数据到otu
############################################################
bash

folder3="/mnt/DATA1/TCGA/08-ku2"
folder4="/mnt/DATA1/TCGA/09-otu-mpa"
kt_code_path="/root/1-Biosoft/itm_helper"

python3 ${kt_code_path}/kraken2otu.py --extension kraken.report.txt --inputfolder $folder3 --level p --outdir $folder4
python3 ${kt_code_path}/kraken2otu.py --extension kraken.report.txt --inputfolder $folder3 --level f --outdir $folder4
python3 ${kt_code_path}/kraken2otu.py --extension kraken.report.txt --inputfolder $folder3 --level g --outdir $folder4
python3 ${kt_code_path}/kraken2otu.py --extension kraken.report.txt --inputfolder $folder3 --level s --outdir $folder4

mv "${folder3}"/otu_table* $folder4
#############################################################
