#!/bin/bash
#$ -N Bcl2Fastq
#$ -cwd

#Command line example:
# GenoStep1.sh targetpath metadata samplesheet

#variables:
#output_directory where files for analysis will be created
targetpath=$1

#metadata_file
#metadata="metadata_200325_NB501183.txt"
metadata=$2
#SampleSheet_byWell_TARGET.csv
samplesheet=$3

#load required modules

module load bcl2fastq/2.20.0.422
module load fastq-tools
module load fastq-pair

ulimit -n 4000

cd ${targetpath}

mkdir -p fastq_per_well

bcl2fastq -o ${targetpath}fastq_per_well --sample-sheet ${samplesheet} --no-lane-splitting

rm -rf fastq_per_well/Stats
rm -rf fastq_per_well/Reports
rm -rf fastq_per_well/Undetermined*

#list all available barcodes for PCR1 from metadata file

awk -F"\t" '!seen[$2, $3]++' ${metadata} > unique.barcodes.txt

module load fastq-tools

mkdir -p fastq_per_well_plate/

#unzip files before running fastq-grep
gunzip fastq_per_well/*.gz

#iterate through all barcode plates and demultiplex fastq files
sed 1d unique.barcodes.txt | while read line
do
   	plate=`awk -F"\t" '{print $3}' <<< "$line"`
	barcode=`awk -F"\t" '{print $2}' <<< "$line"`

	for file in fastq_per_well/*R1*.fastq
	do
		filename=`basename $file`
		fastq-grep -i "$barcode" $file > fastq_per_well_plate/${plate}_${filename}

	done

	for file in fastq_per_well/*R2*.fastq
	do
		filename=`basename $file`
		fastq-grep -i "$barcode" $file > fastq_per_well_plate/${plate}_${filename}

	done
done

#tidy-up
rm unique.barcodes.txt

#run fastq-pair to confirm that reads are properly paired

mkdir -p fastq_per_cellid/

cd fastq_per_well_plate/

for file1 in *R1*.fastq

do
        prefix=${file1%_R1*}

        file2=$(ls ${prefix}_*R2*.fastq)

        fastq_pair $file1 $file2
        
        mv ${file1}.paired.fq ../fastq_per_cellid/${file1}.paired.fq
        mv ${file2}.paired.fq ../fastq_per_cellid/${file2}.paired.fq
        
        rm -f ${file1}.single.fq
        rm -f ${file2}.single.fq

done

#Rename files with cell_id from metadata file

tail -c1 < "$metadata" | read -r _ || echo >> "$metadata"

cd ../fastq_per_cellid/
#zip files after running fastq-pair
gzip *.fq

sed 1d $metadata | while read line
do
	plate=$(echo $line | awk '{print $3}')
        well=$(echo $line | awk '{print $5}')
        cellid=$(echo $line | awk '{print $1}')

        targetFile1="${plate}_${well}_*R1*.gz"
        targetFile2="${plate}_${well}_*R2*.gz"

        mv ${targetFile1} "${cellid}_R1.fq.gz"
        mv ${targetFile2} "${cellid}_R2.fq.gz"

done

rm -rf Plate*

##Optional; remove fastq_per_well/ and fastq_per_well_plate/ directories:
#rm -rf ${targetpath}fastq_per_well
#rm -rf ${targetpath}fastq_per_well_plate
