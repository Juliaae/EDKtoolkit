#!/bin/sh
#PBS -q core40
#PBS -l walltime=666:00:00,nodes=4:ppn=16,mem=60gb

##multi sample
###Usage: bash REDItools.sh .//bulk/PRJNA341833/1-output

##output pathway

for file in $1/*
do
input=$file
aa=$(echo "$input"| awk -F "/" '{print $NF}')
output=$2/${aa}
mkdir ${output}
echo "$input"
if [ -d $file ]; then 
		
	#Step1. Parallel Call Editing Using REDItools2
	SOURCE_BAM_FILE="${input}/*_star_RSEM.genome.sorted.bam"
	REFERENCE=".//reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	SIZE_FILE=".//reference/Homo_sapiens.GRCh38.dna.primary_assembly2.fa.fai"
	#BED_FILE="/pnas/zhangz_group/zhangysh/reference/hg38/GRCh38/Homo_sapiens.GRCh38.79.bed"
	
	NUM_CORES=4
	OUTPUT_FILE="${output}/parallel_table.txt.gz"
	TEMP_DIR="${output}/REDItools_temp/"
	COVERAGE_FILE="${output}/REDItools_coverage/*_star_RSEM.genome.sorted.cov"
	COVERAGE_DIR="${output}/REDItools_coverage/"
	if [ ! -f "${SOURCE_BAM_FILE}.bai" ]; then
		samtools index $SOURCE_BAM_FILE
	fi
	
	#if [ ! -f $COVERAGE_FILE ]; then
		#.//REDItools2-master/extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
	#fi
	
	if [ ! -f "${TEMP_DIR}merge-chronometer.txt" ]; then
		#mpirun -np $NUM_CORES /software/biosoft/software/python/python2.7_2018_12/bin/python .//REDItools2-master/src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -S -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR
		python ./REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 12 -i $SOURCE_BAM_FILE -f $REFERENCE -c 5,5 -q 30,30 -m 30,30 -O 5,5 -p -u -a 11-6 -l -v 1 -n 0.1 -k .//reference/nochr -o $output
	fi
	
	#Step2. Strand correction
	sed -i '1d' ${output}/DnaRna_*/outTable_*
	#gunzip --stdout ${output}/parallel_table.txt.gz > ${output}/parallel_table.txt
	python .//bulk/0_sh/Correct_strand.py -a .//reference/Homo_sapiens.sorted.GRCh38.99.gtf.gz -i ${output}/DnaRna_*/outTable_* -o ${output}/parallel_table_corrected.txt
	
	#Step3. Annotate Positions Using RepeatMasker and dbSNP Annotations
	REDI_OutFolder="${output}/REDItools_Annotation"
	mkdir $REDI_OutFolder
	
	cat ${output}/parallel_table_corrected.txt | sed 's/^/chr/g' > ${output}/parallel_table_corrected2.txt
	python ./REDItools-1.0.3/reditools/AnnotateTable.py -a .//reference/UCSC/rmsk/rmsk.gtf.gz -n rmsk -i ${output}/parallel_table_corrected2.txt -o ${REDI_OutFolder}/parallel_table.rmsk -u
	python ./REDItools-1.0.3/reditools/AnnotateTable.py -a .//reference/dbSNP/dbSNP.151.snp.gtf_sorted.gtf.gz -n snp151 -i ${REDI_OutFolder}/parallel_table.rmsk -o ${REDI_OutFolder}/parallel_table.rmsk.snp -u
	
	#Step4. Select Positions Separately
	python ./REDItools-1.0.3/reditools/selectPositions.py -i ${REDI_OutFolder}/parallel_table.rmsk.snp -c 5 -v 1 -f 0.0 -o ${REDI_OutFolder}/parallel_table.rmsk.snp.sel1
	python ./REDItools-1.0.3/reditools/selectPositions.py -i ${REDI_OutFolder}/parallel_table.rmsk.snp -c 10 -v 3 -f 0.1 -o ${REDI_OutFolder}/parallel_table.rmsk.snp.sel2
	
	##Step4.1. Select ALU Sites From the First Set of Positions and except SNP
	awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-" && $8!="-" && $13=="-") print}' ${REDI_OutFolder}/parallel_table.rmsk.snp.sel1 > ${REDI_OutFolder}/parallel_table.rmsk.snp.alu
	
	##Step4.2. Select REP NON ALU Sites From the Second Set of Positions, Excluding Sites in Simple Repeats or Low Complexity Regions
	awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $9>=0.1) print}' ${REDI_OutFolder}/parallel_table.rmsk.snp.sel2 > ${REDI_OutFolder}/parallel_table.rmsk.snp.nonalu
	
	##Step4.3. Select NON REP Sites From the Second Set of Positions
	awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $9>=0.1) print}' ${REDI_OutFolder}/parallel_table.rmsk.snp.sel2 > ${REDI_OutFolder}/parallel_table.rmsk.snp.nonrep
	
	
	#Step5. Annotate ALU, REP NON ALU and NON REP Sites Using Known Editing Events From REDIportal
	python ./REDItools-1.0.3/reditools/AnnotateTable.py -a .//reference/REDIportal/hg38.REDIportal_sorted.gtf.gz -n ed -c 1 -i ${REDI_OutFolder}/parallel_table.rmsk.snp.alu -o ${REDI_OutFolder}/parallel_table.rmsk.snp.alu.ed -u
	python ./REDItools-1.0.3/reditools/AnnotateTable.py -a .//reference/REDIportal/hg38.REDIportal_sorted.gtf.gz -n ed -c 1 -i ${REDI_OutFolder}/parallel_table.rmsk.snp.nonalu -o ${REDI_OutFolder}/parallel_table.rmsk.snp.nonalu.ed -u
	python ./REDItools-1.0.3/reditools/AnnotateTable.py -a .//reference/REDIportal/hg38.REDIportal_sorted.gtf.gz -n ed -c 1 -i ${REDI_OutFolder}/parallel_table.rmsk.snp.nonrep -o ${REDI_OutFolder}/parallel_table.rmsk.snp.nonrep.ed -u
	
	mv ${REDI_OutFolder}/*.rmsk.snp.alu.ed ${REDI_OutFolder}/alu
	mv ${REDI_OutFolder}/*.rmsk.snp.nonalu.ed ${REDI_OutFolder}/nonalu
	mv ${REDI_OutFolder}/*.rmsk.snp.nonrep.ed ${REDI_OutFolder}/nonrep
	cat ${REDI_OutFolder}/alu ${REDI_OutFolder}/nonalu ${REDI_OutFolder}/nonrep > ${REDI_OutFolder}/alu-nonalu-nonrep
	
	#Step6. Some REDIportal Editing Sites
	awk 'FS="\t" {if ($19!="-") print}' ${REDI_OutFolder}/alu-nonalu-nonrep > ${REDI_OutFolder}/knownEditing
	
	#Step7.1. Convert Editing Candidates in REP NON ALU and NON REP Sites in GFF Format for Further Filtering
	cat ${REDI_OutFolder}/nonalu ${REDI_OutFolder}/nonrep > ${REDI_OutFolder}/nonalu-nonrep
	awk 'FS="\t" {if ($19!="ed") print}' ${REDI_OutFolder}/nonalu-nonrep > ${REDI_OutFolder}/pos.txt
	sed -i 's/^chr//g' ${REDI_OutFolder}/pos.txt
	python ./REDItools-1.0.3/reditools/TableToGFF.py -i ${REDI_OutFolder}/pos.txt -s -o ${REDI_OutFolder}/pos.gff
	
	#Step7.2. Convert Editing Candidates in ALU Sites in GFF Format for Further Filtering
	awk 'FS="\t" {if ($19!="ed") print}' ${REDI_OutFolder}/alu > ${REDI_OutFolder}/posalu.txt
	sed -i 's/^chr//g' ${REDI_OutFolder}/posalu.txt
	python ./REDItools-1.0.3/reditools/TableToGFF.py -i ${REDI_OutFolder}/posalu.txt -s -o ${REDI_OutFolder}/posalu.gff
	
	#Step8. Launch REDItoolDnaRna.py on ALU Sites Using Stringent Criteria to Recover Potential Editing Candidates
	rm -r "${REDI_OutFolder}/firstalu"
	rm -r "${REDI_OutFolder}/first"
	rm -r "${REDI_OutFolder}/second"
	rm "${REDI_OutFolder}/pos.sorted.gff.gz"
	rm "${REDI_OutFolder}/posalu.sorted.gff.gz"
	rm "${REDI_OutFolder}/pos.sorted.gff.gz.tbi"
	rm "${REDI_OutFolder}/posalu.sorted.gff.gz.tbi"
	
	python ./REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 12 -i $SOURCE_BAM_FILE -f $REFERENCE -c 5,5 -q 30,30 -m 30,30 -O 5,5 -p -u -a 11-6 -l -v 1 -n 0.0 -e -T ${REDI_OutFolder}/posalu.sorted.gff -w .//reference/UCSC/mysplicesites_chrDel.ss -k .//reference/nochr -R -o ${REDI_OutFolder}/firstalu
	
	#Step9. Launch REDItoolDnaRna.py on REP NON ALU and NON REP Sites Using Stringent Criteria to Recover RNAseq Reads Harboring Reference Mismatches
	python ./REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 12 -i $SOURCE_BAM_FILE -f $REFERENCE -c 10,10 -q 30,30 -m 30,30 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T ${REDI_OutFolder}/pos.sorted.gff -w .//reference/UCSC/mysplicesites_chrDel.ss -k .//reference/nochr --reads -R --addP -o ${REDI_OutFolder}/first
	
	#Step10. Launch pblat on RNAseq Reads Harboring Reference Mismatches From Step 9 and Select Multimapping Reads
	/pnas/zhangz_group/zhutt/software/icebert-pblat-652d3b3/pblat -threads=12 -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 .//reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${REDI_OutFolder}/first/DnaRna_*/outReads_* ${REDI_OutFolder}/reads.psl
	./REDItools/accessory/readPsl.py ${REDI_OutFolder}/reads.psl ${REDI_OutFolder}/badreads.txt
	#Step11. Extract RNAseq Reads Harboring Reference Mismatches From Step 9 and Remove Duplicates
	sort -k1,1 -k2,2n -k3,3n ${REDI_OutFolder}/first/DnaRna_*/outPosReads_* | /software/biosoft/software/bedtools-2.25/bedtools2/bin/mergeBed > bed
	samtools view -@ 12 -L bed -h -b $SOURCE_BAM_FILE > ${REDI_OutFolder}/remove_duplicates_bed.bam
	samtools sort -@ 12 -n ${REDI_OutFolder}/remove_duplicates_bed.bam -o ${REDI_OutFolder}/remove_duplicates_bed_ns.bam
	samtools fixmate -@ 12 -m ${REDI_OutFolder}/remove_duplicates_bed_ns.bam ${REDI_OutFolder}/remove_duplicates_bed_ns_fx.bam
	samtools sort -@ 12 ${REDI_OutFolder}/remove_duplicates_bed_ns_fx.bam -o ${REDI_OutFolder}/remove_duplicates_bed_ns_fx_st.bam
	samtools markdup -r -@ 4 ${REDI_OutFolder}/remove_duplicates_bed_ns_fx_st.bam ${REDI_OutFolder}/remove_duplicates_bed_dedup.bam
	samtools index ${REDI_OutFolder}/remove_duplicates_bed_dedup.bam
	#Step12. Re-run REDItoolDnaRna.py on REP NON ALU and NON REP Sites Using Stringent Criteria, Deduplicated Reads and Mis-mapping Info
	python ./REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 12 -i ${REDI_OutFolder}/remove_duplicates_bed_dedup.bam -f $REFERENCE -c 10,10 -q 30,30 -m 30,30 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T ${REDI_OutFolder}/pos.sorted.gff.gz -w .//reference/UCSC/mysplicesites_chrDel.ss -k .//reference/nochr -R -b ${REDI_OutFolder}/badreads.txt -o ${REDI_OutFolder}/second
	
	#Step13. Collect filtered ALU, REP NON ALU and NON REP sites
	cd $REDI_OutFolder
	sed -i 's/chr//g' knownEditing
	sed -i '1d' knownEditing
	python .//REDItools2-master/NPscripts/collect_editing_candidates.py 
	sort -k1,1 -k2,2n editing.txt > editing_sorted.txt
	#Gene Symbol Annotation of Editing Sites
	python ./REDItools-1.0.3/reditools/Geneanno.py -a .//reference/Homo_sapiens.sorted.GRCh38.99.gtf.gz -i editing_sorted.txt -o result_anno.txt -c 1,2,3 -n Gene -S
	#Step14. Inspect the distribution of editing candidates to look at A-to-I enrichment
	python .//REDItools2-master/NPscripts/get_Statistics.py

fi
done