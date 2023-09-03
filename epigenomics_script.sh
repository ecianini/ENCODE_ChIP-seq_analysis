#### EPIGENOMICS PROJECT#####

#1) 
# After the download, we want to obtain the percentage of mapped reads, multiple-mapped reads (unf mapped - flt mapped)
# A good percentage of mapped reads is  = > 75%
#FLAGSTAT
samtools flagstat encode/unf1.bam > encode/stats/unf1.txt #28689347 (96.85%) 
samtools flagstat encode/flt1.bam > encode/stats/flt1.txt # 23039346 (100%) multiple mapped = 19.69% -> 23040415 80.31% unique
samtools flagstat encode/unf2.bam > encode/stats/unf2.txt # 40983039 (96.79%)
samtools flagstat encode/flt2.bam > encode/stats/flt2.txt # 32323595 (100%)  multiple mapped = 21.13% -> 32323323 78.87% unique 

#multiple mapped reads:  
# unf1 - flt1 = 28689347 - 23039346 = 5650001 
# unf2 - flt2 = 40983039-32323595 = 8659444 

#inspect the control file 
samtools flagstat encode/ctrl.unf.bam > encode/stats/ctrl.unf.txt #56538503 total 46184772(81.69%)
samtools flagstat encode/ctrl.flt.bam > encode/stats/ctrl.flt.txt #30799798 (100.00%)


#2) 
#PEAK CALLING
macs2 callpeak -t encode/unf1.bam -c encode/ctrl.unf.bam -g hs -q 0.01 --outdir peaks -n unf1 2> peaks/logs/peak.unf1.log 
macs2 callpeak -t encode/unf2.bam -c encode/ctrl.unf.bam -g hs -q 0.01 --outdir peaks -n unf2 2> peaks/logs/peak.unf2.log 
macs2 callpeak -t encode/flt1.bam -c encode/ctrl.flt.bam -g hs -q 0.01 --outdir peaks -n flt1 2> peaks/logs/peak.flt1.log 
macs2 callpeak -t encode/flt2.bam -c encode/ctrl.flt.bam -g hs -q 0.01 --outdir peaks -n flt2 2> peaks/logs/peak.flt2.log 

macs2 callpeak -t encode/unf1.bam encode/unf2.bam -c encode/ctrl.unf.bam -g hs -q 0.01 --outdir peaks -n mrg.unf 2> peaks/logs/peak.mrg.unf.log 
macs2 callpeak -t encode/flt1.bam encode/flt2.bam -c encode/ctrl.flt.bam -g hs -q 0.01 --outdir peaks -n mrg.flt 2> peaks/logs/peak.mrg.flt.log 

#results peak calling RR < 0.2 
#unf1 = RR = 0.05 trt e 0.05 ctrl, fragment size = 136 bps (predicted)
#unf2 = RR = 0.07 trt e 0.05 ctrl, fragment size = 135 bps (predicted)
#flt1 = RR = 0.00 trt e 0.00 ctrl, fragment size = 145 bps (predicted)
#flt2 = RR = 0.00 trt e 0.00 ctrl, fragment size = 142 bps (predicted)

#mrg.unf = RR 0.07 trt e 0.05 ctrl, fragment size = 137 bps (predicted)
#mrg.flt = RR 0.01 trt e 0.00 ctrl, fragment size = 144 bps (predicted)

# For clearer results, we move the files into another directory
mv peaks/*.xls other/
mv peaks/*model.r other/


# 3)
# Low complexity regions (i.e, BLACKLISTED regions) must be removed from the peaks

#for the narrowpeaks
bedtools intersect -a peaks/unf1_peaks.narrowPeak -b encode/blaclist.bed -v > peaks/unf1.clean_peaks.narrowPeak ; mv peaks/unf1.clean_peaks.narrowPeak peaks/unf1_peaks.narrowPeak
bedtools intersect -a peaks/unf2_peaks.narrowPeak -b encode/blaclist.bed -v > peaks/unf2.clean_peaks.narrowPeak ; mv peaks/unf2.clean_peaks.narrowPeak peaks/unf2_peaks.narrowPeak
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b encode/blaclist.bed -v > peaks/flt1.clean_peaks.narrowPeak ; mv peaks/flt1.clean_peaks.narrowPeak peaks/flt1_peaks.narrowPeak
bedtools intersect -a peaks/flt2_peaks.narrowPeak -b encode/blaclist.bed -v > peaks/flt2.clean_peaks.narrowPeak ; mv peaks/flt2.clean_peaks.narrowPeak peaks/flt2_peaks.narrowPeak

bedtools intersect -a peaks/mrg.unf_peaks.narrowPeak -b encode/blaclist.bed -v > peaks/mrg.unf.clean_peaks.narrowPeak ; mv peaks/mrg.unf.clean_peaks.narrowPeak peaks/mrg.unf_peaks.narrowPeak
bedtools intersect -a peaks/mrg.flt_peaks.narrowPeak -b encode/blaclist.bed -v > peaks/mrg.flt.clean_peaks.narrowPeak ; mv peaks/mrg.flt.clean_peaks.narrowPeak peaks/mrg.flt_peaks.narrowPeak

# for the summits
bedtools intersect -a peaks/unf1_summits.bed -b encode/blaclist.bed -v > peaks/unf1.clean_summits.bed ; mv peaks/unf1.clean_summits.bed peaks/unf1_summits.bed
bedtools intersect -a peaks/unf2_summits.bed -b encode/blaclist.bed -v > peaks/unf2.clean_summits.bed ; mv peaks/unf2.clean_summits.bed peaks/unf2_summits.bed
bedtools intersect -a peaks/flt1_summits.bed -b encode/blaclist.bed -v > peaks/flt1.clean_summits.bed ; mv peaks/flt1.clean_summits.bed peaks/flt1_summits.bed
bedtools intersect -a peaks/flt2_summits.bed -b encode/blaclist.bed -v > peaks/flt2.clean_summits.bed ; mv peaks/flt2.clean_summits.bed peaks/flt2_summits.bed

bedtools intersect -a peaks/mrg.unf_summits.bed -b encode/blaclist.bed -v > peaks/mrg.unf.clean_summits.bed ; mv peaks/mrg.unf.clean_summits.bed peaks/mrg.unf_summits.bed
bedtools intersect -a peaks/mrg.flt_summits.bed -b encode/blaclist.bed -v > peaks/mrg.flt.clean_summits.bed ; mv peaks/mrg.flt.clean_summits.bed peaks/mrg.flt_summits.bed


#4) 
#PEAKS COUNTS
#narrow peaks
wc -l peaks/unf1_peaks.narrowPeak # 16995
wc -l peaks/unf2_peaks.narrowPeak # 23430	
wc -l peaks/flt1_peaks.narrowPeak # 15426
wc -l peaks/flt2_peaks.narrowPeak # 19698
wc -l peaks/mrg.unf_peaks.narrowPeak # 21572
wc -l peaks/mrg.flt_peaks.narrowPeak # 16938

#summits
wc -l peaks/unf1_summits.bed #16996
wc -l peaks/unf2_summits.bed #23432
wc -l peaks/flt1_summits.bed #15426
wc -l peaks/flt2_summits.bed #19699
wc -l peaks/mrg.unf_summits.bed #21574
wc -l peaks/mrg.flt_summits.bed #16939

#IDR regions
wc -l encode/IDR12.bed #22236
wc -l encode/IDR1.bed #13039
wc -l encode/IDR2.bed #17310


# 5) OVERLAP
# we decide to compare the replicates and to consider as "overlapping summits" those having a distance <= 100 bp
bedtools closest -a peaks/unf1_summits.bed -b peaks/flt1_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #14536
bedtools closest -a peaks/unf1_summits.bed -b peaks/unf2_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #14668
bedtools closest -a peaks/unf2_summits.bed -b peaks/flt2_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #19080
bedtools closest -a peaks/flt1_summits.bed -b peaks/flt2_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #13218

bedtools closest -a peaks/unf1_summits.bed -b peaks/mrg.unf_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #15711
bedtools closest -a peaks/unf2_summits.bed -b peaks/mrg.unf_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #19896
bedtools closest -a peaks/flt1_summits.bed -b peaks/mrg.flt_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #13794
bedtools closest -a peaks/flt2_summits.bed -b peaks/mrg.flt_summits.bed -d | awk '$11<=100{c++} END{print c+0}' #16115


#6) 
#Create the IDR files with the filtered bam files using the ENCODE parameters
#call less stringent peaks
macs2 callpeak -t encode/flt1.bam -c encode/ctrl.flt.bam -g hs -n peaks/myIDR1peak -q 0.1  2> peaks/logs/myIDR1peak.log
macs2 callpeak -t encode/flt2.bam -c encode/ctrl.flt.bam -g hs -n peaks/myIDR2peak -q 0.1  2> peaks/logs/myIDR2peak.log

# move what is not a peak in directory other/
mv peaks/*.xls other/
mv peaks/*model.r other/

# Sort peaks by -log10(q-value)
sort -nk9 peaks/myIDR1peak_peaks.narrowPeak > peaks/myIDR1peak2_peaks.narrowPeak ; mv peaks/myIDR1peak2_peaks.narrowPeak peaks/myIDR1peak_peaks.narrowPeak
sort -nk9 peaks/myIDR2peak_peaks.narrowPeak > peaks/myIDR2peak2_peaks.narrowPeak ; mv peaks/myIDR2peak2_peaks.narrowPeak peaks/myIDR2peak_peaks.narrowPeak

# create our IDR file
# run the command on the server because it has the command idr installed
idr -s peaks/myIDR1peak_peaks.narrowPeak peaks/myIDR2peak_peaks.narrowPeak -i 0.2 -o myIDR/myIDR12.bed 2> myIDR/myIDR12.log
# Number of peaks passing IDR cutoff of 0.2 - 13608/17873 (76.1%) <- info in the log file

wc -l myIDR/myIDR12.bed #13608 
# number of peaks from encode 22236 


#7) 
#INTERSECTION RESULTS VS ENCODE IDR files
# we take into consideration only the mrg.file compared to the narrow peaks of ENCODE
bedtools intersect -a peaks/mrg.flt_peaks.narrowPeak -b encode/IDR12.bed | wc -l #15157
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b encode/IDR12.bed | wc -l #13875
bedtools intersect -a peaks/flt2_peaks.narrowPeak -b encode/IDR12.bed | wc -l #16221
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b peaks/flt2_peaks.narrowPeak | bedtools intersect -a - -b encode/IDR12.bed | wc -l #13139
bedtools intersect -a myIDR/myIDR12.bed -b encode/IDR12.bed | wc -l #12678

#intersect flt1 and flt2 
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b peaks/flt2_peaks.narrowPeak | wc -l #13834
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b peaks/flt2_peaks.narrowPeak | bedtools intersect -a - -b encode/IDR12.bed | wc -l #13139


# 8)
# extract the -log10(qvalues) for visualization purposes (boxplots)
bedtools intersect -a peaks/mrg.flt_peaks.narrowPeak -b encode/IDR12.bed | awk '{print $9}' > other/enr.mrg.IDR12.csv ; bedtools intersect -a peaks/mrg.flt_peaks.narrowPeak -b encode/IDR12.bed -v | awk '{print $9}' | paste other/enr.mrg.IDR12.csv - > other/enr.mrg.IDR12.final.csv ; mv other/enr.mrg.IDR12.final.csv other/enr.mrg.IDR12.csv
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b encode/IDR12.bed | awk '{print $9}' > other/enr.flt1.IDR12.csv ; bedtools intersect -a peaks/flt1_peaks.narrowPeak -b encode/IDR12.bed -v | awk '{print $9}' | paste other/enr.flt1.IDR12.csv - > other/enr.flt1.IDR12.final.csv ; mv other/enr.flt1.IDR12.final.csv other/enr.flt1.IDR12.csv
bedtools intersect -a peaks/flt2_peaks.narrowPeak -b encode/IDR12.bed | awk '{print $9}' > other/enr.flt2.IDR12.csv ; bedtools intersect -a peaks/flt2_peaks.narrowPeak -b encode/IDR12.bed -v | awk '{print $9}' | paste other/enr.flt2.IDR12.csv - > other/enr.flt2.IDR12.final.csv ; mv other/enr.flt2.IDR12.final.csv other/enr.flt2.IDR12.csv
bedtools intersect -a peaks/flt1_peaks.narrowPeak -b peaks/flt2_summits.bed > temp.bed ; bedtools intersect -a temp.bed -b encode/IDR12.bed | awk '{print $9}' > other/enr.fltINT12.IDR12.csv ; bedtools intersect -a temp.bed -b encode/IDR12.bed -v | awk '{print $9}' | paste other/enr.fltINT12.IDR12.csv - > other/enr.fltINT12.IDR12.final.csv ; mv other/enr.fltINT12.IDR12.final.csv other/enr.fltINT12.IDR12.csv; rm temp.bed
bedtools intersect -a myIDR/myIDR12.bed -b encode/IDR12.bed |  awk '{print $7}' > other/enr.myIDR12.csv ; bedtools intersect -a myIDR/myIDR12.bed -b encode/IDR12.bed -v | awk '{print $7}' | paste other/enr.myIDR12.csv - > other/enr.myIDR12.final.csv ; mv other/enr.myIDR12.final.csv other/enr.myIDR12.csv


# 9) 
# download the chromHMM file for the chromatine state of the entire genome in order to find the cromatine state with seqminer
bedtools intersect -a peaks/mrg.flt_summits.bed -b encode/chromHMM.bed -wa -wb | awk '{print $9}' > other/cromHMM_mrg.states.csv


#10) 
# crossing with ATAC-seq results
bedtools intersect -a peaks/mrg.flt_summits.bed -b encode/ATAC.bed -wa | sort | uniq | wc -l
# 16939 initial summits 
# Overlapping: 13418 -> 80% overlapping
# not overlapping: 3521 -> 20% not overlapping

#11) 
# Determine the TF' target genes with GREAT using the 30K nearest neighbour rule
#sort the file
sort -nk5 peaks/mrg.flt_summits.bed | awk '{print $1,$2,$3,$4}' | tail -n 5000 - > peaks/mrg.flt.GREAT_summits.bed


#12)
# As a final step we use PscanChIP to determine the enrichment motifs
