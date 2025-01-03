#!/bin/bash


#Software used
#fasta_formatter Part of FASTX Toolkit 0.0.14
#VIBRANT v1.2.1
#bwa 0.7.17-r1188
#samtools 1.13
#bedtools v2.31.0
#R
#trim_galore 0.6.7
#samclip 0.4.0


#This script needs a reference genome file in fasta format and two Illumina paired-end fastq files per treatment. These should be trimmed for adapters with trim_galore or similar.

#Make symlinks to raw data files


#Create folders for various samples and symlinks to raw data
for i in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        mkdir $i $i/$i$time
        ln -s /home/tkn/work/EUCAST_mob19/coli_pilot/raw_data/$i$time/$i$time"_R1_val_1.fq" $i/$i$time/$i$time"_R1_val_1.fq"
        ln -s /home/tkn/work/EUCAST_mob19/coli_pilot/raw_data/$i$time/$i$time"_R2_val_2.fq" $i/$i$time/$i$time"_R2_val_2.fq"
    done
done

#Start looking for prophages - takes a while
nohup python3 /home/tkn/Programs/VIBRANT/VIBRANT_run.py -i /home/tkn/work/EUCAST_mob19/coli_pilot/MG1655_pOLA52.fasta -t 30 &

#Make a function to loop through samples to find peaks in mobilome coverage.



peak_calling () {
treatment=$1
time=$2

genomefile=/DATA_1/tkn/Mobilome_tool/U00096/U00096_3.fasta

#Make sure that linewidth is unlimited in genome file
fasta_formatter -i $genomefile -w0 -o genomefile.fna
samtools faidx genomefile.fna

mob_fwd=*_val_1.fq
mob_rev=*_val_2.fq


#First screening of highly covered mob regions
replicons=$(grep '>' genomefile.fna | sed 's/>//')

#Map mob data to genomefile
bwa-mem2 index genomefile.fna

bwa-mem2 mem -t 30 genomefile.fna $mob_fwd $mob_rev | samtools view -S -b --threads 30 | samtools sort --threads 30 -o mob_mapped.bam
#Only allow "proper pairs"
samtools view -@30 -f 0x2 -F 256 -F 2048 -b -h mob_mapped.bam > mob_mapped_proper.bam
#Get discordantly mapped reads
samtools view -@30 -F 1294 -b -h mob_mapped.bam > discordant.bam
mv mob_mapped_proper.bam mob_mapped.bam

samtools index mob_mapped.bam
samtools index discordant.bam

bedtools genomecov -d -ibam mob_mapped.bam > coverage.txt
bedtools genomecov -d -ibam discordant.bam > disc_coverage.txt
awk '{print $1,$2-1,$2,$3}' coverage.txt | sed 's/ /\t/g' > coverage.bed


#Split bam file into replicons
replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}' | sort -u)

rm -fr splitbam/* genome.bed genome.bed
mkdir splitbam
#samtools index mob_mapped.bam
for chrom in $replicons
do
    echo $chrom
    samtools view -bh mob_mapped.bam $chrom > splitbam/$chrom.bam
    samtools index splitbam/$chrom.bam
    rep_length=$(grep -A1 -w ">$chrom" genomefile.fna | grep -v '>' | awk '{print length($1)}')
    echo "$chrom $rep_length" >> genome.bed
done

sed -i 's/ /\t/' genome.bed

#bedtools makewindows -g genome.bed -w 10 > genome10.bed
bedtools makewindows -g genome.bed -w 1 > genome1.bed
bedtools makewindows -g genome.bed -w 100 > genome100.bed

replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}' | sort -u)
for chrom in $replicons
do
    rep_length=$(grep -A1 -w ">$chrom" genomefile.fna | grep -v '>' | awk '{print length($1)}')
    bedtools map -a genome100.bed -b coverage.bed -c 4 -o mean | grep "$chrom" > "$chrom"_treat_pileup_100int.bdg
done

cat *_treat_pileup_100int.bdg > all_100int.bed


#Get (soft) clipped reads only. These are not used for peak calling but for expanding called regions.samtools faidx genomefile.fna
samtools view -h mob_mapped.bam | samclip --ref genomefile.fna --invert --max 75 | samtools view -b >  samclip.bam
bedtools genomecov -d -ibam samclip.bam > clipped_coverage.txt
rm -rf high_clipped.bed
for repl in $replicons
do
    #Use an arbitrary cutoff of 3 clipped reads for a position. 
    grep -w "$repl" clipped_coverage.txt | sort -nk2,2 | awk -v cutoff=3 '{if ($3 > cutoff) print $1,$2,$2+1,$3,$4}' | sed -e 's/ /\t/g' -e 's/\t$//' > tempcov$repl
    #A region with clipped reads should be at least 10 bases long. Avoid spurious regions.
    bedtools merge -i tempcov$repl -c 4 -o mean -d 50 | awk -v min_size=10 '{if ($3-$2 > min_size) print $0"\t"$3-$2"\t"median_cov"\t"len}' >> high_clipped.bed
done

echo '
library(dplyr)
indat <- read.csv("NA.cov.bed", header=F, sep="\t")

#Remove sites with 20 times higher coverage than median - obviously interesting but interrupts with kernel density estimates for finding less obvious peaks. 
indat_filt <- indat %>%
    filter(V4 < (median(V4))+1*20)
    
indat_filt <- indat_filt %>%
    filter(V4 > (median(V4))+1*2)

repl <- unique(indat$V1)
outfile <- paste("antimode_treat_pileup_",repl, sep ="")

#Make the cutoff at the 97.5% quantile of coverage.
write(data.frame(quantile(indat_filt$V4, 0.975))[,1],outfile)
' > temp.R

replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}' | sort -u)
for repl in $replicons
do
    grep "$repl" coverage.bed > $repl.cov.bed
    sed "s/NA/$repl/" temp.R > temp2.R
    chmod +x temp2.R && Rscript ./temp2.R
done

replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}' | sort -u)
rm -rf high_cov_intervals.bed
for repl in $replicons
do
    echo $repl
    cutoff=$(cat antimode_treat_pileup_$repl)
    #Sometimes cutoff is extremely low. Coverage should be at least 5x
    min_cutoff=5
    if [ "$cutoff" -eq "$min_cutoff" ]; 
    then
    cutoff=5
    fi
    echo $cutoff
    grep -w "$repl" "$repl".cov.bed | sort -nk2,2 | awk -v cutoff=$cutoff '{if ($4 > cutoff) print $1,$2,$3,$4}' | sed 's/ /\t/g' > tempcov$repl
    bedtools merge -i tempcov$repl -c 4 -o mean -d 1000 | awk -v min_size=100 '{if ($3-$2 > min_size) print $0"\t"$3-$2"\t"median_cov"\t"len}' >> high_cov_intervals.bed
done


num_regions=$(wc -l high_cov_intervals.bed)
echo "There are $num_regions preliminary highly covered regions"

#If there is an abnormal number of regions passing coverage cutoff, the filter is too loose. Subset to 30 most abundant regions. From manual inspection of results from all samples, this will capture all regions of interest. This number should not be assumed to work for other genomes/samples.
if [ "$num_regions" > 30 ]; then
    sort -nk4,4 high_cov_intervals.bed | tail -n 30 > temp
    echo "Subsetting to the 30 most higly covered regions"
fi
mv temp high_cov_intervals.bed



#Merge high coverage regions with regions with high numbers of soft clipped reads (>3) - extends regions if there is lower coverage in the ends, but still clipped reads indicating circularity
replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}')
rm -rf replicon_coverage.txt
for repl in $replicons
do
    grep -w "$repl" "$repl".cov.bed | sort -nk4,4 | awk -v rep=$repl '{ a[i++]=$4; } END { print rep,a[int(i/2)]; }' >> replicon_coverage.txt
    grep -w "$repl" "$repl".cov.bed | sort -nk4,4 | awk -v rep=$repl '{total += $4} END {print rep,total/NR}' >> replicon_coverage_mean.txt
done

rm -rf high_cov_intervals_clip.bed
awk '{print $1"_"NR,$2,$3,$4,$5}' high_cov_intervals.bed | sed 's/ /\t/g'  > temp
intervals=$(awk '{print $1}' temp)

#Loop through intervals and extend start/end positions based on clipped reads. Takes a while.
for i in $intervals
do
    replicon=$(echo $i | sed 's/_[0-9]*$//')
    echo $i
    median_cov=$(grep -w "$replicon" "$replicon".cov.bed | sort -nk4,4 | awk ' { a[i++]=$4; } END { print a[int(i/2)]; }')
    if [[ "$median_cov" < 0.1 ]]
    then
        median_cov=0.1
    fi

    repl_length=$(grep -w -A1 "$replicon" genomefile.fna | grep -v '>' | awk '{print length($1)}')
    
    s=$(grep -w "$i" temp | awk '{print $2}')
    new_start=$(awk -v start_pos=$s '{if ($2 > start_pos-100 && $2 < start_pos+100) print $2}' high_clipped.bed | head -n1)
    if [ -z "$new_start" ]
    then
        new_start=$s
    else
        new_start=$new_start
    fi
    
    e=$(grep -w "$i" temp | awk '{print $3}')
    new_end=$(awk -v end_pos=$e '{if ($3 > end_pos-200 && $3 < end_pos+200) print $3}' high_clipped.bed | tail -n1)
    if [ -z "$new_end" ]
    then
        new_end=$e
    else
        new_end=$new_end
    fi
    grep -w "$i" temp | grep -w "$s" | grep -w "$e" | \
    awk -v repl=$replicon -v new_start=$new_start -v new_end=$new_end -v median_cov=$median_cov -v len=$repl_length '{print repl,new_start,new_end,$4,$5,median_cov,len}' | \
    sed 's/ /\t/g' >> high_cov_intervals_clip.bed
done
        

#Find the first antimode in the coverage distribution of regions passing the above filter.
#Remove regions that fall below the first antinode coverage
echo '
library("multimode")
library("dplyr")
high_cov = read.csv("high_cov_intervals_clip.bed", sep = "\t", header = F)
replicons=unique(high_cov$V1)
for (repl in replicons){
    #Subset to current replicon; filter out those who have extremely high increase in coverage (those are surely true positives). 
    temp=subset(high_cov, grepl(repl,high_cov$V1))
    temp = temp %>%
    filter(temp$V4/temp$V6 < 5)
    #Only perform antinode filtering, if there are more than 3 entries
    if (nrow(temp) >= 3) {
        print(nrow(temp))
        #Define the number of modes in the data. Use the median coverage of high regions as scaling factor.
        n_modes = nmodes(temp$V4/temp$V6,median(temp$V4))
        if(n_modes <= 2){
            n_modes = 3
            }
            if(n_modes > nrow(temp)) {
                n_modes = nrow(temp)
                }
    modes = locmodes(temp$V4/temp$V6, mod0=n_modes)
    outfile=paste("antimode_",repl, sep ="")
    write(data.frame(modes$locations)[2,],outfile)
    pdffile=paste(repl,"high_cov_relative_increase.pdf",sep = "")
    pdf(file = pdffile)
    plot(modes)
    dev.off()
    print(data.frame(modes$locations)[2,])
    } else {
        temp=subset(high_cov, grepl(repl,high_cov$V1))
        outfile=paste("antimode_",repl, sep ="")
        write(min(temp$V4),outfile)
        print(min(temp$V4))}
}
' > temp0.R

chmod +x temp0.R && Rscript temp0.R

#Filter regions based on antimodes for increase in depth.
replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}')
rm -rf high_cov_intervals2.bed
for repl in $replicons
do
    antimode=$(cat antimode_$repl)
    awk -v antimode=$antimode '{if ($4/$6 > antimode) print $0}' high_cov_intervals_clip.bed >> high_cov_intervals2.bed
done

#Remove regions that match the entire length of its replicon
cp high_cov_intervals2.bed temp
sort -u temp | awk '{if ($5 < $7) print $0}' > high_cov_intervals2.bed

#Get masked coverage (excluding high-covered regions)
bedtools intersect -a coverage.bed -b high_cov_intervals2.bed -v -wa > coverage_nonscaled_masked.txt

#Extract fasta of the regions
bedtools getfasta -fi genomefile.fna -bed high_cov_intervals2.bed > high_cov_intervals.fasta
replicons=$(grep '>' genomefile.fna | sed 's/>//' | awk '{print $1}')

rm -rf replicon_coverage_masked.median.txt replicon_coverage_masked.mean.txt
for repl in $replicons
do
    grep -w "$repl" coverage_nonscaled_masked.txt | sort -nk4,4 | awk -v rep=$repl '{ a[i++]=$4; } END { print rep,a[int(i/2)]; }' >> replicon_coverage_masked.median.txt
    grep -w "$repl" coverage_nonscaled_masked.txt | awk -v rep=$repl  '{x[NR]=$4; s+=$4} {total += $4} END{a=s/NR; for (i in x){ss += (x[i]-a)^2} sd = sqrt(ss/NR); print rep,total/NR,sd}' >> replicon_coverage_masked.mean.txt
done

#Get the non-scaled depths for later copy number calculations
intervals=$(awk '{print $1":"$2"-"$3}' high_cov_intervals2.bed)
rm -rf nonscaled_cov
for i in $intervals
do
    echo $i | sed -e 's/:/\t/g' -e 's/-/\t/g' > temp
    cov_ns=$(bedtools intersect -a temp -b coverage.bed -wa -wb | awk '{x+=$7; next} END{print x/NR}')
    echo $cov_ns >> nonscaled_cov
done

paste high_cov_intervals2.bed nonscaled_cov > temp
mv temp high_cov_intervals2.bed

#Find the expected coverage of intervals, based on their presence of potentially multiple replicons. Some IS elements occur on multiple replicons
makeblastdb -in genomefile.fna -out genomefile -title genomefile -dbtype 'nucl'
#Filtering blastn ID is tricky, since peak calling is based on read mapping. 
blastn -query high_cov_intervals.fasta -db genomefile -outfmt 6 -perc_identity 98 -qcov_hsp_perc 80 -num_threads 50 -out high_cov.blastn
intervals=$(grep '>' high_cov_intervals.fasta | sed -e 's/>//' | sort -u)


rm -f interval_cutoffs
for i in $intervals
do
    echo $i
    grep -w "$i" high_cov.blastn | awk '{print $2}' | sort | uniq -c > temp
    replicons=$(awk '{print $2}' temp | sort -u)
    rm -f temp2
    for a in $replicons
    do
        count=$(grep -w "$a" temp | awk '{print $1}')
        cov=$(grep -w "$a" replicon_coverage_masked.mean.txt | awk '{print $2}')
        stderr=$(grep -w "$a" replicon_coverage_masked.mean.txt | awk '{print $3}')
        sum_rep_cov=$(awk -v c=$count -v cov=$cov -v stderr=$stderr 'BEGIN{print cov+stderr}')
        echo $a $sum_rep_cov >> temp2
    done
    cov_cut=$(awk '{sum+=$2;} END{print sum;$2}' temp2)
    echo $i $cov_cut >> interval_cutoffs
done

#NZ_CP020538.1   2356390 2395720 0.3925214237    39330   0.0360394       2919834
#Filter intervals based on identified cutoffs (sum of background coverage for all replicons they occur on).
intervals=$(awk '{print $1":"$2"-"$3}' high_cov_intervals2.bed)
sed -e 's/\t/:/' -e 's/\t/-/' high_cov_intervals2.bed > temp
rm -f filtered_intervals_pass filtered_intervals_fail
for i in $intervals
do
    echo $i
    cutoff=$(grep -w "$i" interval_cutoffs | awk '{print $2}')
    cov=$(grep -w "$i" temp | awk -F'\t' '{print $2}')
    #cutoff_ns=$(grep -w "$i" interval_cutoffs | awk '{print $3}')
    echo $cutoff
    echo $cov
    if (( $(echo "$cov > $cutoff" |bc -l) ))
    then
        #echo "$i is higher"
        grep -w "$i" temp | awk -v cutoff=$cutoff '{print $0,cutoff}' >> filtered_intervals_pass
    else
        #echo "$i is lower"
        grep -w "$i" temp | awk -v cutoff=$cutoff '{print $0,cutoff}' >> filtered_intervals_fail
    fi
done

sed -i -e 's/:/\t/' -e 's/-/\t/' -e 's/ /\t/g' filtered_intervals_????

#Columns are:
#Replicon start end scaled_cov size main_rep_scaled_cov main_rep_size nonscaled_cov scaled_cutoff nonscaled_cutoff


#Make 100 bp bed file
bedtools intersect -a /DATA_1/tkn/Mobilome_tool/U00096/all_regions_merged_curated.bed -b all_100int.bed -wa -wb | awk -F'\t' 'BEGIN {OFS = FS} {print $1,$2,$3,$4,$5,$6,$7,$1":"$2"-"$3}' > 100int_roi.bed
}



}



#Loop through samples and call the peak_calling function
for treat in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        (cd $treat/$treat$time && peak_calling $treat $time) &
    done
done







#Once all samples have been processed, merge the passing regions of interest (ROI) and inspect the read mappings manually in a genome browser. I used CLC Genomics Workbench.
#Go to the main work dir with the samples in them
cat */*/filtered_intervals_pass > all_regions_merged
#A few ROIs were missed by the peak_calling function. I added these manually.
#The resulting file with curated regions and manual entries is called "all_regions_merged_curated"
#The ROIs were manually inspected and gene annotations were added in the file "all_regions_merged_data_annot_curated.csv"

#Process the manually curated ROIs.
#Upload manual all_regions_merged_curated to thoth and convert to bed
tail -n +2 all_regions_merged_curated | awk -F',' '{print $2,$3}' | sed -e 's/-/ /' -e 's/ /\t/g' > all_regions_merged_curated.bed


#Loop through samples and make 100 bp windowed coverage of ROIs after manually curating ROI coordinates.
for treat in UV Cntr Cet Nal Tet SDS Cu
do
    cd $treat
    for time in 05 20 60
    do
        cd $treat$time
        echo $treat$time
        bedtools intersect -a /DATA_1/tkn/Mobilome_tool/all_regions_merged_curated.bed -b all_100int.bed -wa -wb | awk -F'\t' 'BEGIN {OFS = FS} {print $1,$2,$3,$4,$5,$6,$7,$1":"$2"-"$3}' > 100int_roi.bed
        cd /DATA_1/tkn/Mobilome_tool/clean_scripts/$treat/
    done #time
    cd /DATA_1/tkn/Mobilome_tool/clean_scripts/
done #treat





#Get the clipped coverage of ROIs

clip_cov_fun () {
intervals=$(sed '/chr,range/d' ../../all_regions_merged_curated | awk -F',' '{print $2":"$3}')
treat=$1
time=$2

for int in $intervals 
do
    echo $int | sed -e 's/:/\t/' -e 's/-/\t/' > temp.bed  
    start=$(awk '{print $2}' temp.bed)  
    end=$(awk '{print $3}' temp.bed)  
    grep "U00096.3" clipped_coverage.txt |  grep -w -A100 -B100 "$start" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tstart"}' >> clip_cov_roi
    grep "U00096.3" clipped_coverage.txt |  grep -w -A100 -B100 "$end" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tend"}' >> clip_cov_roi
done #interval
}

for treat in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        rm -f $treat/$treat$time/clip_cov_roi 
        (cd $treat/$treat$time && clip_cov_fun $treat $time) &
    done
done

cat */*/clip_cov_roi > clip_cov_roi_all


#Control regions not expected to have high clipped reads
rm -rf control_regions
#rpoB
echo "U00096.3:4181245-4185273" >> control_regions
#glnA
echo "U00096.3:4056625-4058034" >> control_regions
#ydiJ
echo "U00096.3:1765629-1768685" >> control_regions
#Cryptic prophage Qin
echo "U00096.3:1632277-1652745" >> control_regions
#araB
echo "U00096.3:68348-70048" >> control_regions


clip_cov_cntr_fun () {
cntr_intervals=$(cat ../../control_regions)
treat=$1
time=$2

for int in $cntr_intervals 
do
    echo $int | sed -e 's/:/\t/' -e 's/-/\t/' > temp.bed  
    start=$(awk '{print $2}' temp.bed)  
    end=$(awk '{print $3}' temp.bed)  
    grep "U00096.3" clipped_coverage.txt |  grep -w -A100 -B100 "$start" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tstart"}' >> clip_cov_roi_cntr
    grep "U00096.3" clipped_coverage.txt |  grep -w -A100 -B100 "$end" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tend"}' >> clip_cov_roi_cntr
done #interval
}

for treat in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        rm -f $treat/$treat$time/clip_cov_roi_cntr 
        (cd $treat/$treat$time && clip_cov_cntr_fun $treat $time) &
    done
done


cat */*/clip_cov_roi_cntr > clip_cov_roi_all_cntr



#Get the total mapped bases to normalize clipped reads with #Is not used. There are many factors biasing this number. Especially exonuclease treatment. 
rm -f total_mapped_bases 
for treat in UV Cntr Cet Nal Tet SDS Cu 
do
    for time in 05 20 60
    do
    echo $treat$time 
    bases=$(samtools stats $treat/$treat$time/splitbam/U00096.3.bam | grep 'accurate' | awk '{print $5}') 
    echo $bases 
    echo $treat$time $bases >> total_mapped_bases 
    done #time 
done #treat

#Get the total number of mapped clipped bases. Use this as normalization for clipped reads near ROIs (% of total clipped for chr) 
rm -f total_mapped_clip_bases 
for treat in UV Cntr Cet Nal Tet SDS Cu 
    do
    for time in 05 20 60
    do
       echo $treat$time 
       samtools index $treat/$treat$time/samclip.bam  
       clip_bases=$(samtools view -bh $treat/$treat$time/samclip.bam U00096.3 | samtools stats | grep 'accurate' | awk '{print $5}')
       echo $clip_bases 
       echo $treat$time $clip_bases >> total_mapped_clip_bases 
       done #time
done #treat





#The same but with discordant reads instead of clipped

disc_cov_fun () {
intervals=$(sed '/chr,range/d' ../../all_regions_merged_curated | awk -F',' '{print $2":"$3}')
treat=$1
time=$2

for int in $intervals 
do
    echo $int | sed -e 's/:/\t/' -e 's/-/\t/' > temp.bed  
    start=$(awk '{print $2}' temp.bed)  
    end=$(awk '{print $3}' temp.bed)  
    grep "U00096.3" disc_coverage.txt |  grep -w -A100 -B100 "$start" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tstart"}' >> disc_cov_roi
    grep "U00096.3" disc_coverage.txt |  grep -w -A100 -B100 "$end" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tend"}' >> disc_cov_roi
done #interval
}

for treat in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        rm -f $treat/$treat$time/disc_cov_roi 
        (cd $treat/$treat$time && disc_cov_fun $treat $time) &
    done
done

cat */*/disc_cov_roi > disc_cov_roi_all



#Get the total number of mapped clipped bases. Use this as normalization for clipped reads near ROIs (% of total clipped for chr) 
rm -f total_mapped_disc_bases 
for treat in UV Cntr Cet Nal Tet SDS Cu 
    do
    for time in 05 20 60
    do
       echo $treat$time 
       samtools index $treat/$treat$time/discordant.bam  
       disc_bases=$(samtools view -bh $treat/$treat$time/discordant.bam U00096.3 | samtools stats | grep 'accurate' | awk '{print $5}')
       echo $disc_bases 
       echo $treat$time $disc_bases >> total_mapped_disc_bases 
       done #time
done #treat

#Control regions with discordant reads - not only clipped
disc_cov_cntr_fun () {
cntr_intervals=$(cat ../../control_regions)
treat=$1
time=$2

for int in $cntr_intervals 
do
    echo $int | sed -e 's/:/\t/' -e 's/-/\t/' > temp.bed  
    start=$(awk '{print $2}' temp.bed)  
    end=$(awk '{print $3}' temp.bed)  
    grep "U00096.3" disc_coverage.txt |  grep -w -A100 -B100 "$start" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tstart"}' >> disc_cov_roi_cntr
    grep "U00096.3" disc_coverage.txt |  grep -w -A100 -B100 "$end" | awk -v time=$time -v treat=$treat -v interval=$int '{print $0"\t"treat"\t"time"\t"treat""time"\t"interval"\tend"}' >> disc_cov_roi_cntr
done #interval
}

for treat in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        rm -f $treat/$treat$time/disc_cov_roi_cntr
        (cd $treat/$treat$time && disc_cov_cntr_fun $treat $time) &
    done
done

cat */*/disc_cov_roi_cntr > disc_cov_roi_all_cntr








#From the merged bedfile from R, extract all regions from all samples for a very large dataset.

#Manually defined control regions (no increased mobilome coverage in any samples) were defined.
#cntr_intervals=$(cat control_regions)

merge_data_fun () {
treat=$1
time=$2

intervals=$(sed '/,chr,range/d' ../../all_regions_merged_curated | awk -F',' '{print $2":"$3}')

for int in $intervals
do

echo $int | sed -e 's/:/\t/' -e 's/-/\t/' > temp.bed
#Get fasta of manual entries for blasting
bedtools getfasta -fi genomefile.fna -bed temp.bed > temp.fasta
blastn -query temp.fasta -db genomefile -outfmt 6 -perc_identity 98 -qcov_hsp_perc 80 -num_threads 100 -out temp.blastn
grep -w "$int" temp.blastn | awk '{print $2}' | sort | uniq -c > temp
replicons=$(awk '{print $2}' temp)
rm -f temp2 all_regions_cutoffs
#rm -f temp2 all_regions_cutoffs_cntr
for a in $replicons
do
    count=$(grep -w "$a" temp | awk '{print $1}')
    cov=$(grep -w "$a" replicon_coverage_masked.mean.txt | awk '{print $2}')
    stderr=$(grep -w "$a" replicon_coverage_masked.mean.txt | awk '{print $3}')
    sum_rep_cov=$(awk -v c=$count -v cov=$cov -v stderr=$stderr 'BEGIN{print cov+stderr}')
    echo $a $sum_rep_cov >> temp2
done #replicons
cov_cut=$(awk '{sum+=$2;} END{print sum;$2}' temp2)
echo $int $cov_cut >> all_regions_cutoffs
#Get the scaled depths
echo $int | sed -e 's/:/\t/g' -e 's/-/\t/g' > temp
size=$(echo $int | sed -e 's/:/\t/g' -e 's/-/\t/g' | awk '{print $3-$2}')
replicon=$(echo $int | awk -F':' '{print $1}')
main_rep_size=$(grep "$replicon" genome.bed | awk '{print $2}')
nonscaled_cov=$(bedtools intersect -a temp -b coverage.bed -wa -wb | awk '{x+=$7; next} END{print x/NR}')
nonscaled_cutoff=$(grep $int all_regions_cutoffs | awk '{print $2}')
echo $int $treat $time $treat$time $size $main_rep_size $nonscaled_cov $nonscaled_cutoff >> all_regions_data
done
}


for treat in UV Cntr Cet Nal Tet SDS Cu
do
    for time in 05 20 60
    do
        rm -f ./*/*/all_regions_data
        (cd $treat/$treat$time && merge_data_fun $treat $time) &
    done
done


cat */*/all_regions_data > all_regions_merged_data
#cat */*/all_regions_data_cntr > all_regions_merged_data_cntr


#pOLA52 coverage for calculating relative copy numbers
rm -f pOLA52_cov
for treat in UV Cntr Cet Nal Tet SDS Cu
do
    echo $treat
    for time in 05 20 60
    do
        echo $time
        cov=$(grep NC_010378.1 $treat/$treat$time/replicon_coverage.txt | awk '{print $2}')
        echo $treat$time $cov >> pOLA52_cov
    done
done

rm -f chr_cov
for treat in UV Cntr Cet Nal Tet SDS Cu
do
    echo $treat
    for time in 05 20 60
    do
        echo $time
        cov=$(grep U00096 $treat/$treat$time/replicon_coverage.txt | awk '{print $2}')
        echo $treat$time $cov >> chr_cov
    done
done



#Combined 100-windowed coverage over regions of interest
rm -f all_treat_100int_roi.bed
for treat in UV Cntr Cet Nal Tet SDS Cu
do
    echo $treat
    for time in 05 20 60
    do
        echo $time
        cat $treat/$treat$time/100int_roi.bed | awk -F'\t' -v var=$treat$time 'BEGIN {OFS = FS} {print var,$0}' | sort -u >> all_treat_100int_roi.bed
    done
done



