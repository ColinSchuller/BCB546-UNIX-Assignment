#UNIX Assignment

##Data Inspection

###Attributes of `fang_et_al_genotypes`

    wc fang_et_al_genotypes.txt 
    file fang_et_al_genotypes.txt 
    du -h fang_et_al_genotypes.txt
    tail -n +1 fang_et_al_genotypes.txt | awk -F "\t" '{print NF; exit}'

    tail -n +6 Mus_musculus.GRCm38.75_chr1.gtf | awk -F "\t" '{print NF; exit}'

By inspecting this file I learned that:

1.  There are 2783 lines, 2744038 words, and 11051939 bytes in the file.
2.  The file is in ASCII text, with very long lines.
3.  The file is 6.7M in human file terms.
4.  And the file has 986 columns.

###Attributes of `snp_position.txt`

    wc snp_position.txt 
    file snp_position.txt 
    du -h snp_position.txt
    awk -F "\t" '{print NF; exit}' snp_position.txt

    tail -n +6 Mus_musculus.GRCm38.75_chr1.gtf | awk -F "\t" '{print NF; exit}'

By inspecting this file I learned that:

1.  There are 984 lines, 13198 words, and 82763 bytes in the file.
2.  The file is ASCII text
3.  The file is 49k in human file terms.
4.  The file has 15 columns.

##Data Processing

###Maize Data

    #Grep the important things into a new file
        grep -E "(Sample_ID|ZMMIL|ZMMLR|ZMR)" fang_et_al_genotypes.txt > maize_genotypes.txt 

    #Transposed the data first so that the columns become rows
        awk -f transpose.awk maize_genotypes.txt > transposed_maize_genotypes.txt

    #sort the new transposed file
        sort -k1,1 transposed_maize_genotypes.txt > transposed_maize_genotypes_sorted.txt

    #Cut the columns out of the snp file before 
        cut -f 1,3,4 snp_position.txt > snp_position_cut.txt

    #Remove the column title from the first row of the snp file
         tail -n +2 snp_position_cut.txt > snp_position_cut_removedheader

    #sort the snp file
         sort -k1,1 snp_position_cut_removedheader.txt > snp_position_cut_removedheader_sorted.txt

    #Merged the two files together
        join -1 1 -2 1 snp_position_cut_removedheader_sorted.txt transposed_maize_genotypes_sorted.txt > genotypes_snp_merged.txt

    #sorted the merged file again by position 
         sort -k3,3n genotypes_snp_merged.txt > genotypes_snp_merged_sortedposition.txt

    #Performed a for loop to seperate the increasing position order file into the 10 chromosomes
    for i in {1..10}; do awk -v chr="$i" '$2 == chr { print $0 > "maize_chr"chr".txt"}' genotypes_snp_merged_sortedposition.txt; done

    #sorted the merged file again by descending position 
         sort -r -k3,3n genotypes_snp_merged.txt > genotypes_snp_merged_sortedposition2.txt

    #replaced the ? with - using sed command
        sed 's/?/-/g' genotypes_snp_merged_sortedposition2.txt > genotypes_snp_merged_sortedposition2_replaced.txt

    #Performed a for loop to seperate the increasing position order file into the 10 chromosomes
    for i in {1..10}; do awk -v chr="$i" '$2 == chr { print $0 > "maize_2chr"chr".txt"}' genotypes_snp_merged_sortedposition2_replaced.txt; done

    #1 file with all snps with unknown positions
        awk '$3 ~ /^unknown$/' genotypes_snp_merged_sortedposition.txt > genotypes_snp_merged_unknownposition.txt

    #1 file with all the snps with multiple locations
        awk '$3 ~ /^multiple$/' genotypes_snp_merged_sortedposition.txt > genotypes_snp_merged_multipleposition.txt

Isolates maize samples specified by the user, then transposes the data so the columns become rows. Following that it sorts the data by sample ID. Then to join it to the snp file we first have to remove the header row of the snp file so it does not interfere with the join command. Sort the snp file by sample ID, then we join the two files. From there it sorts in ascending order by snp position and performs a for loop to separate out the data by chromosome into 10 different files (number can be specified). After that it repeats the same process but this time it replaces all ? in the file with -. Then it sorts the merged file with - in descending order. Then performs the for loop like before. Finally the code takes the merged file and separates out the snps in unknown positions into a new file, and snps with multiple positions into a new file.

###Teosinte Data

    #Grep the important things into a new file
        grep -E "(Sample_ID|ZMPBA|ZMPIL|ZMPJA)" fang_et_al_genotypes.txt > teosinte_genotypes.txt 

    #Transposed the data first so that the columns become rows
        awk -f transpose.awk teosinte_genotypes.txt > transposed_teosinte_genotypes.txt

    #sort the new transposed file
        sort -k1,1 transposed_teosinte_genotypes.txt > transposed_teosinte_genotypes_sorted.txt

    #Merged the two files together
        join -1 1 -2 1 snp_position_cut_removedheader_sorted.txt transposed_teosinte_genotypes_sorted.txt > genotypes_snp_merged_teosinte.txt

    #sorted the merged file again by position 
         sort -k3,3n genotypes_snp_merged_teosinte.txt > genotypes_snp_merged_teosinte_sortedposition.txt

    #Performed a for loop to seperate the increasing position order file into the 10 chromosomes
    for i in {1..10}; do awk -v chr="$i" '$2 == chr { print $0 > "teosinte_chr"chr".txt"}' genotypes_snp_merged_teosinte_sortedposition.txt; done

    #sorted the merged file again by descending position 
         sort -r -k3,3n genotypes_snp_merged_teosinte.txt > genotypes_snp_merged_teosinte_sortedposition2.txt

    #replaced the ? with - using sed command
        sed 's/?/-/g' genotypes_snp_merged_teosinte_sortedposition2.txt > genotypes_snp_merged_teosinte_sortedposition2_replaced.txt

    #Performed a for loop to seperate the increasing position order file into the 10 chromosomes
    for i in {1..10}; do awk -v chr="$i" '$2 == chr { print $0 > "teosinte_2chr"chr".txt"}' genotypes_snp_merged_teosinte_sortedposition2_replaced.txt; done

    #1 file with all snps with unknown positions
        awk '$3 ~ /^unknown$/' genotypes_snp_merged_teosinte_sortedposition.txt > genotypes_snp_merged_teosinte_unknownposition.txt

    #1 file with all the snps with multiple locations
        awk '$3 ~ /^multiple$/' genotypes_snp_merged_teosinte_sortedposition.txt > genotypes_snp_merged_teosinte_multipleposition.txt

Isolates maize samples specified by the user, then transposes the data so the columns become rows. Following that it sorts the data by sample ID. Then we join the two files because the snp file has already been cleaned up and sorted from the previous step. From there it sorts in ascending order by snp position and performs a for loop to separate out the data by chromosome into 10 different files (number can be specified). After that it repeats the same process but this time it replaces all ? in the file with -. Then it sorts the merged file with - in descending order. Then performs the for loop like before. Finally the code takes the merged file and separates out the snps in unknown positions into a new file, and snps with multiple positions into a new file.

