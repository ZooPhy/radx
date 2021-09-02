bam_files = [x for x in os.listdir(PATH_TO_BAM) if x[-17:]=="masked.sorted.bam"]

with open (r"metrics_masked.tsv", "w", newline = "") as file:
    tsv_file = csv.writer(file, delimiter ="\t")
    tsv_file.writerow(["File", "Breadth of coverage", "Total read count", "Mean reads"]
    for file in bam_files:
        breadth = run_cmd(["samtools", "depth", "-a", file, "|", "awk" , "{c++; if($3>10) total+=1}END{print (total/(c))*100}"])
        count = run_cmd(["samtools", "view", "-c", "-F", "4", file])
        mean = run_cmd(["samtools", "depth", "-a", file, "|", "awk", '{c++;s+=$3}END{print s/c}')
        tsv_file.writerow([File, breadth, count, mean])
      
"""
# bash version
echo -e "File \tBreadth_of_Coverage \tRead_count \tMean_reads" >> metrics_masked_new.txt
for file in ~/*.masked.sorted.bam;
do
	breadth=$(samtools depth -a "${file}" | awk '{c++; if($3>10) total+=1}END{print (total/(c))*100}')
	count=$(samtools view -c -F 4 "${file}")
	mean=$(samtools depth -a "${file}" | awk '{c++;s+=$3}END{print s/c}')
	echo -e "${file} \t${mean} \t${breadth} \t${count}" >> metrics_masked_new.txt
done
"""
