# NortheasternCoOp
Using bash script to analyze illumina NGS data to leverage the 16S ribosomal RNA gene of prokaryote present in the samples analyzed. The 16S rRNA gene was focused on to help used to understand the phylogenic structure of species content and diversity because it offers a proxy/heuristic to be able to conduct analysis, as opposed to needing to sequence the entire genome of each and every organism present in each of the hundreds of samples that require analysis.



Developing a 16S Metagenomics Pipeline for Microbial Analyses and Testing with Diverse Coral Data


Patrick Choi


Introduction 
Metagenomic analysis is broadly defined as the study of genomic DNA material recovered directly form environmental samples. In this particular application of metagenomics, we leveraged the 16S ribosomal RNA gene of prokaryote present in the samples analyzed. The 16S rRNA gene was focused on to help used to understand the phylogenic structure of species content and diversity because it offers a proxy/heuristic to be able to conduct analysis, as opposed to needing to sequence the entire genome of each and every organism present in each of the hundreds of samples that require analysis. The 16S region is approximately 1,500 base-pairs (bp) long and contains nine variable regions that are located between more easily identifiable conserved regions of DNA. These variable regions are used in phylogenetic classification such as genus or species in microbial populations. 
	The Illumina HiSEQ platform produces shorter reads length than other technologies, such as PacBIO however, this type of Next Generation Sequencing (NGS) technology produces data which is more accurate and can be analyzed with software tools such as QIIME (Quantitative Insights into Microbial Ecology), USEARCH, and CLUSTAL, among others. QIIME is an open-source bioinformatics pipeline used to perform microbe analysis from sequencing data. The open sourced platform, while rudimentary, can be used to analyze and interpret the structural dynamics of these communities, and provide insights as to how the community diversity and population dynamics affect larger systems, such as the interactions between hosts (corals) and symbionts (dinoflagellate algae). The process of creating a pipeline requires that the results of one step to be used as the input for the next step. This type of automation requires that strict information control across of temporary files produced be maintained and that the python script developed adhere to specific parameters that remain universal. Resulting outputs can subsequently be analyzed using statistical methods via python or R. 

Background 
	From the 16S ribosomal RNA gene, we were given data regarding the Variable 3 (V3) and Variable 4 (V4) region of coral bacteria as well as potential 18S contaminants from Eukaryotes present in the extractions. The targeted region was sequenced using the Illumina HiSeq2500 platform, resulting in 300 million Paired-End (PE) reads, that were 250bp long, that were produced from DNA fragments of less than 500bp. These particular parameters should have results in significant overlap between the forward facing and reverse sequenced reads, that then could be overlapped to create a single read. We used QIIME, which was pre-installed on Defiance for the analysis, to conduct a majority of the work, in addition to, custom scripts that were developed. The primary focus was to modify the data, through the development of a pipeline, so that it would to adhere to the requirements for QIIME for further processing. We initially started by demultiplexing the PE data, overlapping the reads to produce a single longer read, and then creating a mapping file that would be used to analyze the data via QIIME. Some of the last and experimental steps included performing statistical analyses using R. 
	The sequencing data were in the fastq format whereas the metadata was in a xlsx format. The R1.fastq and the R2.fastq files were provided to us on the server. There were two lanes of sequencing data with overlapping barcodes that required addressing. Because of the overlap in barcodes, as well as potentials for different trimming lengths, we needed to keep the two lanes of sequencing data separated. Additionally, there were two xlsx files, one for lane 1, which contained 4 sheets, and one for lane 2, which contained 5 sheets (there was an additional sheet with forward and reverse primer sequence, which was used later to trim the sequence out as well as identify appropriately facing reads). We needed to modify these files into a new file in order to use QIIME, but converting it from a xlsx file to a Unix based text file via custom programing. QIIME provides scripts to convert sequencing data into usable data. Because of the nature of the mistakes made in the experimental design as well as mismatching heads produced by the sequencing facility, we had to customize the scripts to adapt them to our specific needs. Below are some of the scripts used.
	
Validate_mapping_file.py
This script validates our mapping files to ensure that all the required fields are met and is formatted correctly. The results were given to us in html format, which allowed for easier viewing in order to understand the indicated error or warning. This particular feature was not utilized, as all analysis was conducted in the Command Line Interface (CLI)

Extract_barcode.py
This script extracted the barcodes, as in the actual data, from the main sequencing fastq format files, and stored it in a separate file. The barcodes located in separate files is a  necessary feature to be able to use with the split_libraries_fastq.py, which was the script we used in the next step. 

Split_libraries_fastq.py
 This script demultiplexed the fastq sequence. This script requires the barcode file from the extract_barcode.py script. If the mapping file is not validated then this split_libraries script will not run. It is important to note that this script produces only one file, and redesigned the header of each sequence in order to retain sample information. This is important later during read mapping, since it allows up to determine where each read came from. 

Count_seqs.py
Once we have the seqs.fna file generated from the split_libraries step, we need to count the sequences that are in the file, using SampleID as a counter. This is important so that we can identify how much sequencing data we have that is usable versus unassigned. Count_seqs.py does this and provides the output on the terminal or we can write the results into a separate field. 
Operation Taxonomical Units (OTUs) are the basis of identifying species within our communities/samples. An OTU is broadly defined as a species or group of species that is often used when only DNA sequence data is available. There are three protocols in QIIME to pick OTUs, and we depended on de novo OTU picking since the DNA regions we focused on we not of a sufficient length to map against existing reference with high confidence. 

1.)	pick_de_novo_otus.py- this can be used for de novo OTU picking. 

2.)	Pick_closed_reference_otus.py- to be used for closed reference OTU picking

3.)	Pick_open_reference_otus.py- to be used for open reference OUT picking.

Core_diversity_analysis.py – This script is used for diversity analysis. It requires a BIOM file, mapping file, and an output path.
Additionally, we used the FastQC, a tool that allows quality control check on raw sequencing data before starting the analysis on QIIME. This provides an overview of potential problems.  Tools like Trimmomatic, help to trim the reads for Illumina NGS data. We can remove the leading or the tail end or the N bases as well as control for quality of the reads. 


Development 
In processing the data and getting it ready for QIIME we leveraged several command line tools that are present in most UNIX environments, such as sed, awk, cut, paste, among many others. The particular environment in which we developed tools was BASH, Bourne Again SHell, which comes preinstalled with CentOS Linux. As an example, awk has its own programming language built into the editor program, which is a lot like other text file processing utilities such as vim, and allow for compilation into pipeline with other UNIX commands via the “|”. This facilitated the development of BASH scripts to process data streams more efficiently then opening several large files in memory using Perl, Python, or R.
Initial Steps:

1.)	The xlsx files (from both lane files) were separately converted into UNIX compatable a csv files and white space was converted into tabs [\t] so as to not confuse it with multiple spaces [\s] during subsequent processing.

2.)	For the purposes of QIIME analysis, we then converted the csv files to a tsv file using the command: “cat filename.csv | sed’s/,/\t/g’>filename.tsv” 

3.)	Lastly, we converted the .tsv file to a UNIX text file format using dos2unix function: dos2unix –n a.tsv e.txt

Generating the Barcode file

From the text files we extracted the “Forward Barcode”, the “Reverse Barcode”, “SampleID”, which were tab separated. In order to extract the required column, we used the following line of code:
Awk’{print $6”\t”$5”\t”$2}’ Lane1_IndexFile__A_T0.txt>Lane1_BarcodeFile1.txt
This extracted the column numbers and printed them into a new temporary file in the order we provided, using the streaming function of joining command in UNIX via [>]. 

Next, we extracted all four of the sheet’s txt files into a new txt file titled Lane1_BarcodFile. We then combined all of the 4 barcode files into one file with:
Awk ‘FNR==1 &&NR!{next;}{print}’Lane_BarcodeFile*.txt > Lane1_BarcodeFile.txt
The resulting single Barcode file will have all Lane1 sheets data with the required output. We repeated the same procedure for the Lane 2 files. 

Generating the Mapping Files:
The mapping file headers need to specifically contain these header names. If the header names are not correctly specified, then the program will return an error: 
A.	#sampleID
B.	BarcodeSequence
C.	LinkerPrimerSequence
D.	Metadata columns- Optional
E.	Description

It is the same procedure to generate Barcode files but there needs to be different columns and different combinations
awk ‘{print $2”\t”$6$5”\t””””\t””+A_T0_W”$1”_”$2”_””F”$4”R”$$3}’Lan1_IndexFile_+A_T0.txt>Lane1_MappingFile1.txt
In order to remove the headers of the file
- echo “$(tail –n+2Lane1_MappingFile1.txt):>Lane1_MappingFile1.txt

The next command is used to create a new header:
sed #SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription’Lane1_MappingFile1.txt
We did this for all four of the lane mapping files, Lane1_MappingFile1 – 4. After that we combined all of the mapping files using the following lines of code:
 awk’FNR==1 && NR!{next;}{print:}’Lane1_MappingFile*.txt>Lane1_MappingFile.txt
The result was a single mapping file with all of the sheet data and desired output. This format can be used for the validation process as well as subsequent demultiplexing. We repeat this procedure for all of the Lane2 files in order to generate a single Lane2 mapping file, that has off of the pertinent data.

According QIIME, the headers need to be as the following: SampleID, Barcode Sequence, Linker Primer, Sequence, and Description. The data is tab separated. The SampleID field has to be alphanumeric with a “.”. The other fields will accept alphanumeric, periods”.”, underscore “_”, percent “%”, plus “+”, space “ “, semicolon “;”, colon “:”, comma “,”, or forward slash “/”. Because of the format, we need to convert the “_” in the SampleID to “.”. We did so with:

Awk –F’\t’-vOFS=’\t’’{gsub(“_”,”.”,$1);print}’Lane1_MappingFile.txt>Lane1_MappingFile_new.txt
Next, we needed to validate the generated Mapping file. To do this we can use the validate_mapping_file.py script. This will check to see if there are any errors. If there is, then the script will generate an error log. Here is the code for validating the mapping file:
Validate_mapping_file.py –m Lane1_MappingFile_new.txt –p –o validate_mapping_output/

We also need to generate the barcode file in a fastq format in order to use the rest of the scripts needed to process the data. To generate the barcode in fastq format we used:
extract_barcodes.py –f R1.fastq –r R2.fastq –c barcode_paired_end –bc1_len 12 –bc2_len 12 –o processed_seqs

We used both the R1 and the R2 file with the paired end barcode. We need to provide a 12 bp length for both files. The output of the barcode.fastq file will generate a sequence of 24 bp. This is due to a combination of forward and reverse primer.
The next step is to demultiplex the data. We can demultiplex the data with the split_libraries_fastq.py script. In order for this script to run we need to provide the 24 barcode file generated from the previous sequence. 
split_libraries_fastq.py –I R1.fastq –m Lane1_MappingFile_new.txt –b processed_seqs/barcodes.fastq –barcode_type24 –o slout_lane1
This also was required for the R2.fastq file as well, so that reads in both directions were processed. The split_libraries script generated three different files: 
1.)	Histogram.txt
2.)	Seqs.fna
3.)	Split_library_log.txt

The code that was used to count the sequence number from the seqs.fna file was: Count_seqs.py –I seqs.fna

Conduct OTU picking
OTU picking uses an open-reference OTU picking protocol that searches the reads against a reference database. We can take the demultiplexed sequences and cluster these sequences into OTUs. There are three ways to do this in QIIME. We can use de novo, closed reference, or open reference OTU picking. Even though open reference OTU picking is the more popular method, we used de novo picking. We did this for several reasons, most importantly, we did not have a closely related reference map. Additionally, a lot of our reads were shorter than the references available and did not come from samples collected and processed in a similar enough manner. 

Diversity analyses.
We can use the core_diversity_analyses.py script to generate index.html file which provides us with different analysis. We might get a runtime warning, but this isn’t something to worry about because it is related to mismatches in field that are not relevant to out analysis, such as sequence descriptions, because we omitted this data from the mapping file that was made. It is important to note that the results are a function of analysis that treats samples as independent of eachother. We can categorize samples as a function of their metadata by passing the headers from our mapping file with a –c parameter. 

Pipeline creation
We can use R for downstream analysis, supervised learning, and distance matrix comparison, although this was not something that we did during this particular set of analysis. These analysis’ can be performed either with a single script or using multiple scripts that interact with each other, which allows the user to provide the parameters and the use of the results form previous command. The pipeline processes the input sequence and generates the output. This entire pipeline also can accommodate a master parameters file, which is specified in a centrally accessible location, and can provide specific details and command parameter for each step of the pipeline in QIIME. This type of automation can prove to be difficult to implement however, since each command has to be thoroughly reviewed to ensure that inputs match outputs, and that file formats specified upstream do not adversely affect inputs generated in downstream analysis. 

CoOp experience
It was interesting to work on NGS data and to learn how to process the data. This Co-op has given me that exposure. More importantly I was able to understand the scripts used by QIIME and utilize R for data analysis, which is a very useful skillset. I was constantly anxious that my education might not be sufficient for a job in bioinformatics, but this Co-op has given me some courage to pursue my career. 
















Bibliography

Caporaso, Gregory, Christian L. Lauber, and Williams A. Walters. "Ultra-high-throughput microbial community analysis on the Illumina HiSeq and MiSeq platforms ." International Society for Microbial Ecology, 2012: 1621-1624.

Fadrosh, Douglas W., Bing Ma, and Naomi Sengamalay. "An improved dual-indexing approach for multiplexed 16s rRNA gene sequencing on the Illumina MiSeq platform." Microbiome , 2014: 2-6.

Gignoux-Wolfsohn, Sarah A., and Steven V. Vollmer. "Identification of Candidate Coral Pathogens on White Band Disease-Infected Staghorn Coral." PLoS ONE, 2015.

Illumina. 16sMetagenomic Sequencing Library Preparation. Manual, Illumina, San Diego: Illumina.

QIIME. QIIME Quantitative Insights into Microbial Ecology. August 22, 2016. QIIME.org (accessed August 22, 2016).






