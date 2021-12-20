# KARGA
Multi-platform Toolkit for K-mer-based Antibiotic Resistance Gene (ARG) Analysis of High-Throughput Sequencing Data (v.1.01 Dec 20 2021)

# Installation
KARGA requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
- KARGA can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. 
- Without other parameters, the MEGARes database (https://megares.meglab.org/) is used with a default value of k=17. Please download the latest MEGARes release here: (curl) https://megares.meglab.org/download/index.php; https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta.
- If you wish to classify mobile genetic elements (MGEs), we recommend to enable the multinomial classification option "m:y" which outputs multiple weighted hits. Several MGE reference databases can be used, e.g. ICEBerg (https://db-mml.sjtu.edu.cn/ICEberg/; curl https://db-mml.sjtu.edu.cn/ICEberg2/download/ICE_seq_all.fas).
- By default, the program outputs individual read classification as well as mapping of the resistome database given in input.
- If you want to classify antibiotic resistance in genes due to specific point mutations (e.g. in housekeeping genes), we recommend to use the dedicated KARGVA module (https://github.com/DataIntellSystLab/KARGVA).
- Type "java KARGA readfile.fastq" for the default execution.
The java class accepts the following optional parameters: "k:your_k_value" (positive integer for k-mer length); "d:your_db_fasta" (any ARG/MGE database in FASTA format where resistance annotation is specified in the header); "f:your_read_fastq" (read file in FASTQ format with any file extension); "r:[y,yes,n,no]" (if you want to print or omit individual read classification, as the program is slightly faster when this print is omitted); "m:[y,yes,n,no]" (if you want to print the top-scoring ARG hits for each read, not only the best one); "i:your_value" (number of iterations to calculate frequency threshold from customized random string hit distribution, default is 125,000); "s:seed" (integer seed value for random number generator to ensure code reproducibility). For large databases, we recommend to use -Xmx16GB or larger.

# Output
- inputFileName_KARGA_mappedReads.csv : a CSV file --one line per read-- with the following fields: Read_Idx, GeneProbability/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal, GeneAnnotation.
- inputFileName_KARGA_mappedGenes.csv : a CSV --one line per ARG-- with the following fields: GeneIdx, PercentGeneCovered, AverageKMerDepth. Note that ARGs with coverage below 1% are not printed; recommended ARG coverage is 80%.

# Citation
M. Prosperi and S. Marini, "KARGA: Multi-platform Toolkit for k-mer-based Antibiotic Resistance Gene Analysis of High-throughput Sequencing Data," 2021 IEEE EMBS International Conference on Biomedical and Health Informatics (BHI), 2021, pp. 1-4, doi: 10.1109/BHI50953.2021.9508479.
Available at: https://ieeexplore.ieee.org/document/9508479
