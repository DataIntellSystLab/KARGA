# KARGA
Multi-platform Toolkit for k-mer-based Antibiotic Resistance Gene (ARG) Analysis of High-throughput Sequencing Data

# Installation
KARGA requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
KARGA can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. Without other parameters, the MEGARes v.2.0 (https://megares.meglab.org/) is used with a default value of k=17. By default, the program outputs individual read classification as well as mapping of the resistome database given in input.
- Type "java KARGA readfile.fastq" for the default execution.
The java class accepts the following optional parameters: "k:your_k_value" (k-mer length(; "d:your_arg_db_fasta" (any ARG database in FASTA format where resistance annotation is specified in the header); "f:your_read_fastq" (read file in FASTQ format with any file extension); "r:[y,yes,n,no]" (if you want to print or omit individual read classification, as the program is slightly faster when this print is omitted).

# Output
- inputFileName_KARGA_mappedResistome.csv : a CSV file --one line per read-- with the following fields: read_idx, ARG_hits/mapped/total, ARG_annotation.
- inputFileName_KARGA_mappedReads.csv : a CSV --one line per ARG-- with the following fields: ARG_idx, percent_ARG_coverage, median_kmer_depth. Note that ARGs with coverage below 1% are not printed.

