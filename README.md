# KARGA
Multi-platform Toolkit for k-mer-based Antibiotic Resistance Gene (ARG) Analysis of High-throughput Sequencing Data

# Installation
KARGA requires the Java Virtual Machine (https://www.java.com/en/). The .class files available on this GitHub have been compiled on MS Windows 10 using 64-bit javac v.15.

# Usage
KARGA can be launched from the command line. The minimum input is a read file in (optionally gzipped) FASTQ format, which is automatically detected if the extension is .fastq or .gz. Without other parameters, the MEGARes v.2.0 (https://megares.meglab.org/) is used with a default value of k=17.
- Individual read classification can be done typing "java KARGA_ReadMapper readfile.fastq"
- Resistome mapping can be done typing "java KARGA_ResistomeMapper readfile.fastq"
Both java classes accept the following optional parameters: "k:your_k_value"; "d:your_arg_db_fasta"; "f:your_read_fastq". Note that if another ARG database is chosen, the corresponding FASTA file needs to contain headers consistent with MEGARes ARG ontology.

# Output
- KARGA_RadMapper prints a CSV file --one line per read-- with the following fields: read_idx, arg_type, arg_class, arg_mechanism, arg_group, type_kmer_hits/mapped/total, class_kmer_hits/mapped/total, mechanism_kmer_hits/mapped/total, group_kmer_hits/mapped/total.
- KARGA_ResistomeMapper prints a CSV --one line per ARG-- with the following fields: arg_idx, percent_gene_coverage, median_kmer_depth. Note that ARGs with coverage below 1% are not printed.

