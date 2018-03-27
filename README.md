# protein-express
Compare proteomic data to metagenomic bins to understand patterns of protein expression across taxonomic groups

## Dependencies
Include these in the path
- blastp
- makeblastdb
- prodigal

## Arguments
- -p, --prot_dirs <br />
List of directories for different proteomic datasets <br />
Directory name must be the proteome name found in filenames, e.g., proteome1 <br />
Necessary files in each directory: <br />
  - Postnovo output table, reported_df.tsv
  - Database search results (e.g., proteome1.metagenome1.reads.tsv, proteome1.metagenome1.DBGraphPep2Pro.tsv)
  - Corresponding fastas of the translated genes that produced PSMs (e.g., proteome1.metagenome1.reads.fasta, proteome1.metagenome1.DBGraphPep2Pro.tsv)
- -b, --bin_dir <br />
Directory solely containing metagenomic bin fastas <br />
- -o, --out <br />
Output directory