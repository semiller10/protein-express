# protein-express
Compare proteomic data to metagenomic bins to understand differential protein expression by taxonomic groups

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
  - Postnovo output table
  - Database search results (e.g., proteome1.metagenome1.tsv)
  - Corresponding fastas of the translated genes that produced PSMs (e.g., proteome1.metagenome1.fasta)
- -b, --bin_dir <br />
Directory solely containing metagenomic bin fastas <br />
- -o, --out <br />
Output directory
