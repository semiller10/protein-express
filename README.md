# protein-express
Compare proteomic data to metagenomic bins to understand patterns of protein expression across taxonomic groups and environments

## Dependencies
Include these in the path
- blastp
- makeblastdb
- prodigal

## Arguments
### Required
- -p, --prot_dirs <br />
List of directories for different proteomic datasets <br />
Directory name must be the proteome name found in filenames, e.g., proteome1 <br />
The directory names can be placed in a file and inserted into the command line: <br />
`$(cat prot_dirs.txt | tr '\n' ' ')` <br />
Necessary files in each directory: <br />
  - Postnovo output table, reported_df.tsv
  - Database search results (e.g., proteome1.metagenome1.reads.tsv, proteome1.metagenome1.DBGraphPep2Pro.tsv)
  - Corresponding fastas of the translated genes that produced PSMs (e.g., proteome1.metagenome1.reads.fasta, proteome1.metagenome1.DBGraphPep2Pro.tsv) <br />
- -b, --bin_dir <br />
Directory solely containing metagenomic bin fastas <br />
- -o, --out <br />
Output directory
### Optional
- -s, --state <br />
Table (tsv with header) relating each proteome name to a labeled state , e.g., <br />

| Proteome  | Environment |
|:---------:|:-----------:|
| proteome1 | org_soil    |
| proteome2 | min_soil    |

- -t, --threads <br />
Number of threads, default = 1
