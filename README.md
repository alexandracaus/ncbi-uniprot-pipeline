# Connecting NCBI to UniProt: A Python Pipeline for Protein Mapping

An automated cascading computational pipeline designed to maximize protein identifier mapping between NCBI (GenBank/RefSeq) and UniProtKB repositories. This framework is specifically engineered to rescue functional metadata for custom clinical bacterial isolates facing repository downsizing policies.

## Features
- **Cascading Architecture:** Sequentially interrogates direct cross-references (Case 0), NCBI Identical Protein Groups (Case 1), and local ClusteredNR databases (Case 2).
- **Data Rescue:** Bypasses broken online indexes using programmatic UniParc and UniRef90 translation bridges.
- **FAIR Compliance:** Generates structured logging to guarantee absolute data provenance and reproducibility.

## Requirements & Installation
The pipeline requires Python 3.x and the `requests` library to handle remote REST API communications.

```bash
pip install requests
```

## Usage
The script supports execution via input files containing lists of protein identifiers or by passing individual protein accessions directly through the command line.

### Option 1: Processing an Input File (-f)
To process a file containing unmappable structural isolates or historical failure lists (e.g., in .lst or .txt format):

```bash
python ncbi_uniprot_pipeline.py -f file.txt
```

### Option 2: Processing Individual Protein Codes (-p)
To test or execute a quick mapping lookup for specific protein accessions directly:

```bash
python ncbi_uniprot_pipeline.py -p code
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.
