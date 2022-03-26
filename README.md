# Platanus_read_correction README.md

## Description
Read error correction tool derived from Platanus assembler.

## Installation
```sh
make
cp platanus <installation_path>
```

## Synopsis
```sh
platanus correct -f reads.fastq 2>log
```

## Options
    -o STR               : suffix of output file (def .corrected, length <= 200)
    -f FILE1 [FILE2 ...] : reads file (fasta or fastq, number <= 100)
    -k INT [INT ...]     : k values (<= 32, def 18)
    -e FLOAT             : max_edit_distance / read_length (<= 1, def 0.03)
    -q INT               : quality value cutof (def 0)
    -t INT               : number of threads (<= 100, def 1)
    -m INT               : memory limit(GB, >= 1, def 4)
    -h, -help, --help    : print usage

## Outputs
INPUT_FILE1.corrected [INPUT_FILE2.corrected ...]

## Author
Rei Kajitani at the Tokyo Institute of Technology wrote key source codes.
