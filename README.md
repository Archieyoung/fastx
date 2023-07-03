# Fastx

Fastx is an ultrafast toolkit for manipulating FASTA/Q file.

## Features

- Ultrafast
- Seamlessly parsing both FASTA and FASTQ formats
- Support both single-end and paired-end FASTA/Q files
- Parsing by number of bases or number of reads(head, sample, split)


## Build

```sh
git clone --recursive git@github.com:Archieyoung/fastx.git
cd fastx
mkdir build && cd build
cmake .. && make -j
```

You can find the executable `fastx` in build/bin.

## Usage

```
./fastx 
fastx v0.3.2

  fasta/fastq tool kit.

Usage: fastx <command> <arguments>

Commands:
  head           head sequences
  sample         subsample sequences
  split          split fasta/fastq files.
  subseq         extract subsequences of fasta/fastq
```

License

[MIT License](https://github.com/Archieyoung/fastx/blob/master/LICENSE)
