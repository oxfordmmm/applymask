# applymask

## Description

Apply mask to a given fasta file

## Supported formats

- fasta format e.g.
```
>mask header
010101010101
```
where 1 means to mask

- wyllie format e.g.
```
3
5
10
```
where the numbers are the positions to mask (0 indexed)

## Command line use

```
$ python3 applymask.py
usage: applymask.py [-h] mask_filepath mask_format fasta_filepath use_gzip
```