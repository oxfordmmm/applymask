#! /usr/bin/env python3

import math
import gzip
import argparse
import pathlib

# fasta format;
# >header
# 000101010010101010 ...
# where 1 means mask
# 0 means don't mask
def load_mask_fasta(mask_filepath):
    with open(mask_filepath) as f:
        lines = f.readlines()
    mask = "".join([line.strip() for line in lines[1:]])
    return mask

# positiion format
# one position per line
# 0 index
def load_mask_positon(mask_filepath, length):
    with open(mask_filepath) as f:
        lines = f.readlines()
    return load_mask_position_aux(lines, length)

def load_mask_position_aux(lines, length):
    masked_positions = set()
    for line in lines:
        try:
            masked_positions.add(int(line))
        except ValueError:
            print(f"masked position '{line}' is not an integer")
    mask = list()
    for i in range(0, length):
        if i in masked_positions:
            mask.append('1')
        else:
            mask.append('0')
    return mask
    
# range format
# tsv range per line
# 2000 3000 
# 12000 130000
def load_mask_range(mask_filepath, length):
    with open(mask_filepath) as f:
        lines = f.readlines()
    return load_mask_range_aux(lines, length)   

def load_mask_range_aux(lines, length):
    masked_ranges = list()
    for line in lines:
        try:
            begin, end = line.strip().split('\t')
            masked_ranges.append((int(begin), int(end)))
        except ValueError:
            print(f"masked range '{line}' couldn't be parsed")

    mask = list()
    for i in range(0, length):
        mask.append('0')

    for begin, end in masked_ranges:
        if begin == end:
            mask[begin] = '1'
        else:  
            for i in range(begin, end+1):
                mask[i] = '1'
    return mask

# [1,2,3,4,5,7] -> "1-5,7"
def lst_to_range_str(lst):
    ret = list()
    start = None
    end = None
    for elem in lst:
        if start is None:
            start = elem
            end = elem
        if elem == end or elem == end + 1:
            end = elem
        else:
            if start == end:
                ret.append((start, start))
            else:
                ret.append((start, end))
            start = elem
            end = elem
    if start:
        if start == end:
            ret.append((start, start))
    return ret

def get_mask_ranges(fasta_mask):
    mask_lst = list()
    for i, x in enumerate(fasta_mask):
        if x == '1':
            mask_lst.append(i)
    mask_str = lst_to_range_str(mask_lst)
    mask_ranges = [f"{x}\t{y}" for (x, y) in mask_str]
    return "\n".join(mask_ranges)

def load_fasta(fasta_filepath):
    with open(fasta_filepath) as f:
        lines = f.readlines()
    header = lines[0].strip()
    chunk_len = len(lines[1]) - 1 #remove the new line
    sequence = "".join([line.strip() for line in lines[1:]])
    return (header, sequence, chunk_len)

def load_fasta_gzip(fasta_filepath):
    with gzip.open(fasta_filepath, "rb") as f:
        lines = f.read().decode().split('\n')
    header = lines[0].strip()
    chunk_len = len(lines[1]) # no new line in gzip
    sequence = "".join([line.strip() for line in lines[1:]])
    return (header, sequence, chunk_len)

def string_insert_newlines(in_str, chunk_len):
    ret = list()
    for i in range(0, math.ceil(len(in_str) / chunk_len)):
        chunk = in_str[i * chunk_len:(i + 1) * chunk_len]
        ret.append(chunk)
    new_str = "\n".join(ret)
    return new_str

def save_fasta(fasta_filepath, header, sequence, chunk_len):
    with open(fasta_filepath, "w") as f:
        f.write(header + '\n')
        f.write(string_insert_newlines(sequence, chunk_len))
        print(f'file {fasta_filepath} created.')

def save_fasta_gzip(fasta_filepath, header, sequence, chunk_len):
    with gzip.open(fasta_filepath, "wb") as f:
        f.write(bytearray(header + '\n', 'utf-8'))
        f.write(bytearray(string_insert_newlines(sequence, chunk_len), 'utf-8'))
        print(f'file {fasta_filepath} created.')

def apply_mask(mask, sequence):
    new_sequence = list()
    for i, _ in enumerate(sequence):
        if mask[i] == '1':
            new_sequence.append('N')
        else:
            new_sequence.append(sequence[i])
    new_sequence_str = "".join(new_sequence)
    return new_sequence_str

def arg_is_true(arg):
    return arg in ["true", "yes", "oui"]

def main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges):
    
    #Load mask file, transform to fasta format
    if arg_is_true(use_gzip):
        (fasta_header, fasta_sequence, fasta_chunk_len) = load_fasta_gzip(fasta_filepath)
    else:
        (fasta_header, fasta_sequence, fasta_chunk_len) = load_fasta(fasta_filepath)

    if mask_format == "fasta":
        mask = load_mask_fasta(mask_filepath)
    elif mask_format == "position":
        mask = load_mask_positon(mask_filepath, len(fasta_sequence))
    elif mask_format == "range":
        mask = load_mask_range(mask_filepath, len(fasta_sequence))
    else:
        print(f"unknown mask format: {mask_format}")
        print("expected one of: fasta, position, range")
        return

    new_sequence = apply_mask(mask, fasta_sequence)

    fasta_filename = pathlib.Path(fasta_filepath).stem
    if arg_is_true(use_gzip):
        new_fasta_filepath = fasta_filename + '.masked.gz'
    else:
        new_fasta_filepath = fasta_filename + '.masked.fasta'

    if arg_is_true(use_gzip):
        save_fasta_gzip(new_fasta_filepath, fasta_header, new_sequence, fasta_chunk_len)
    else:
        save_fasta(new_fasta_filepath, fasta_header, new_sequence, fasta_chunk_len)

    if arg_is_true(print_mask_ranges):
        print(get_mask_ranges(mask))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mask_filepath")
    p.add_argument("mask_format")
    p.add_argument("fasta_filepath")
    p.add_argument("use_gzip")
    p.add_argument("print_mask_ranges")
    args = p.parse_args()
    main(args.mask_filepath, args.mask_format, args.fasta_filepath, args.use_gzip, args.print_mask_ranges)
