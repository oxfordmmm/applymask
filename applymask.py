import math
import gzip
import argparse
import pathlib

# format 1, similar to fasta;
# >header
# 000101010010101010 ...
# where 1 means mask
# 0 means don't mask
def load_mask_format1(mask_filepath):
    with open(mask_filepath) as f:
        lines = f.readlines()
    mask = "".join([line.strip() for line in lines[1:]])
    return mask

# dwyllie format
# one position per line
def load_mask_format2(mask_filepath, length):
    with open(mask_filepath) as f:
        lines = f.readlines()
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
def load_mask_format3(mask_filepath, length):
    with open(mask_filepath) as f:
        lines = f.readlines()
    masked_ranges = list()
    for line in lines:
        try:
            begin, end = line.split('\t')
            masked_ranges.append((int(begin), int(end)))
        except ValueError:
            print(f"masked range '{line}' couldn't be parsed")
    mask = list()
    for i in range(0, length):
        mask.append('0')
    for begin, end in masked_ranges:
        for i in range(begin, end):
            mask[i] = '1'
    return mask

def load_fasta(fasta_filepath):
    with open(fasta_filepath) as f:
        lines = f.readlines()
    header = lines[0]
    chunk_len = len(lines[1]) - 1
    sequence = "".join([line.strip() for line in lines[1:]])
    return (header, sequence, chunk_len)

def load_fasta_gzip(fasta_filepath):
    with gzip.open(fasta_filepath, "rb") as f:
        lines = f.read().decode().split('\n')
    header = lines[0]
    chunk_len = len(lines[1]) - 1
    sequence = "".join([line.strip() for line in lines[1:]])
    return (header, sequence, chunk_len)

def string_insert_newlines(in_str, chunk_len):
    ret = list()
    for i in range(0, math.ceil(len(in_str) / chunk_len)):
        chunk = in_str[i * chunk_len:(i + 1) * chunk_len]
        ret.append(chunk)
        #print(i, i * chunk_len, (i + 1) * chunk_len, chunk)
    #print(ret)
    new_str = "\n".join(ret)
    return new_str

def save_fasta(fasta_filepath, header, sequence, chunk_len):
    with open(fasta_filepath, "w") as f:
        f.write(header)
        f.write(string_insert_newlines(sequence, chunk_len))
        f.write("\n")

def save_fasta_gzip(fasta_filepath, header, sequence, chunk_len):
    with gzip.open(fasta_filepath, "wb") as f:
        f.write(bytearray(header, 'utf-8'))
        f.write(bytearray(string_insert_newlines(sequence, chunk_len), 'utf-8'))
        f.write(bytearray("\n", 'utf-8'))

def apply_mask(mask, sequence):
    new_sequence = list()
    for i, _ in enumerate(sequence):
        if mask[i] == '1':
            new_sequence.append('N')
        else:
            new_sequence.append(sequence[i])
    new_sequence_str = "".join(new_sequence)
    return new_sequence_str

def main():
    p = argparse.ArgumentParser()
    p.add_argument("mask_filepath")
    p.add_argument("mask_format")
    p.add_argument("fasta_filepath")
    p.add_argument("use_gzip")
    args = p.parse_args()

    if args.use_gzip == "true":
        (fasta_header, fasta_sequence, fasta_chunk_len) = load_fasta_gzip(args.fasta_filepath)
    else:
        (fasta_header, fasta_sequence, fasta_chunk_len) = load_fasta(args.fasta_filepath)

    if args.mask_format == "fasta":
        mask = load_mask_format1(args.mask_filepath)
    elif args.mask_format == "dwyllie":
        mask = load_mask_format2(args.mask_filepath, len(fasta_sequence))
    elif args.mask_format == "ranges":
        mask = load_mask_format3(args.mask_filepath, len(fasta_sequence))
    else:
        print(f"unknown mask format: {args.mask_format}")
        print("expected one of: fasta, dwyllie, ranges")
        return

    new_sequence = apply_mask(mask, fasta_sequence)

    fasta_filename = pathlib.Path(args.fasta_filepath).stem
    if args.use_gzip == "true":
        new_fasta_filepath = fasta_filename + '.masked.gz'
    else:
        new_fasta_filepath = fasta_filename + '.masked.fasta'


    if args.use_gzip == "true":
        save_fasta_gzip(new_fasta_filepath, fasta_header, new_sequence, fasta_chunk_len)
    else:
        save_fasta(new_fasta_filepath, fasta_header, new_sequence, fasta_chunk_len)

if __name__ == "__main__":
    main()
