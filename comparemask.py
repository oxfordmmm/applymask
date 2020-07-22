# compare the masks
# return the differen positions
# python3 compare.py
def load_position_mask(mask_filepath):
    with open(mask_filepath) as f:
        lines = f.readlines()
    masked_positions = []
    for line in lines:
        try:
            masked_positions.append(int(line))
        except ValueError:
            print(f"masked position '{line}' is not an integer")
    return masked_positions

def compare_position_mask(mask1, mask2):
    mask1_set = set(mask1)
    mask2_set = set(mask2)
    diff1 = mask1_set - mask2_set
    diff2 = mask2_set - mask1_set
    diff_all = diff1.union(diff2)
    return sorted(list(diff1)), sorted(list(diff2)), sorted(list(diff_all))

def fasta_to_posistions(fasta_mask_file):
    with open(fasta_mask_file) as f:
        orig = f.readlines()
        orig = [x.strip() for x in orig[1:]]
        orig = ''.join(orig)
    
    out = []
    for i,x in enumerate(orig):
        if x == '1':
            out.append(i)
    return out

def write_position_to_file(mask, output_file):
    mask_str = [str(i) for i in mask]
    mask_str_lines = '\n'.join(mask_str)
    with open(output_file, "w") as f:
        f.write(mask_str_lines)
    return mask_str_lines

if __name__ == '__main__':
    inputfile1 = "tb/TB-exclude.txt"
    masked_positions1 = load_position_mask(inputfile1)

    inputfile2 = "tb/TB-exclude-adaptive.txt"
    masked_positions2 = load_position_mask(inputfile2)

    diff1_2, diff2_1, diff_all = compare_position_mask(masked_positions1, masked_positions2)
    
    print(len(diff1_2))
    print(len(diff2_1))
    print(len(diff_all))

    inputfile3 = "tb/NC_000962_2_repmask.array"
    masked_positions3 = fasta_to_posistions(inputfile3)
    result = write_position_to_file(masked_positions3, 'tb/repmask.txt')
