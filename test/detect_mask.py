# Mask detect simulation 
# Detect the mask with all possible 1-SNP variants
# if the distance are not consistent with the variant, that could be the mask position

def catwalk_simulator(seq1, seq2, mask, use_mask):
    position = []
    if len(seq1) == len(seq2):
        for i, (x, y) in enumerate(zip(seq1, seq2)):
            if x != y:
                if i not in mask or not use_mask:
                    position.append(i)
    return position

if __name__ == '__main__':
    ref  =  'A' * 10
    mask = [2, 5]

    # create sequences
    seqs = list()
    for i, _ in enumerate(ref):
        seq = ref[:i] + 'C' + ref[i + 1:]
        seqs.append(seq)

    derived_mask = set()

    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            with_mask = catwalk_simulator(seq1, seq2, mask, True)
            without_mask = catwalk_simulator(seq1, seq2, mask, False)
            #print(f"{seq1} {seq2} {with_mask} {without_mask}")

            if with_mask != without_mask:
                #print(f"masked position: {set(with_mask).symmetric_difference(set(without_mask))}")
                derived_mask.add(list(set(with_mask).symmetric_difference(set(without_mask)))[0])

    print(derived_mask)





