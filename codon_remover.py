import math
from common import *


def _remove_codons(dna_seq, seq_to_remove, first_codon):
    after_cleanup = dna_seq
    search_from_index = 0
    print(f"Trying to remove [{seq_to_remove}] from DNA seq [{trim_string_for_log(dna_seq)}] (length: {len(dna_seq)})")
    while True:
        index_in_dna_seq = after_cleanup.find(seq_to_remove, search_from_index)
        if index_in_dna_seq == -1:
            break
        if search_from_index == index_in_dna_seq:
            break
        print(f"Found [{seq_to_remove}] in index: [{index_in_dna_seq}] (search from index: {search_from_index})")
        start_from = get_closest_start_of_aa(first_codon, index_in_dna_seq)
        end_after = int(math.ceil((start_from + len(seq_to_remove)) / 3.0) * 3)
        relevant_sub_seq = dna_seq[start_from:end_after]
        index_in_relevant_sub_seq = index_in_dna_seq - start_from
        print(f"Trying to remove [{seq_to_remove}] from DNA sub seq {relevant_sub_seq} (index in DNA seq [{start_from}, {end_after}])")
        after_cleanup_sub_seq = try_to_remove(relevant_sub_seq, index_in_relevant_sub_seq, seq_to_remove)
        # print(f"after_cleanup_sub_seq={after_cleanup_sub_seq}")
        after_cleanup = after_cleanup[:start_from] + after_cleanup_sub_seq + after_cleanup[end_after:]
        search_from_index = index_in_dna_seq
    print(f"Finished removing [{seq_to_remove}] from DNA seq.")
    return after_cleanup


def remove_gggg_codons(dna_seq, first_codon=0):
    return _remove_codons(dna_seq, "gggg", first_codon)


def remove_cccc_codons(dna_seq, first_codon=0):
    return _remove_codons(dna_seq, "cccc", first_codon)
