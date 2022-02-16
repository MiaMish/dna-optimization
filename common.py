# Common functions
import random

STOP_CODONS = ['tag', 'taa', 'tga']
AA_TO_CODON = {'A': ['gct', 'gcc', 'gca', 'gcg'], 'C': ['tgt', 'tgc'], 'D': ['gat', 'gac'], 'E': ['gaa', 'gag'], 'F': ['ttt', 'ttc'], 'G': ['ggt', 'ggc', 'gga', 'ggg'], 'H': ['cat', 'cac'], 'I': ['att', 'atc', 'ata'], 'K': ['aag', 'aaa'], 'L': ['tta', 'ttg', 'ctt', 'ctc', 'ctg', 'cta'], 'M': ['atg'], 'N': ['aat', 'aac'], 'P': ['cct', 'ccc', 'cca', 'ccg'], 'Q': ['caa', 'cag'], 'R': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'], 'S': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'], 'T': ['act', 'acc', 'aca', 'acg'], 'V': ['gtt', 'gtc', 'gta', 'gtg'], 'W': ['tgg'], 'Y': ['tat', 'tac']}
CODON_TO_AA = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C', 'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C', 'tta': 'L', 'tca': 'S', 'ttg': 'L', 'tcg': 'S', 'tgg': 'W', 'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R', 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R', 'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R', 'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R', 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S', 'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S', 'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R', 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R', 'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G', 'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G', 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G', 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
DNA_SAMPLE = 'atgaccagccagattcgccagaactatagcaccgaagtggaagcggcggtgaaccgcctggtgaacctgcatctgcgcgcgagctatacctatctgagcctgggctttttttttgatcgcgatgatgtggcgctggaaggcgtgggccatttttttcgcgaactggcggaagaaaaacgcgaaggcgcggaacgcctgctggaatttcagaacgatcgcggcggccgcgcgctgtttcaggatgtgcagaaaccgagccaggatgaatggggcaaaacccaggaagcgatggaagcggcgctggcgatggaaaaaaacctgaaccaggcgctgctggatctgcatgcgctgggcagcgcgcgcaccgatccgcatctgtgcgattttctggaaagccattatctggataaagaagtgaaactgattaaaaaaatgggcaaccatctgaccaacctgcgccgcgtggcgggcccgcagccggcgcagaccggcgcgccgcagggcagcctgggcgaatatctgtttgaacgcctgaccctgaaacatgat'
MAX_DNA_SEQ_LEN_TO_LOG = 20  # always log only the first 20 chars of DNA seq to prevent long logs


def normalize_dna(dna_seq):
    normalized_dna_seq = dna_seq
    normalized_dna_seq = normalized_dna_seq.lower()
    normalized_dna_seq = normalized_dna_seq.replace('\n', ' ').replace(' ', '')
    return normalized_dna_seq


def validate_optimized_seq(optimized_seq, input_dna_seq, start_codon):
    optimized_aa = to_aa_representation(optimized_seq, start_codon)
    initial_aa = to_aa_representation(input_dna_seq, start_codon)
    if initial_aa != optimized_aa:
        print(f"Oh no! Something went wrong! AA representation of optimized seq is diff from the input seq!")
        return False
    return True


def analyze_optimized_seq(optimized_seq, input_dna_seq, start_codon):
    different_chars = 0
    codon_count_original = {"g": 0, "c": 0, "a": 0, "t": 0}
    codon_count_optimized = {"g": 0, "c": 0, "a": 0, "t": 0}
    for i in range(len(optimized_seq)):
        if optimized_seq[i] != input_dna_seq[i]:
            different_chars += 1
        codon_count_original[input_dna_seq[i]] += 1
        codon_count_optimized[optimized_seq[i]] += 1
    print(f"Replaced {different_chars} characters during optimization (DNA seq length is {len(optimized_seq)})")
    print(f"In original seq: "
          f"G - {(codon_count_original['g'] * 100) / len(input_dna_seq):.3f}% "
          f"C - {(codon_count_original['c'] * 100) / len(input_dna_seq):.3f}% "
          f"A - {(codon_count_original['a'] * 100) / len(input_dna_seq):.3f}% "
          f"T - {(codon_count_original['t'] * 100) / len(input_dna_seq):.3f}%")
    print(f"In optimized seq: "
          f"G - {(codon_count_optimized['g'] * 100) / len(optimized_seq):.3f}% "
          f"C - {(codon_count_optimized['c'] * 100) / len(optimized_seq):.3f}% "
          f"A - {(codon_count_optimized['a'] * 100) / len(optimized_seq):.3f}% "
          f"T - {(codon_count_optimized['t'] * 100) / len(optimized_seq):.3f}%")


def randomize_dna(start_codon):
    return ''.join(random.choice(["g", "c", "a", "t"]) for _ in range(999 + start_codon))


def validate_dna_seq(normalized_dna_seq, start_from=0):
    # TODO
    # Validate: only A, C, T, G
    # Validate length (divided by 3)
    return True


def to_aa_representation(codon_representation, start_from=0):
    index_in_codon_representation = start_from
    aa_representation = ""
    while index_in_codon_representation + 3 <= len(codon_representation):
        curr_codon = codon_representation[index_in_codon_representation:index_in_codon_representation + 3]
        if curr_codon in CODON_TO_AA:
            aa_representation = aa_representation + CODON_TO_AA[curr_codon]
        elif curr_codon in STOP_CODONS:
            # print(f"Translation to stop codons (as Z)")
            aa_representation = aa_representation + 'Z'
        index_in_codon_representation = index_in_codon_representation + 3
    return aa_representation


def trim_string_for_log(string):
    if string is None:
        return string
    if len(string) <= MAX_DNA_SEQ_LEN_TO_LOG:
        return string
    return f"{string[0:MAX_DNA_SEQ_LEN_TO_LOG]}..."


def get_closest_start_of_aa(first_codon, index_in_dna_seq):
    return int((index_in_dna_seq - first_codon) / 3) * 3 + first_codon


def try_to_remove(dna_sub_seq, start_index_in_relevant_sub_seq, to_remove):
    for index_in_relevant_sub_seq in range(start_index_in_relevant_sub_seq, len(to_remove) + start_index_in_relevant_sub_seq):
        replacement = try_to_replace(dna_sub_seq, index_in_relevant_sub_seq)
        if replacement != dna_sub_seq:
            print(f"[ Yay! ] Found replacement for {dna_sub_seq} without {to_remove}!")
            print(f"Replacement: [{replacement}]")
            return replacement
    print(f"[ Sorry! ] Could not find replacement for {dna_sub_seq} without {to_remove}!")
    return dna_sub_seq


def try_to_replace(dna_seq, index_to_replace):
    for candidate_codon in ['a', 'c', 'g', 't']:
        if dna_seq[index_to_replace] == candidate_codon:
            continue
        candidate_sub_seq_as_list = list(dna_seq)
        candidate_sub_seq_as_list[index_to_replace] = candidate_codon
        candidate_sub_seq = "".join(candidate_sub_seq_as_list)
        if to_aa_representation(candidate_sub_seq) == to_aa_representation(dna_seq):
            return candidate_sub_seq
    return dna_seq

