import json
import os

import boto3

s3 = boto3.resource("s3")

# DEFAULT SETTINGS AND CONSTANTS
STOP_CODONS = ['TAG', 'TAA', 'TGA']

AA_TO_CODON = {'A': ['GCT', 'GCC', 'GCA', 'GCG'],
               'C': ['TGT', 'TGC'],
               'D': ['GAT', 'GAC'],
               'E': ['GAA', 'GAG'],
               'F': ['TTT', 'TTC'],
               'G': ['GGT', 'GGC', 'GGA', 'GGG'],
               'H': ['CAT', 'CAC'],
               'I': ['ATT', 'ATC', 'ATA'],
               'K': ['AAG', 'AAA'],
               'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTG', 'CTA'],
               'M': ['ATG'],
               'N': ['AAT', 'AAC'],
               'P': ['CCT', 'CCC', 'CCA', 'CCG'],
               'Q': ['CAA', 'CAG'],
               'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
               'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
               'T': ['ACT', 'ACC', 'ACA', 'ACG'],
               'V': ['GTT', 'GTC', 'GTA', 'GTG'],
               'W': ['TGG'],
               'Y': ['TAT', 'TAC']}

CODON_TO_AA = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
               'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
               'TTA': 'L', 'TCA': 'S', 'TTG': 'L', 'TCG': 'S',
               'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H',
               'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H',
               'CGC': 'R', 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q',
               'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q',
               'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N',
               'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N',
               'AGC': 'S', 'ATA': 'I', 'ACA': 'T', 'AAA': 'K',
               'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K',
               'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D',
               'GGT': 'G', 'GTC': 'V', 'GCC': 'A', 'GAC': 'D',
               'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E',
               'GGA': 'G', 'GTG': 'V', 'GCG': 'A', 'GAG': 'E',
               'GGG': 'G'}

DNA_SAMPLE = 'atgaccagccagattcgccagaactatagcaccgaagtggaagcggcggtgaaccgcctggtgaacctgcatctgcgcgcgagctatacctatctgagcctgggctttttttttgatcgcgatgatgtggcgctggaaggcgtgggccatttttttcgcgaactggcggaagaaaaacgcgaaggcgcggaacgcctgctggaatttcagaacgatcgcggcggccgcgcgctgtttcaggatgtgcagaaaccgagccaggatgaatggggcaaaacccaggaagcgatggaagcggcgctggcgatggaaaaaaacctgaaccaggcgctgctggatctgcatgcgctgggcagcgcgcgcaccgatccgcatctgtgcgattttctggaaagccattatctggataaagaagtgaaactgattaaaaaaatgggcaaccatctgaccaacctgcgccgcgtggcgggcccgcagccggcgcagaccggcgcgccgcagggcagcctgggcgaatatctgtttgaacgcctgaccctgaaacatgat'
s3 = boto3.resource("s3")


def read_input(event):
    if is_local_env():
        return DNA_SAMPLE, "hard coded sample"
    else:
        bucket_name = event.get('Records', [{}])[0].get("s3", {}).get("bucket", {}).get("name")
        s3_path = event.get('Records', [{}])[0].get("s3", {}).get("object", {}).get("key")
        print(f"Reading input from bucket=[{bucket_name}] path=[{s3_path}]")
        s3_object = s3.Object(bucket_name, s3_path).get()
        print(s3_object)
        return s3_object['Body'].read().decode('utf-8'), f"s3://{bucket_name}/{s3_path}"


def normalize_dna(dna_seq):
    normalized_dna_seq = dna_seq
    normalized_dna_seq = normalized_dna_seq.upper()
    normalized_dna_seq = normalized_dna_seq.replace('\n', ' ').replace('\r', '')
    return normalized_dna_seq


def validate_dna_seq(normalized_dna_seq, start_from=0):
    # TODO
    # Validate: only A, C, T, G
    # Validate length (divided by 3)
    return True, None


def to_aa_representation(codon_representation, start_from=0):
    index_in_codon_representation = start_from
    aa_representation = ""
    while index_in_codon_representation + 3 < len(codon_representation):
        curr_codon = codon_representation[index_in_codon_representation:index_in_codon_representation + 3]
        aa_representation = aa_representation + CODON_TO_AA[curr_codon]
        index_in_codon_representation = index_in_codon_representation + 3
    return aa_representation


def remove_codons(dna_seq, to_remove, first_codon=0):
    # TODO impl
    # after_cleanup = dna_seq
    # aa_representation = to_aa_representation(dna_seq, start_from=first_codon)
    # search_from_index = 0
    # while True:
    #     index_in_dna_seq = after_cleanup.find(to_remove, search_from_index)
    #     if index_in_dna_seq == -1:
    #         break
    #     start_from  = ((index_in_dna_seq - first_codon) / 3) * 3 + first_codon
    #     end_after = start_from + ((len(to_remove) / 3) + 1) * 3
    #     relevant_sub_seq = dna_seq[start_from, end_after]
    #     for i in range(0, len(to_remove)):
    return dna_seq


def lambda_handler(event, context):
    response = {}
    print(f"Event=[{event}]\nContext=[{context}]")
    response['input'], response['input_source'] = read_input(event)
    print(response)
    response['normalized_input'] = normalize_dna(response['input'])
    is_valid_input, error_reason = validate_dna_seq(response['normalized_input'])
    if not is_valid_input:
        response['error'] = 'invalid input'
        response['error_reason'] = error_reason
        return response
    response['aa_representation'] = to_aa_representation(response['normalized_input'])
    response['steps'] = []
    response['steps'].append({
        'name': 'GGGG removal',
        'input': response['normalized_input'],
        'output': remove_codons(response['normalized_input'], 'GGGG')
    })

    response['output_path'] = write_output()
    return response


def write_output():
    if 'S3_OUTPUT_BUCKET' not in os.environ.keys():
        print(f"Not writing output to file, missing S3_OUTPUT_BUCKET config")
        return 'no output file'

    encoded_string = "Hello".encode("utf-8")

    bucket_name = os.environ.get('S3_OUTPUT_BUCKET')
    s3_path = "hello.txt"

    print(f"Writing output to bucket=[{bucket_name}] path=[{s3_path}]")
    s3.Bucket(bucket_name).put_object(Key=s3_path, Body=encoded_string)
    return f"s3://{bucket_name}/{s3_path}"


def is_local_env():
    return 'AWS_EXECUTION_ENV' not in os.environ.keys()


if is_local_env():
    print(json.dumps(lambda_handler(None, None), indent=4))

