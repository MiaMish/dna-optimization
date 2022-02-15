import json
import os
import math
import boto3
import requests
import re
import click
import random
from datetime import datetime


s3 = boto3.resource("s3")

# DEFAULT SETTINGS AND CONSTANTS
STOP_CODONS = ['tag', 'taa', 'tga']
AA_TO_CODON = {'A': ['gct', 'gcc', 'gca', 'gcg'], 'C': ['tgt', 'tgc'], 'D': ['gat', 'gac'], 'E': ['gaa', 'gag'], 'F': ['ttt', 'ttc'], 'G': ['ggt', 'ggc', 'gga', 'ggg'], 'H': ['cat', 'cac'], 'I': ['att', 'atc', 'ata'], 'K': ['aag', 'aaa'], 'L': ['tta', 'ttg', 'ctt', 'ctc', 'ctg', 'cta'], 'M': ['atg'], 'N': ['aat', 'aac'], 'P': ['cct', 'ccc', 'cca', 'ccg'], 'Q': ['caa', 'cag'], 'R': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'], 'S': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'], 'T': ['act', 'acc', 'aca', 'acg'], 'V': ['gtt', 'gtc', 'gta', 'gtg'], 'W': ['tgg'], 'Y': ['tat', 'tac']}
CODON_TO_AA = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C', 'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C', 'tta': 'L', 'tca': 'S', 'ttg': 'L', 'tcg': 'S', 'tgg': 'W', 'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R', 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R', 'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R', 'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R', 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S', 'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S', 'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R', 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R', 'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G', 'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G', 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G', 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
DNA_SAMPLE = 'atgaccagccagattcgccagaactatagcaccgaagtggaagcggcggtgaaccgcctggtgaacctgcatctgcgcgcgagctatacctatctgagcctgggctttttttttgatcgcgatgatgtggcgctggaaggcgtgggccatttttttcgcgaactggcggaagaaaaacgcgaaggcgcggaacgcctgctggaatttcagaacgatcgcggcggccgcgcgctgtttcaggatgtgcagaaaccgagccaggatgaatggggcaaaacccaggaagcgatggaagcggcgctggcgatggaaaaaaacctgaaccaggcgctgctggatctgcatgcgctgggcagcgcgcgcaccgatccgcatctgtgcgattttctggaaagccattatctggataaagaagtgaaactgattaaaaaaatgggcaaccatctgaccaacctgcgccgcgtggcgggcccgcagccggcgcagaccggcgcgccgcagggcagcctgggcgaatatctgtttgaacgcctgaccctgaaacatgat'
MAX_DNA_SEQ_LEN_TO_LOG = 20  # always log only the first 20 chars of DNA seq to prevent long logs


# IO functions for lambda


def read_input(event):
    bucket_name = event.get('Records', [{}])[0].get("s3", {}).get("bucket", {}).get("name")
    s3_path = event.get('Records', [{}])[0].get("s3", {}).get("object", {}).get("key")
    print(f"Reading input from bucket=[{bucket_name}] path=[{s3_path}]")
    s3_object = s3.Object(bucket_name, s3_path).get()
    print(s3_object)
    return s3_object['Body'].read().decode('utf-8'), f"s3://{bucket_name}/{s3_path}"


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


# Pre processing


def normalize_dna(dna_seq):
    normalized_dna_seq = dna_seq
    normalized_dna_seq = normalized_dna_seq.lower()
    normalized_dna_seq = normalized_dna_seq.replace('\n', ' ').replace(' ', '')
    return normalized_dna_seq


# Common functions


def validate_dna_seq(normalized_dna_seq, start_from=0):
    # TODO
    # Validate: only A, C, T, G
    # Validate length (divided by 3)
    return True


def to_aa_representation(codon_representation, start_from=0):
    index_in_codon_representation = start_from
    aa_representation = ""
    while index_in_codon_representation + 3 < len(codon_representation):
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


# GGGG/CCCC Removal


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


# Splicing Optimization


def req_splicing(dna_seq, first_codon=0):
    after_cleanup = dna_seq

    url = "https://www.fruitfly.org/cgi-bin/seq_tools/splice.pl"
    payload = f"organism=human&which=both&reverse=no&min_donor=0.4&min_acc=0.4&text={dna_seq}"
    headers = {
        "Content-Type": "application/x-www-form-urlencoded"
    }
    # print(f"Creating request: url: [{url}] headers: [{headers}] payload: [{payload}]")
    response = requests.request("POST", url, headers=headers, data=payload)
    if response.status_code != 200:
        print(f"Got response: statusCode: [{response.status_code}] body: [{response.text}]")
        print(f"Request failed")
        return after_cleanup

    rows = re.findall(
        r'(<br>[\s]+([\d.]+)[\s\d.]+([gatc]+)<font size="\+2">([gatc]+)</font><font size="\+2">([gatc]+)</font>([gatc]+)+<br>)',
        response.text)
    print(f"Found {len(rows)} splice sites")
    for row in rows:
        # print(f"Parsing row [{row}]")
        start_index = int(row[1])
        intron = row[2]
        left_boundary = row[3]
        right_boundary = row[4]
        exon = row[5]

        print(f"Trying to handle splice site that starts in index [{start_index}]: intron [{intron}] left_boundary [{left_boundary}] right_boundary [{right_boundary}] exon [{exon}]")
        start_sub_seq = get_closest_start_of_aa(first_codon, start_index)
        end_sub_seq = int(math.ceil((start_sub_seq + len(intron) + len(exon) + 2) / 3.0) * 3)
        relevant_sub_seq = dna_seq[start_sub_seq:end_sub_seq]
        index_in_relevant_sub_seq = start_index + len(intron) + 1 - start_sub_seq
        print(f"Trying to remove left boundary of splice site ({left_boundary} in index {index_in_relevant_sub_seq}) from DNA sub seq {relevant_sub_seq} (index in DNA seq [{start_sub_seq}, {end_sub_seq}])")
        after_cleanup_sub_seq = try_to_replace(relevant_sub_seq, index_in_relevant_sub_seq)
        if after_cleanup_sub_seq != relevant_sub_seq:
            print(f"[ Yay! ] Found replacement for {relevant_sub_seq} without left boundary of splice site!")
            print(f"Replacement: [{after_cleanup_sub_seq}]")
            after_cleanup = after_cleanup[:start_sub_seq] + after_cleanup_sub_seq + after_cleanup[end_sub_seq:]
            continue
        after_cleanup_sub_seq = try_to_replace(relevant_sub_seq, index_in_relevant_sub_seq + 1)
        if after_cleanup_sub_seq != relevant_sub_seq:
            print(f"[ Yay! ] Found replacement for {relevant_sub_seq} without right boundary of splice site!")
            print(f"Replacement: [{after_cleanup_sub_seq}]")
            after_cleanup = after_cleanup[:start_sub_seq] + after_cleanup_sub_seq + after_cleanup[end_sub_seq:]
            continue
        print(f"[ Sorry! ] Could not find replacement for {relevant_sub_seq} without left or right boundary!")
    return after_cleanup


# Lambda handler

# init S3 client outside lambda to prevent cold start
# s3 = boto3.resource("s3")


# def lambda_handler(event, context):
#     response = {}
#     print(f"Event=[{event}]\nContext=[{context}]")
#     response['input'], response['input_source'] = read_input(event)
#     print(response)
#     response['normalized_input'] = normalize_dna(response['input'])
#     is_valid_input, error_reason = validate_dna_seq(response['normalized_input'])
#     if not is_valid_input:
#         response['error'] = 'invalid input'
#         response['error_reason'] = error_reason
#         return response
#     response['aa_representation'] = to_aa_representation(response['normalized_input'])
#     response['steps'] = []
#     response['steps'].append({
#         'name': 'GGGG removal',
#         'input': response['normalized_input'],
#         'output': remove_codons(response['normalized_input'], 'GGGG')
#     })
#
#     response['output_path'] = write_output()
#     return response


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


# Run locally

@click.command()
@click.option('--dna-seq', help='DNA sequence', prompt='DNA seq (leave empty for random DNA generation)', default="")
@click.option('--optimization', type=click.Choice(['gggg_removal', 'cccc_removal', 'splicing']), prompt='Optimization type')
@click.option('--start-codon', type=click.Choice(['0', '1', '2']), prompt='Start codon', default='0')
@click.option('--output-file', prompt='Output file path', default='out.txt')
@click.option('--seq-name', prompt='Name for seq', default='')
def cli_main(dna_seq, optimization, start_codon, output_file, seq_name):
    start_codon_as_int = int(start_codon)
    if dna_seq == "":
        print("Randomizing DNA seq")
        dna_seq = randomize_dna(start_codon_as_int)
    print(f"Using input DNA seq [{trim_string_for_log(dna_seq)}] (length: {len(dna_seq)}) with start codon {start_codon_as_int}")
    normalized_input = normalize_dna(dna_seq)
    is_valid_input = validate_dna_seq(normalized_input)
    if not is_valid_input:
        return
    if optimization == 'gggg_removal':
        optimized_seq = remove_gggg_codons(normalized_input, start_codon_as_int)
    elif optimization == 'cccc_removal':
        optimized_seq = remove_cccc_codons(normalized_input, start_codon_as_int)
    elif optimization == 'splicing':
        optimized_seq = req_splicing(normalized_input, start_codon_as_int)
    is_optimized_valid = validate_optimized_seq(optimized_seq, normalized_input, start_codon_as_int)
    if not is_optimized_valid:
        return
    analyze_optimized_seq(optimized_seq, normalized_input, start_codon_as_int)
    if seq_name == "":
        seq_name = f"DNA seq {datetime.today().strftime('%Y-%m-%d-%H:%M:%S')}"
        print(f"Seq name was not provided, using auto generated name: {seq_name}")
    print(f"Writing results to {output_file}")
    with open(output_file, 'a') as fout:
        fout.write(f"> Original {seq_name}\n")
        fout.write(f"{dna_seq}\n\n")
        fout.write(f"> Optimized {seq_name}\n")
        fout.write(f"{optimized_seq}\n\n")
    return


def is_local_env():
    return 'AWS_EXECUTION_ENV' not in os.environ.keys()


if is_local_env():
    cli_main()
