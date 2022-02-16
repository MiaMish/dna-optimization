# Splicing Optimization
import math
import re

import requests

from common import *


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
