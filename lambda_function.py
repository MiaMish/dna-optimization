import os
import boto3
import os

import boto3

s3 = boto3.resource("s3")

# DEFAULT SETTINGS AND CONSTANTS

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







# def is_local_env():
#     return 'AWS_EXECUTION_ENV' not in os.environ.keys()
#
#
# if is_local_env():
#     cli_main()
