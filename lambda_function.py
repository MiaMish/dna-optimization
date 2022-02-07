import json
import os
import boto3


s3 = boto3.resource("s3")


def lambda_handler(event, context):
    print(f"Event=[{event}]\nContext=[{context}]")
    write_output()
    return {
        'statusCode': 200,
        'body': json.dumps('Hello from Lambda!')
    }


def write_output():
    encoded_string = "Hello".encode("utf-8")

    bucket_name = os.environ.get('S3_OUTPUT_BUCKET')
    s3_path = "hello.txt"

    print(f"Writing output to bucket=[{bucket_name}] path=[{s3_path}]")
    s3.Bucket(bucket_name).put_object(Key=s3_path, Body=encoded_string)

