#!/usr/bin/env python3

import subprocess
import argparse
import json

parser = argparse.ArgumentParser(description='Encode a genome using base64 and post this to PubMLST API for analysis.')

parser.add_argument('--genome', dest='genome', required=True, type=str, help='Genome assembly to type')
parser.add_argument('--api-url', dest='api', required=True, type=str, help='API URL of the typing scheme')

args = parser.parse_args()

def construct_command(genome, api_url):
  commandstring = '(echo -n \'{"base64":true,"sequence": "\'; base64 ' + genome + '; echo \'"}\') | curl -s -H "Content-Type: application/json" -X POST "' + api_url + '" -d @-'
  return commandstring

def post_pubmlst_api(command):
  api_result = subprocess.run(command, shell=True, capture_output=True)
  return api_result

def parse_pubmlst_api_output(api_result):
  cleaned_output = api_result.stdout.decode('utf-8').rstrip('\n')
  json_dict = json.loads(cleaned_output)
  print(json.dumps(json_dict, indent=2))

def main(args):
  commandstring = construct_command(args.genome, args.api)
  api_result = post_pubmlst_api(commandstring)
  parse_pubmlst_api_output(api_result)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Encode a genome using base64 and post this to PubMLST API for analysis.')
  parser.add_argument('--genome', dest='genome', required=True, type=str, help='Genome assembly to type')
  parser.add_argument('--api-url', dest='api', required=True, type=str, help='API URL of the typing scheme')
  args = parser.parse_args()

  main(args)
