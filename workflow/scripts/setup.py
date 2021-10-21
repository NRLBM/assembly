#!/usr/bin/env python3

# Import necessary libraries
import subprocess
import os
import re


def print_version(tool):
  '''
  Run command with --version flag.

  Parameters
  ----------
  tool : str
    tool to run with --version flag

  Returns
  -------
  cleaned_output : str
    Parsed output from version command
  '''
  command = tool + ' --version'
  output = subprocess.run([command], capture_output=True, shell=True)
  cleaned_output = output.stdout.decode('utf-8').rstrip('\n')
  return cleaned_output

def check_version():
  '''
  Check versions for dependencies.

  Loops through three tools and parses version commands.

  Parameters
  ----------
  None

  Returns
  -------
  None

  Raises
  ------
  ValueError
    If the --version command does not produce STDOUT, suggesting the tool is not installed.

  '''
  list_dependencies = ['snakemake', 'sinfo', 'conda']
  for tool in list_dependencies:
    print('Checking version for ' + tool + ':')
    output = print_version(tool)
    if output == '':
      errstring = tool + ' does not seem to be available currently. Please install and try again.'
      raise ValueError(errstring)
    else:
      print('OK: ' + output)

def check_bashrc():
  '''
  Check whether bin directory is in PATH.

  Checks os.environ whether the bin directory is already in PATH and adds an export statement to .bashrc if this is not the case.

  Parameters
  ----------
  None

  Returns
  -------
  None

  '''
  print('Checking whether assembly/bin is in path')
  cwd = os.getcwd()
  bin_path = os.path.join(cwd, 'bin')
  export_path_line = 'export PATH=$PATH:' + bin_path
  PATH = os.environ['PATH'].split(':')
  if (bin_path in PATH):
    print('assembly/bin is already in PATH')
  else:
    bashrc_path = os.path.expanduser('~/.bashrc')
    print('Adding export statement: "' + export_path_line  + '" to ' + bashrc_path)
    f = open(bashrc_path, 'a')
    f.write('\n')
    f.write(export_path_line)
    f.close()

def check_email_file(person):
  '''
  Check whether the .email_* files contain valid email addresses.

  Checks whether .email_* files exist. If a file does not exist, do nothing and continue to next check. If the file exists, open content and check whether it contains a single line with a single valid email address.

  Parameters
  ----------
  person : str
    Person to parse email address for. Usually either "user" or "bioinformatician".

  Returns
  -------
  True
    If the email file is valid.

  Raises
  ------
  ValueError
    If the email file contains multiple lines.
  ValueError
    If the email address in the email file is invalid.

  '''
  email_file_path = '.email_' + str(person)
  if os.path.exists(email_file_path):
    f = open(email_file_path)
    email_raw = f.readlines()
    f.close()
    email = email_raw[0].rstrip('\n')
    if len(email_raw) != 1:
      errstring = 'Email file ' + str(email_file_path) + ' contains multiple lines (should be a single line with a single email address).'
      raise ValueError(errstring)
    regex = '[^@]+@[^@]+\.[^@]+'
    if not re.match(regex, email):
      errstring = str(email) + ' does not seem to be a valid email address. Please run the script again or enter manually'
      raise ValueError(errstring)
    else:
      print(email + ' is saved in ' + str(email_file_path) + ' a valid email address.')
      return True

def input_email(person):
  '''
  Read input to write email address to email file.

  Parameters
  ----------
  person : str
    Person to parse email address for. Usually either "user" or "bioinformatician".

  Returns
  -------
  None

  Raises
  ------
  ValueError
    If the entered email address is not valid.

  '''
  email = input('Please enter email address for ' + str(person) + ': ')
  regex = '[^@]+@[^@]+\.[^@]+'
  if not re.match(regex, email):
    errstring = str(email) + ' does not seem to be a valid email address. Please run the script again or enter manually'
    raise ValueError(errstring)
  else:
    email_file_path = '.email_' + str(person)
    print('Writing ' + str(email) + ' to hidden file ' + str(email_file_path))
    f = open(email_file_path, 'w')
    f.write(email)
    f.close()

def check_email():
  '''
  Combine email check functions for multiple persons.

  Parameters
  ----------
  None

  Returns
  -------
  None

  Raises
  ------
  None

  '''
  list_persons = ['user', 'bioinformatician']
  for person in list_persons:
    checkpoint = check_email_file(person)
    if not checkpoint:
      input_email(person)

def ln_analyse_runs():
  '''
  Create a symlink in PATH to analyse_runs.py

  Parameters
  ----------
  None

  Returns
  -------
  None

  Raises
  ------
  None

  '''
  full_path = os.path.abspath('workflow/scripts/analyse_runs.py')
  executable_target = 'bin/analyse_runs'

  subprocess.run(['ln', '-s', full_path, executable_target])
  subprocess.run(['chmod', 'u+x', executable_target])


if __name__ == '__main__':
  check_version()
  check_bashrc()
  check_email()
  ln_analyse_runs()
  print('All checks passed')
