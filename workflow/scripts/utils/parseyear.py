#!/usr/bin/env python3

import re
import datetime

def parse_year(isolate):
  current_year = datetime.datetime.now().year + 1
  matches_long = re.findall('^[0-9]{8}', isolate)
  if len(matches_long) == 0:
    matches_short = re.findall('^[0-9]{7}', isolate)
    if len(matches_short) == 0:
      year = 'NA'
    else:
      year_raw = matches_short[0][0:3]
      year = year_raw[0] + '0' + year_raw[1:3]
  else:
    year = matches_long[0][0:4]
  if year != 'NA':
    if (int(year) < 2000) or (int(year) > current_year):
      year = 'NA'
  return year
