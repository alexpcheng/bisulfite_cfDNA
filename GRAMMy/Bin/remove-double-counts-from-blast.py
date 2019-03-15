#!/usr/bin/env python
# Selects only one read from a pair from BLAST results
# Usage: ./<this>.py <blast file> <orphan list>

import sys
from collections import defaultdict

blast1 = defaultdict(lambda: defaultdict(list))
blast2 = defaultdict(lambda: defaultdict(list))
orphans = []
f = open(sys.argv[2])
for line in f:
  orphans.append(line.strip())

def rm_suffix(name):
  return name[:-2]

def get_suffix(name):
  return name[-1]

f = open(sys.argv[1])
for line in f:
  parts = line.split()
  query = parts[0]
  target = parts[1]
  suffix = get_suffix(query)
  if suffix == "C":
    print line,
  else:
    if query in orphans:
      print line,
    elif suffix == "1":
      blast1[target][rm_suffix(query)].append(line)
    elif suffix == "2":
      blast2[target][rm_suffix(query)].append(line)

f.close()

for target in blast1:
  for query in blast1[target]:
    if len(blast2[target][query]) >= 1:
      print blast1[target][query][0],
