#!/usr/bin/env python

doHT = 4 # 0-6

import subprocess
import glob
import os

def run_cmd(command):
  print "executing command = '%s'" % command
  p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE) # , stderr = subprocess.PIPE
  stdout, stderr = p.communicate()
  return stdout

HTparts = ["HT200to300", "HT300to500", "HT500to700", "HT700to1000", "HT1000to1500", "HT1500to2000","HT2000toInf"]
part1 ='./DelphesCMSFWLite cards/gen_card.tcl /eos/user/a/acarvalh/VLQNLO/QCD_'+HTparts[doHT]+'_'
part2 = '.root '

files = glob.glob('/eos/cms/store/mc/RunIISpring18MiniAOD/QCD_'+HTparts[doHT]+'_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/*000/*.root')

print  len(files)
#print part1+str(1)+part2+files[1]

for line in range(0, len(files)): run_cmd(part1+str(line)+part2+files[line])
#    proc=subprocess.Popen([part1+str(line)+part2+files[line]], shell=True,stdout=subprocess.PIPE)
#    out=proc.stdout.read()

print "processed "+ str(len(files))
