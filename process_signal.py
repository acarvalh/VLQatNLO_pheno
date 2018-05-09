#!/usr/bin/env python

import subprocess
import glob
import os

def run_cmd(command):
  print "executing command = '%s'" % command
  p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE) # , stderr = subprocess.PIPE
  stdout, stderr = p.communicate()
  return stdout

part1 ='./DelphesHepMC cards/gen_card.tcl '

files = procP1=glob.glob('/eos/user/a/acarvalh/VLQNLO_files/T_W_W*.hepmc.gz')

print  len(files)
#print part1+str(1)+part2+files[1]

for line in range(2, len(files)):
    print files[line]
    run_cmd('gunzip '+files[line])
    hepmcFile = files[line].replace('.hepmc.gz', '.hepmc')
    outputFile = files[line].replace('.hepmc.gz', '.root')
    #print part1+outputFile+' '+hepmcFile
    outputlog = run_cmd(part1+outputFile+' '+hepmcFile)
    file = open(files[line].replace('.hepmc.gz', '.log'),"w")
    file.write(outputlog)
    file.close()
    run_cmd('gzip '+hepmcFile)

print "processed "+ str(len(files))
