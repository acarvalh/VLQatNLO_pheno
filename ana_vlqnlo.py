#!/usr/bin/env python
from __future__ import print_function

DEBUG = False

bashscripttemplate='''#!/bin/bash

cd RUNDIR
source /afs/cern.ch/work/d/dm2/Analysis/VLQNLOPheno/VLQatNLO_pheno/setup.bash
cd RUNDIR

python /afs/cern.ch/work/d/dm2/Analysis/VLQNLOPheno/VLQatNLO_pheno/ana_vlqnlo.py -f FNAME

#find -type l -delete

'''

condor_template = """universe              = vanilla
executable            = EXEC
arguments             = $(ClusterID) $(ProcId)
output                = OUTPUT/job_JOB_NUMBER.$(ClusterId).$(ProcId).out
error                 = OUTPUT/job_JOB_NUMBER.$(ClusterId).$(ProcId).err
log                   = OUTPUT/job_JOB_NUMBER.$(ClusterId).log
+JobFlavour           = "QUEUE"
queue
"""

def Pt(jet): return jet.PT




def Eta(jet): return abs(jet.Eta)




def isBTagged(jet):
  btageff_b = 0.7
  btageff_c = 0.2
  btageff_l = 0.01

  from random import random

  if jet.Flavor == 5:
    if random() > btageff_b: return True
    else: return False
  elif jet.Flavor == 4:
    if random() > btageff_c: return True
    else: return False
  else:
    if random() > btageff_l: return True
    else: return False




def analyze_gen(fname, maxEvents=-1):

  import os

  if fname[0] == "#": return

  import ROOT
  import numpy as np
  ROOT.gROOT.SetBatch()
  print('Analyzing {}'.format(fname))

  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
  ROOT.gInterpreter.Declare('#include "DelphesClasses.h"')

  ROOT.gSystem.Load("libDelphes")

  chain = ROOT.TChain("Delphes")
  try: chain.Add(fname)

  except IOError as e:
    print('Couldnt open the file (%s).' % e)
    return -1

  treeReader = ROOT.ExRootTreeReader(chain)
  numberOfEntries = treeReader.GetEntries()

  if DEBUG: print("The tree have {0} from file {1}".format(numberOfEntries, fname))

  branchEvent    = treeReader.UseBranch("Event")
  branchPart     = treeReader.UseBranch("Particle")

  if maxEvents == -1:
    toprocess = numberOfEntries
  else:
    toprocess = maxEvents

  wtfile = open('/eos/cms/store/user/acarvalh/VLQNLO_files/withPDFunc/Results_20181106/Results_T_L_Qjq_4FNS_NLO_muonsWZHdecay/Events/ASCII/T_W_W_4FNS_Qjq_L_muonsWZHdecay_M1200GeV_Kprod0.1_Kdecay0.1_LHCEnergy13TeV_systematics.dat', 'read')
  wtlines = wtfile.readlines()

  for event in range(0, toprocess): #
    treeReader.ReadEntry(event)
    for part in branchPart:
      print('part id = {}'.format(part.PID))
      wts = wtlines[event+1]
      #print(len(wts))




def analyze(fname, maxEvents=-1):

  import os

  if fname[0] == "#": return

  import ROOT
  import numpy as np
  ROOT.gROOT.SetBatch()
  print('Analyzing {}'.format(fname))

  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
  ROOT.gInterpreter.Declare('#include "DelphesClasses.h"')

  ROOT.gSystem.Load("libDelphes")

  chain = ROOT.TChain("Delphes")
  try: chain.Add(fname)
  except IOError as e:
    print('Couldnt open the file (%s).' % e)
    return -1

  treeReader = ROOT.ExRootTreeReader(chain)
  numberOfEntries = treeReader.GetEntries()
  if DEBUG: print("The tree have {0} from file {1}".format(numberOfEntries, fname))

  branchEvent    = treeReader.UseBranch("Event")
  branchJet      = treeReader.UseBranch("GenJet")
  branchFatJet   = treeReader.UseBranch("GenJetAK8")
  branchParticle = treeReader.UseBranch("Particle")
  branchMET      = treeReader.UseBranch("GenMissingET")

  fout = ROOT.TFile("out_"+os.path.basename(fname), "recreate")
  fout.cd()

  h_cutflow = ROOT.TH1D("h_cutflow", "", 5, -0.5, 4.5)
  h_cutflow.GetXaxis().SetBinLabel(1, "All")
  h_cutflow.GetXaxis().SetBinLabel(2, " N(W jet)>0")
  h_cutflow.GetXaxis().SetBinLabel(3, "N(b jet)>0")
  h_cutflow.GetXaxis().SetBinLabel(4, "p_{T}(1st b jet)>200")
  h_cutflow.GetXaxis().SetBinLabel(5, "Fwd jet")

  h_HT              = ROOT.TH1D("h_HT"            , ";H_{T};Events/100 GeV;",20, 1000, 3000)
  h_T_mass          = ROOT.TH1D("h_T_mass"        , ";M(T);Events/100 GeV;",20, 1000, 3000)
  h_wjets_pt        = ROOT.TH1D("h_wjets_pt"      , ";p_{T}(W jets);Events/50 GeV;",20, 200., 1200.)
  h_wjet1_pt        = ROOT.TH1D("h_wjet1_pt"      , ";p_{T}(leading-p_{T} W jet);Events/50 GeV;",20, 200., 1200.)
  h_bjets_pt        = ROOT.TH1D("h_bjets_pt"      , ";p_{T}(b jets);Events/50 GeV;",24, 0., 1200.)
  h_bjet1_pt        = ROOT.TH1D("h_bjet1_pt"      , ";p_{T}(leading-p_{T} b jet);Events/50 GeV;",24, 0., 1200.)
  h_bjet1_fl        = ROOT.TH1D("h_bjet1_fl"      , ";Leading-p_{T} b jet flavour;Events/1 unit;",3,-0.5,5.5)
  h_nwjets          = ROOT.TH1D("h_nwjets"        , ";Number of W jets;Events/ 1 unit;", 6, -0.5, 5.5)
  h_nbjets          = ROOT.TH1D("h_nbjets"        , ";Number of b jets;Events/ 1 unit;", 6, -0.5, 5.5)
  h_forwardjet_pt   = ROOT.TH1D("h_forwardjet_pt" , ";p_{T}(forward jet);Events/10 GeV;",20, 0., 200.)
  h_forwardjet_eta  = ROOT.TH1D("h_forwardjet_eta", ";#eta(forward jet);Events/ 0.2 units;",50, -5, 5)
  h_mu_pt   = ROOT.TH1D("h_mu_pt" , ";p_{T}(mu);Events/10 GeV;",20, 0., 200.)
  h_mu_eta  = ROOT.TH1D("h_mu_eta", ";#eta(mu);Events/ 0.2 units;",50, -5, 5)

  if maxEvents == -1:
    toprocess = numberOfEntries
  else:
    toprocess = maxEvents

  wtfile = open(fname.replace('root', 'dat').replace('Events', 'Events/ASCII')\
      .replace('_pythia8_events', '_systematics'), \
      'read')
  wtlines = wtfile.readlines()

  for event in range(0, toprocess):
    #print ("=================")
    if event > 100 : break
    treeReader.ReadEntry(event)

    wts = wtlines[event+1]
    wts = wts.split()
    #print(len(wts))
    wt_scaleDouble = wts[4]
    wt_scaleHalf = wts[8]

    h_cutflow.Fill(0)

    ### muons
    muons = []
    #electrons = []
    for part in range(0, branchParticle.GetEntries()):
       genparticle =  branchParticle.At(part)
       #status = genparticle.M1
       #########################
       if genparticle.IsPU == 0 and (abs(genparticle.PID) == 13) and genparticle.Status == 1 : #
           # print (event, genparticle.PID, genparticle.Status, genparticle.M1)
           if ( genparticle.PT > 25 and abs(genparticle.Eta) < 2.5 ) : muons.append(genparticle)
       #if genparticle.IsPU == 0 and (abs(genparticle.PID) == 11) and genparticle.Status == 23  :
       #     print (event, pdgCode, genparticle.Status)
       #     if ( genparticle.PT > 25 and abs(genparticle.Eta) < 2.5 ) : electrons.append(genparticle)

    ### W jets
    wjets = []
    for jet in branchFatJet:
      if (jet.PT > 400. and \
          abs(jet.Eta) < 2.5 \
          and jet.Tau[0] != 0 \
          and jet.Tau[1]/jet.Tau[0] < 0.6 \
          and jet.SoftDroppedJet.M() > 65 \
          and jet.SoftDroppedJet.M() < 105 \
          ):
        wjets.append(jet)
        h_wjets_pt.Fill(jet.PT)

    bjets = []
    for jet in branchJet:
      if (jet.PT > 50. \
          and abs(jet.Eta) < 2.5 \
          ):
        p4_jet = jet.P4()
        for wjet in wjets:
          p4_wjet = wjet.P4()
          if p4_jet.DeltaR(p4_wjet) < 1.2:
            continue
        for muon in muons :
            p4_muon = muon.P4()
            if p4_jet.DeltaR(p4_muon) < 0.4:
              continue
        jet.Flavor = 0
        for genparticle in branchParticle:
          ### Match to B
          if genparticle.PID//100%10 == 5 or \
              genparticle.PID//1000%10 == 5 :
                p4_b = genparticle.P4()
                if p4_jet.DeltaR(p4_b) < 0.4:
                  jet.Flavor = 5
                  break
          ### Match to D
          if genparticle.PID//100%10 == 4 or \
              genparticle.PID//1000%10 == 4 :
                p4_c = genparticle.P4()
                if p4_jet.DeltaR(p4_c) < 0.4:
                  jet.Flavor = 4
                  break
        if DEBUG: print("jet flavor = {}".format(jet.Flavor))
        ### Do b-tagging
        if isBTagged(jet):
          bjets.append(jet)
          h_bjets_pt.Fill(jet.PT)

    forwardjets = []
    HT = 0
    for jet in branchJet:
      if (jet.PT > 20 \
          and abs(jet.Eta) > 1.5 \
          and abs(jet.Eta) < 5 \
          ):
        p4_jet = jet.P4()
        HT += p4_jet.Pt()
        for wjet in wjets:
          p4_wjet = wjet.P4()
          if p4_jet.DeltaR(p4_wjet) < 1.2:
            continue
        for bjet in bjets:
          p4_bjet = bjet.P4()
          if p4_jet.DeltaR(p4_bjet) < 0.8:
            continue
          forwardjets.append(jet)

    wjets.sort(key=Pt)
    bjets.sort(key=Pt)
    forwardjets.sort(key=Eta)

    h_nwjets.Fill(len(wjets))

    if len(wjets) < 1: continue
    h_cutflow.Fill(1)

    h_wjet1_pt.Fill(wjets[0].PT)

    h_nbjets.Fill(len(bjets))

    if len(bjets) < 2: continue
    h_cutflow.Fill(2)

    h_bjet1_pt.Fill(bjets[0].PT)
    h_bjet1_fl.Fill(bjets[0].Flavor)

    if bjets[0].PT < 200.: continue
    h_cutflow.Fill(3)

    if len(forwardjets) > 0:
      h_forwardjet_pt .Fill(forwardjets[0].PT)
      h_forwardjet_eta.Fill(forwardjets[0].Eta)
    else:
      h_forwardjet_pt .Fill(0)
      h_forwardjet_eta.Fill(-100)

    if len(forwardjets) < 1: continue
    if abs(forwardjets[0].Eta) < 3.: continue
    h_cutflow.Fill(4)

    ### Reconstruct T->Wb
    h_T_mass.Fill( (wjets[0].P4() + bjets[0].P4()).Mag() )
    h_HT.Fill(HT)

    if len(muons) :
        h_mu_pt.Fill(muons[0].PT)
        h_mu_eta.Fill(muons[0].Eta)
    print (event, len(muons), len(forwardjets), len(bjets))

  fout.Write()
  fout.Close()




def submit(sample, files,  maxEvents, queue):
  import os, re

  print(files)

  workdir = sample
  if not re.search("^/", workdir):
    workdir = os.path.join(os.getcwd(),workdir)

  if not os.path.exists(workdir):
    os.mkdir(sample)
    os.mkdir(os.path.join(workdir,'input'))
    os.mkdir(os.path.join(workdir,'output'))

  rundir = os.path.join(workdir, 'output')

  for f in files:
    bashscript = open(os.path.join(workdir,'input', os.path.basename(f.replace('.root', '.sh'))), 'w')
    bashscriptcontent = re.sub('RUNDIR', rundir, bashscripttemplate)
    bashscriptcontent = re.sub('FNAME', f, bashscriptcontent)
    bashscript.write(bashscriptcontent)
    bashscript.close()

    condor_script = open(os.path.join(workdir,'input', os.path.basename(f.replace('.root', '.condor'))), 'w')
    condor_script_content = re.sub('EXEC', os.path.join(workdir,'input', os.path.basename(f.replace('.root', '.sh'))), condor_template)
    try:
      jobnum = f.rstrip('root').split('_')[-1]
    except:
      jobnum = ''
    condor_script_content = re.sub('JOB_NUMBER', jobnum, condor_script_content)
    condor_script_content = re.sub('OUTPUT',os.path.join(workdir,'output'),condor_script_content)
    condor_script_content = re.sub('QUEUE',queue,condor_script_content)
    condor_script.write(condor_script_content)
    condor_script.close()

    cmd = ' '.join(['condor_submit',  os.path.join(workdir,'input', os.path.basename(f.replace('.root', '.condor')))])
    print(cmd)
    os.system(cmd)




def main():

  from argparse import ArgumentParser
  parser = ArgumentParser(description="Do -h to see usage")

  parser.add_argument("-g", "--analyze_gen",
      dest="analyze_gen",
      action="store_true",
      help="Process only GEN info")

  parser.add_argument("-b", "--batch",
      dest="batch",
      action="store_true",
      help="Process files in batch mode")

  parser.add_argument("-q", "--queue",
      dest="queue",
      action="store",
      default="longlunch",
      help="Condor job flavours")

  parser.add_argument("-f", "--file",
      dest="infile",
      action="store",
      default="filestoprocess.txt",
      type=str,
      help="Name of text file containing names of ROOT files to process")

  parser.add_argument("-n", "--maxEvents",
      dest="maxEvents",
      action="store",
      default=-1,
      type=int,
      help="Maximum number of events to process")

  args = parser.parse_args()

  print(args)

  if args.batch == True:
    print("Setting to batch mode")
    import json
    import functools
    with open(args.infile, 'r') as infile:
      samples = json.load(infile)
      for sample in samples:
        print(sample)
        submit(sample, samples[sample]["files"], args.maxEvents, args.queue)
  else:
    if args.analyze_gen == True:
      analyze_gen(args.infile, args.maxEvents)
    else:
      analyze(args.infile, args.maxEvents)




if __name__ == "__main__":
  main()
