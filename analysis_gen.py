#!/usr/bin/env python
# to run: ./analysis.py delphes_output.root
import os, sys, time,math
import ROOT
#from ROOT import TLatex,TPad,TList,TH1,TH1F,TH2F,TH1D,TH2D,TFile,TTree,TCanvas,TLegend,SetOwnership,gDirectory,TObject,gStyle,gROOT,TLorentzVector,TGraph,TMultiGraph,TColor,TAttMarker,TLine,TDatime,TGaxis,TF1,THStack,TAxis,TStyle,TPaveText,TAttFill,TF2, gPad, TGaxis, TChain,TClass
from array import array
import numpy as np
pow = ROOT.TMath.Power
import bisect
from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as plt
# Delphes headers
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
ROOT.gInterpreter.Declare('#include "DelphesClasses.h"')
# fastjet headers - to make prunned mass
#ROOT.gInterpreter.Declare('#include "external/fastjet/PseudoJet.hh"')
#ROOT.gInterpreter.Declare('#include "external/fastjet/ClusterSequence.hh"')
#ROOT.gInterpreter.Declare('#include "external/fastjet/tools/Filter.hh"')
#ROOT.gInterpreter.Declare('#include "external/fastjet/Selector.hh"')
#ROOT.gInterpreter.Declare('#include "external/fastjet/tools/Pruner.hh"')
#ROOT.gInterpreter.Declare('#include "external/fastjet/tools/MassDropTagger.hh"')
#ROOT.gInterpreter.Declare('using namespace fastjet;')
#http://spartyjet.hepforge.org/
#from spartyjet import *

ROOT.gSystem.Load("libDelphes")

parser = OptionParser()
if len(sys.argv) < 1:
    print " Usage: Example1.py <input_file>"
    sys.exit(1)
inputFile = sys.argv[1]
VLQmass = float(sys.argv[2])
order = "LO"

inputpath="/eos/user/a/acarvalh/VLQNLO/"
#########################
# Cuts
#########################
Tau21cut = 1
Tau31cut = 1
PrunMass2 = 100
PrunMass3 = 1000
VLQresolution = 200
Hresolution=50
etab = 2.4
bjetpt = 30
#########################
# categorization
#########################
nFatJets = []
nBs = []
nJets = []
nBsFat = []

def isbtagged(jets, GenBJet) :
    see=0
    #for j in range(0, jets.Particles.GetEntriesFast()) :
    #    #print "constituents flavour " + str(jets.Particles[j].PID)
    #    if(abs(jets.Particles[j].PID) == 5) :
    #        if(jets.Particles[j].PT > bjetpt and jets.Particles[j].Eta < etab) :
    #            see = see + 1
    for bjets in GenBJet :
        if bjets.DeltaR(jets) < 0.3 : see = see + 1
    if see > 0 : return see
    else : return 0

#########################
# Create chain of root trees
#s=ROOT.MyStruct_0_0()
chain = ROOT.TChain("Delphes")
chain.Add(str(inputpath)+str(inputFile))
# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries() # 100 #
print "The tree have "+str(numberOfEntries)+" events"
# Get pointers to branches used in this analysis
branchEvent = treeReader.UseBranch("Event")
branchJet = treeReader.UseBranch("GenJet")
branchFatJet = treeReader.UseBranch("GenJetAK8")
branchParticle = treeReader.UseBranch("Particle")
#############################################################
# Declare histograms
#############################################################

#############################################################
# Loop over all events
#############################################################
for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)
    weight = 1 #branchEvent.At(0).Weight
    #print weight
    #print "entry : "+str(entry)
    if(weight < 0 ) :
        negative+=1
        continue
    #####################
    # Gen-level particles
    #####################
    Ws = []
    Topone = []
    GenBJets = []
    QQ = True
    #print branchParticle.GetEntries()
    for part in range(0, branchParticle.GetEntries()):
       genparticle =  branchParticle.At(part)
       pdgCode = genparticle.PID
       #print pdgCode
       IsPU = genparticle.IsPU
       status = genparticle.M2 # genparticle.Status
       # check if it is the correct status (for QQ the last 25 is 52 and the last topone 62)
       #print " pdgid "+ str(pdgCode)+" status "+str(status)
       if(IsPU == 0 and (pdgCode == 24)):
          mother =  branchParticle.At(genparticle.M1)
          motherPID = mother.PID
          #print "H mother: "+str(motherPID)
          if(branchParticle.At(genparticle.D1).PID != 24) :
             Ws.append(genparticle) # find other way to follow
             print "W decay: "+str(branchParticle.At(genparticle.D1).PID)
             print "W decay 2: "+str(branchParticle.At(genparticle.D2).PID)
       if (IsPU == 0 and (abs(pdgCode) > 6000000)): #and status==statusT ):
          if(abs(branchParticle.At(genparticle.D1).PID) != abs(pdgCode) and len(Topone) <3) :
             Topone.append(genparticle) # the LHE information...
             print "Q decay: "+str(branchParticle.At(genparticle.D1).PID)
             print "Q decay 2: "+str(branchParticle.At(genparticle.D2).PID)
       if (IsPU == 0 and (abs(pdgCode) == 5) and abs(branchParticle.At(genparticle.M1).PID )> 6000000):
          dumb = ROOT.TLorentzVector()
          dumb.SetPtEtaPhiM(genparticle.PT,genparticle.Eta,genparticle.Phi,genparticle.Mass)
          GenBJets.append(dumb)
          #print "b mother: "+str(genparticle.M1)
          #mother =  branchParticle.At(genparticle.M1)
          #motherPID = mother.PID
          print " pdgid "+ str(pdgCode)
    # taking the gen-jets
    RecoFatJets = []
    RecoBFatJets = []
    for part in range(0, branchFatJet.GetEntries()): # add one more collection to the delphes card
        jet =  branchFatJet.At(part) # take the trimed jet
        if( jet.PT > 250 and jet.PT < 1500 and abs(jet.Eta) < 2 and jet.Mass > 50 ) :
           dumb = ROOT.TLorentzVector()
           dumb.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
           RecoFatJets.append(dumb)
           RecoBFatJets.append(isbtagged(dumb, GenBJets))
    RecoJets = []
    RecoBJets = []
    for part in range(0, branchJet.GetEntries()): # add one more collection to the delphes card
        jet =  branchJet.At(part) # take the trimed jet
        if( jet.PT > 250 and jet.PT < 1500 and abs(jet.Eta) < 2 and jet.Mass > 50 ) :
           dumb = ROOT.TLorentzVector()
           dumb.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
           RecoJets.append(dumb)
           ## using the constituents to find out if there is a b-quark
           RecoBJets.append(isbtagged(dumb, GenBJets))

    print "size of FatJet collection " + str(len(RecoFatJets)) + " size of Jet collection " + str(len(RecoJets))
    print RecoBJets
    numbb = 0
    for i in range(0, len(RecoBJets)) : numbb += RecoBJets[i];
    numfatbb = 0
    for i in range(0, len(RecoBFatJets)) : numfatbb += RecoBFatJets[i];
    print " total "+str(numbb)
    nFatJets.append(len(RecoFatJets))
    nJets.append(len(RecoJets))
    nBs.append(numbb)
    nBsFat.append(numbb)

#########################
print "Plotting test histograms"
plt.figure(figsize=(5,5))
plt.hist(nFatJets, bins=10, range=(0,5), normed=1, histtype='bar', label='nFatJets', fill=False, color= 'k', edgecolor='k', lw = 4)
plt.hist(nJets, bins=10, range=(0,5), normed=1, histtype='bar', label='Jets', fill=False, color= 'g', edgecolor='g', lw = 4)
plt.hist(nBs, bins=10, range=(0,5), normed=1, histtype='bar', label='nBJets', fill=False, color= 'y', edgecolor='y', lw = 4)
plt.legend(loc='upper right')
plt.title(" jet collections" )
plt.xlabel("Njets")
plt.ylabel("normalized")
plt.savefig("Njets.pdf")
#########################
# output the efficiencies # see DiHiggs project for template
#########################
