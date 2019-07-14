#!/usr/bin/env python
# cd ../Delphes/CMSSW_7_6_1/src/ ; cmsenv ; cd -
## or the cmssw you had used to install delphes
# to run: ./analysis_gen.py mass_to_cut sample.root
import os, sys, time,math
import ROOT
#from ROOT import TLatex,TPad,TList,TH1,TH1F,TH2F,TH1D,TH2D,TFile,TTree,TCanvas,TLegend,SetOwnership,gDirectory,TObject,gStyle,gROOT,TLorentzVector,TGraph,TMultiGraph,TColor,TAttMarker,TLine,TDatime,TGaxis,TF1,THStack,TAxis,TStyle,TPaveText,TAttFill,TF2, gPad, TGaxis, TChain,TClass
import glob
from array import array
import numpy as np
pow = ROOT.TMath.Power
import bisect
from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as plt
import pandas
#from root_pandas import read_root

import root_numpy
from root_numpy import root2array, rec2array, array2root, tree2array

# Delphes headers
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
ROOT.gInterpreter.Declare('#include "DelphesClasses.h"')
# fastjet headers - if we want to re-do anything
#ROOT.gInterpreter.Declare('#include "fastjet/PseudoJet.hh"')
#ROOT.gInterpreter.Declare('#include "fastjet/ClusterSequence.hh"')
#ROOT.gInterpreter.Declare('#include "fastjet/tools/Filter.hh"')
#ROOT.gInterpreter.Declare('#include "fastjet/Selector.hh"')
#ROOT.gInterpreter.Declare('#include "fastjet/tools/Pruner.hh"')
#ROOT.gInterpreter.Declare('#include "fastjet/tools/MassDropTagger.hh"')
#ROOT.gInterpreter.Declare('using namespace fastjet;')

ROOT.gSystem.Load("libDelphes")

parser = OptionParser()
if len(sys.argv) < 1:
    print " Usage: Example1.py <input_file>"
    sys.exit(1)
inputFile = sys.argv[2]
VLQmass = float(sys.argv[1])
order = "LO"

inputpath="/eos/cms/store/user/acarvalh/VLQNLO_files/"

cx=1.0
nev=1.0
HTpartstoCX =  {
    "HT200to300":   float(cx/nev),
    "HT300to500":   float(cx/nev),
    "HT500to700":   float(cx/nev),
    "HT700to1000":  float(cx/nev),
    "HT1000to1500": float(cx/nev),
    "HT1500to2000": float(cx/nev),
    "HT2000toInf":  float(cx/nev)
    }

nevHTparts = {
    "HT200to300":   0,
    "HT300to500":   0,
    "HT500to700":   0,
    "HT700to1000":  0,
    "HT1000to1500": 0,
    "HT1500to2000": 0,
    "HT2000toInf":  0
    }

nfilesHTparts = {
    "HT200to300":   0,
    "HT300to500":   0,
    "HT500to700":   0,
    "HT700to1000":  0,
    "HT1000to1500": 0,
    "HT1500to2000": 0,
    "HT2000toInf":  0
    }

toProcess = [str(inputpath)+str(inputFile)]
if ".root" not in inputFile :
    toProcess = glob.glob(inputpath+'/QCD_*.root')

#file = open('/eos/user/a/acarvalh/VLQNLO_files/samplesList.txt',"w")
#file.write(str(glob.glob(inputpath+'/*.root')))
#file.close()

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


tree_name = "tree"

tuple = ROOT.TTree(tree_name, tree_name)
br_nFatJets = array('i', [0])
br_nBs = array('i', [0])
br_nJets = array('i', [0])
br_nBsFat = array('i', [0])
br_Weights  = array('d', [0.])
br_Tau21 = array('d', [0.])
br_PrunMass = array('d', [0.])
br_leadSD_subjet_pt = array('d', [0.])
br_FatMass = array('d', [0.])
br_WeightsFatLoop  = array('d', [0.])

tuple.Branch('nFatJets', br_nFatJets, 'nFatJets/I')
tuple.Branch('nBs', br_nBs, 'nBs/I')
tuple.Branch('nJets', br_nJets, 'nJets/I')
tuple.Branch('nBsFat', br_nBsFat, 'nBsFat/I')
tuple.Branch('Weights', br_Weights, 'Weights/D')
tuple.Branch('Tau21', br_Tau21, 'Tau21/D')
tuple.Branch('PrunMass', br_PrunMass, 'PrunMass/D')
tuple.Branch('leadSD_subjet_pt', br_leadSD_subjet_pt, 'leadSD_subjet_pt/D')
tuple.Branch('FatMass', br_FatMass, 'FatMass/D')
tuple.Branch('WeightsFatLoop', br_WeightsFatLoop, 'WeightsFatLoop/D')

def isbtagged(jets, GenB) :
    #print "calculate DR"
    see=0
    for bjets in GenB :
        if not bjets.Pt() > 0 or not jets.Pt() > 0 :
            if bjets.DeltaR(jets) < 0.3 : see = see + 1
        #else : print "The problem is here "+str(bjets.Pt())+" "+str(jets.Pt())
    if see > 0 : return see
    else : return 0

sign = lambda a: 1 if a>0 else -1 if a<0 else 0

#############################################################
# Loop over file list
#############################################################
onlyCount = False
nFatJets = 0
nJets = 0
nBs = 0
nBsFat = 0
PrunMass = -10.
FatMass = -10.
leadSD_subjet_pt = -10.
WeightsFatLoop = -10.
Weights = 1
Tau21 = -1.
for sample in toProcess :
    #for i in range(0,1) :
    #sample = "/eos/user/a/acarvalh/VLQNLO/QCD_HT2000toInf_1.root"
    print sample

    chain = ROOT.TChain("Delphes")
    #chain.Add(str(inputpath)+str(inputFile))
    try: chain.Add(sample)
    except IOError as e:
        print('Couldnt open the file (%s).' % e)
        continue
    # Create object of class ExRootTreeReader
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries() # 100 #
    print "The tree have "+str(numberOfEntries)+" events "+sample
    for dictionary in HTpartstoCX.keys() :
        if dictionary in sample :
            nevHTparts[dictionary] = nevHTparts[dictionary]+numberOfEntries
            nfilesHTparts[dictionary] = nfilesHTparts[dictionary]+1
    #############################################################
    # Loop over all events
    #############################################################
    if not onlyCount :
        # Get pointers to branches used in this analysis
        branchEvent = treeReader.UseBranch("Event")
        branchJet = treeReader.UseBranch("GenJet")
        branchFatJet = treeReader.UseBranch("GenJetAK8")
        branchParticle = treeReader.UseBranch("Particle")
        branchMET = treeReader.UseBranch("GenMissingET")
        weight = 1
        for dictionary in HTpartstoCX.keys() :
            if dictionary in sample :
                print dictionary
                weight = float(HTpartstoCX[dictionary])
        for entry in range(0, numberOfEntries): #
            # Load selected branches with data from specified event
            if entry > 10 : break
            treeReader.ReadEntry(entry)
            #print branchEvent.GetEntries()
            ## check if we have negative weights on the NLO samples, and how to use them
            #if(branchEvent.At(0).Weight < 0) :
            #    print "Weight was negative "+str(branchEvent.At(0).Weight)
            #    weight = weight*sign(branchEvent.At(0).Weight)
            #    negative+=1
            #else : print "Weight "+str(branchEvent.At(0).Weight)
            Weights = sign(branchEvent.At(0).Weight)
            print ("ScalePDF ", str(branchEvent.At(0).ScalePDF), str(branchEvent.At(0).PDF2), str(branchEvent.At(0).Weight) )
            #####################
            # Gen-level particles
            #####################
            Ws = []
            Topone = []
            GenBs = []
            GenMus = []
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
                     #print "W decay: "+str(branchParticle.At(genparticle.D1).PID)
                     #print "W decay 2: "+str(branchParticle.At(genparticle.D2).PID)
               if (IsPU == 0 and (abs(pdgCode) > 6000000) and (abs(pdgCode) < 9000000)): #and status==statusT ):
                    try : branchParticle.At(genparticle.D1).PID
                    except :
                        print "There was a Topone without daughter "+str(pdgCode)
                        continue
                    if(abs(branchParticle.At(genparticle.D1).PID) != abs(pdgCode) and len(Topone) < 3) :
                     Topone.append(genparticle) # the LHE information...
                     print "Q decay: "+str(branchParticle.At(genparticle.D1).PID)
                     print "Q decay 2: "+str(branchParticle.At(genparticle.D2).PID)
               if (IsPU == 0 and (abs(pdgCode) == 13) and abs(branchParticle.At(genparticle.M1).PID ) != 5 ): # > 6000000
                  if genparticle.PT > 10 :
                      dumb = ROOT.TLorentzVector()
                      dumb.SetPtEtaPhiM(genparticle.PT,genparticle.Eta,genparticle.Phi,genparticle.Mass)
                      GenBs.append(dumb)
                      #print "b mother: "+str(genparticle.M1)+" "+str(dumb.Pt())
                  #else : print "b-quark without pt"
                  #mother =  branchParticle.At(genparticle.M1)
                  #motherPID = mother.PID
                  #print " pdgid "+ str(pdgCode)
               if (IsPU == 0 and (abs(pdgCode) == 5) and abs(branchParticle.At(genparticle.M1).PID ) != 5 ): # > 6000000
                  if genparticle.PT > 10 :
                      dumb = ROOT.TLorentzVector()
                      dumb.SetPtEtaPhiM(genparticle.PT,genparticle.Eta,genparticle.Phi,genparticle.Mass)
                      GenMus.append(dumb)
            # taking the gen-jets
            RecoFatJets = []
            RecoBFatJets = []
            #print len(GenBs)
            for part in range(0, branchFatJet.GetEntries()): # add one more collection to the delphes card
                jet =  branchFatJet.At(part) # take the trimed jet
                if( jet.PT > 250 and abs(jet.Eta) < 2.5 and jet.Mass > 50 ) :
                   if (jet.PT > 0 and jet.Mass > 0) :
                       dumb = ROOT.TLorentzVector()
                       dumb.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
                       RecoFatJets.append(dumb)
                       RecoBFatJets.append(isbtagged(dumb, GenBs)) #
                       #print (isbtagged(dumb, GenBs), jet.BTagAlgo , jet.NSubJetsPruned)
                       #print (jet.Tau[1]/jet.Tau[0] , jet.SoftDroppedSubJet1.M())
                       #print "fatjet "+" "+str(dumb.Pt())
                       if jet.Tau[0] > 0 : Tau21 = jet.Tau[1]/jet.Tau[0]
                       #prumass = jet.SoftDroppedJet.M() # (jet.PrunedP4[1]+jet.PrunedP4[0]).M()
                       #print (jet.SoftDroppedJet.M())
                       PrunMass = jet.SoftDroppedJet.M()
                       FatMass = jet.Mass
                       leadSD_subjet_pt = jet.SoftDroppedSubJet1.Pt()
                       WeightsFatLoop = sign(branchEvent.At(0).Weight)
                       # GenJetAK8.SoftDroppedP4[5]
                       # GenJetAK8.Tau[5]
                       # GenJetAK8.SoftDroppedJet
                       # GenJetAK8.SoftDroppedSubJet1 / GenJetAK8.SoftDroppedSubJet2 # TLorentzVec
                       # GenJetAK8.NSubJetsSoftDropped
                       # GenJetAK8.Particles
                       # GenJetAK8.PTD
                       # GenJetAK8.PrunedP4[5]
                       # GenJetAK8.NSubJetsPruned
                       # GenJetAK8.NSubJetsSoftDropped
                   else : print "Fat jet without pt"
            RecoJets = []
            RecoBJets = []
            for part in range(0, branchJet.GetEntries()): # add one more collection to the delphes card
                jet =  branchJet.At(part) # take the trimed jet
                if( jet.PT > 25 ) :
                   dumb = ROOT.TLorentzVector()
                   dumb.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
                   RecoJets.append(dumb)
                   RecoBJets.append(isbtagged(dumb, GenBs))
                   ## using the DR with the genParticles to find out if there is a b-quark
            #print "size of FatJet collection " + str(len(RecoFatJets)) + " size of Jet collection " + str(len(RecoJets))
            numbb = 0
            for i in range(0, len(RecoBJets)) : numbb += RecoBJets[i];
            numfatbb = 0
            for i in range(0, len(RecoBFatJets)) : numfatbb += RecoBFatJets[i];
            nFatJets = len(RecoFatJets)
            nJets = len(RecoJets)
            nBs = numbb
            nBsFat = numfatbb
            #######################
            br_nFatJets[0] = int(nFatJets)
            br_nJets[0] = int(nJets)
            br_nBs[0] = int(nBs)
            br_nBsFat[0] = int(nBsFat)
            br_PrunMass[0] = float(PrunMass)
            br_FatMass[0] = float(FatMass)
            br_leadSD_subjet_pt[0] = float(leadSD_subjet_pt)
            br_WeightsFatLoop[0] = float(WeightsFatLoop)
            br_Weights[0] = float(Weights)
            br_Tau21[0] = float(Tau21)
            tuple.Fill()
    #if not onlyCount : print "Sample had "+str(len(filter(lambda x: x < 0, Weights)))+" negative weight events (total "+str(len(filter(lambda x: x > 0, Weights))+len(filter(lambda x: x < 0, Weights)))+")"
#if not onlyCount : print "Total had "+str(len(filter(lambda x: x < 0, Weights)))+" negative weight events (total "+str(len(filter(lambda x: x > 0, Weights))+len(filter(lambda x: x < 0, Weights)))+")"
print nevHTparts
print nfilesHTparts
#########################

out_file = ROOT.TFile("teste.root", 'RECREATE')
out_file.WriteTObject(tuple, tuple.GetName(), 'Overwrite')
out_file.Close()

#dataAll = pandas.DataFrame() # columns = ["nFatJets" : []]
#dataAll["nFatJets"] = nFatJets
#dataAll.to_root('teste.root')
#array2root(np.array(nFatJets, dtype=np.float32), 'test.root', mode='recreate')
#root2array('test.root')
"""
if not onlyCount :
    print "Plotting test histograms"
    plt.figure(figsize=(5,5))
    plt.hist(nFatJets, weights=Weights, bins=10, range=(0,5), normed=1, histtype='bar', label='nFatJets', fill=False, color= 'k', edgecolor='k', lw = 4)
    plt.hist(nJets, weights=Weights, bins=10, range=(0,5), normed=1, histtype='bar', label='Jets', fill=False, color= 'g', edgecolor='g', lw = 4)
    plt.hist(nBs, weights=Weights, bins=10, range=(0,5), normed=1, histtype='bar', label='nBJets', fill=False, color= 'y', edgecolor='y', lw = 4)
    plt.legend(loc='upper right')
    plt.title(" jet collections" )
    plt.xlabel("Njets_"+str(inputFile.replace(".root",""))+"")
    plt.ylabel("normalized")
    plt.savefig("Njets.pdf")
    plt.clf
    ############################
    plt.figure(figsize=(5,5))
    plt.hist(Tau21, weights=WeightsFatLoop, bins=20, normed=1, histtype='bar',  fill=False, color= 'k', edgecolor='k', lw = 4)
    plt.legend(loc='upper right')
    plt.title(" jet collections" )
    plt.xlabel("Tau21")
    plt.ylabel("normalized")
    plt.savefig("Tau21_"+str(inputFile.replace(".root",""))+".pdf")
    ############################
    plt.figure(figsize=(5,5))
    plt.hist(PrunMass, weights=WeightsFatLoop, bins=20, normed=1, histtype='bar',  fill=False, color= 'k', edgecolor='k', lw = 4)
    plt.legend(loc='upper right')
    plt.title(" jet collections" )
    plt.xlabel("PrunMass")
    plt.ylabel("normalized")
    plt.savefig("PrunMass_"+str(inputFile.replace(".root",""))+".pdf")
    ############################
    plt.figure(figsize=(5,5))
    plt.hist(leadSD_subjet_pt, weights=WeightsFatLoop, bins=20, normed=1, histtype='bar',  fill=False, color= 'k', edgecolor='k', lw = 4)
    plt.legend(loc='upper right')
    plt.title(" jet collections" )
    plt.xlabel("leadSD_subjet_pt")
    plt.ylabel("normalized")
    plt.savefig("leadSD_subjet_pt_"+str(inputFile.replace(".root",""))+".pdf")
    ############################
    plt.figure(figsize=(5,5))
    plt.hist(FatMass, weights=WeightsFatLoop, bins=20, normed=1, histtype='bar',  fill=False, color= 'k', edgecolor='k', lw = 4)
    plt.legend(loc='upper right')
    plt.title(" jet collections" )
    plt.xlabel("AK8 Jet Mass")
    plt.ylabel("normalized")
    plt.savefig("FatMass_"+str(inputFile.replace(".root",""))+".pdf")
    #########################
    # output the efficiencies # see DiHiggs project for template
    #########################
"""
