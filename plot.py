#!/usr/bin/env python

import os, sys, ROOT

schemes = {
    "4F_LO":  {"fname": "T_L_Qjq_4FNS_LO_muonsWZHdecay_M1200GeV/output/out_T_W_W_4FNS_Qjq_L_muonsWZHdecay_M1200GeV_Kprod0.1_Kdecay0.1_LHCEnergy13TeV_pythia8_events.root", "col": 9}, 
    "4F_NLO": {"fname": "T_L_Qjq_4FNS_NLO_muonsWZHdecay_M1200GeV/output/out_events_PYTHIA8_0.root", "col": 29},
    "5F_LO":  {"fname": "T_L_Qj_5FNS_LO_muonsWZHdecay_M1200GeV/output/out_T_W_W_5FNS_Qj_L_muonsWZHdecay_M1200GeV_Kprod0.1_Kdecay0.1_LHCEnergy13TeV_pythia8_events.root", "col": 8},
    "5F_NLO": {"fname": "T_L_Qj_5FNS_NLO_muonsWZHdecay_M1200GeV/output/out_events_PYTHIA8_0.root", "col": 28},
    }

vars = [
    "T_mass",
    "HT",
    "ST",
    ]

for v in vars:

  ROOT.gROOT.SetBatch(1)

  c = ROOT.TCanvas("c_{}".format(v), "", 800, 600)
  c.cd()

  hists = {}

  ymax = 0

  for scheme in schemes:

    print scheme

    f = ROOT.TFile.Open(schemes[scheme]["fname"], "READ")
    col = schemes[scheme]["col"]

    h         = f.Get("h_"+v)
    h_scaleHi = f.Get("h_"+v+"_scaleHi")
    h_scaleLo = f.Get("h_"+v+"_scaleLo")

    g = ROOT.TGraphAsymmErrors(h.GetNbinsX())
    g.SetName("g_{0}_{1}".format(v, scheme))
    g.SetLineColor(col)
    g.SetFillColor(col)
    g.SetFillStyle(3004)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(col)
    for i in range(1, h.GetNbinsX()+1):
      g.SetPoint(i-1, h.GetBinCenter(i), h.GetBinContent(i))
      g.SetPointEXlow(i-1, h.GetBinWidth(i)/2)
      g.SetPointEXhigh(i-1, h.GetBinWidth(i)/2)

      y_scaleHi = h_scaleHi.GetBinContent(i) - h.GetBinContent(i)
      y_scaleLo = h_scaleLo.GetBinContent(i) - h.GetBinContent(i)

      if y_scaleHi > 0:
        eyHi = abs(y_scaleHi)
        eyLo = abs(y_scaleLo)
      else:
        eyHi = abs(y_scaleLo)
        eyLo = abs(y_scaleHi)

      g.SetPointEYlow (i-1, eyLo)
      g.SetPointEYhigh(i-1, eyHi)

      if max(h_scaleHi.GetMaximum(),  h_scaleLo.GetMaximum(), h.GetMaximum()) > ymax: 
        ymax = max(h_scaleHi.GetMaximum(),  h_scaleLo.GetMaximum(), h.GetMaximum())

    hists[scheme] = g

  print hists

  xmin = h.GetBinLowEdge(1)
  xmax = h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())
  ymin = 0

  print xmin, xmax, ymin, ymax
      
  h0 = c.DrawFrame(xmin, ymin, xmax, ymax*1.2, "{}".format(h.GetTitle()))

  leg = ROOT.TLegend(0.6,0.6,0.88,0.88,"","brNDC")
  leg.SetBorderSize(0)

  for scheme in schemes:
    hists[scheme].Draw("E2")
    leg.AddEntry(hists[scheme], scheme, "FP")

  leg.Draw()

  c.SaveAs(c.GetName()+".pdf")
