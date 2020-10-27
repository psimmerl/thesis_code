from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

#gStyle.SetOptStat(1111)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


c0 = TCanvas("c0","c0",1200,450)
c0.cd(1)
hhs['vphimass_pos_helicity_pass_all_exclusivity'].SetLineColor(kRed)
hhs['vphimass_neg_helicity_pass_all_exclusivity'].SetLineColor(kBlue)
hhs['vphimass_neg_helicity_pass_all_exclusivity'].Draw()
hhs['vphimass_pos_helicity_pass_all_exclusivity'].Draw('same')
c0.SaveAs('hist_phi_events_helicity_comparison.png')

