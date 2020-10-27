from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#gStyle.SetOptStat(1111)

hhst = {}
hhsl = {}
fft = TFile(sys.argv[1])
for kk in fft.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhst[obj.GetName()] = obj

ffl = TFile(sys.argv[2])
for kk in ffl.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhsl[obj.GetName()] = obj


c0 = TCanvas("c0","c0",500,1250)
c0.Divide(1,4)
c0.cd(1)
hhst['vphimass_pass_all_exclusivity'].Draw()
hhst['vphimass_pass_all_excl_pass_coplancc0'].SetLineColor(kRed)
hhst['vphimass_pass_all_excl_pass_coplancc0'].Draw('same')
hhsl['vphimass_pass_all_excl_pass_coplancc0'].SetLineColor(kBlue)
hhsl['vphimass_pass_all_excl_pass_coplancc0'].Draw('same')
c0.cd(2)
hhst['vphimass_pass_all_exclusivity'].Draw()
hhst['vphimass_pass_all_excl_pass_coplancc1'].SetLineColor(kRed)
hhst['vphimass_pass_all_excl_pass_coplancc1'].Draw('same')
hhsl['vphimass_pass_all_excl_pass_coplancc1'].SetLineColor(kBlue)
hhsl['vphimass_pass_all_excl_pass_coplancc1'].Draw('same')
c0.cd(3)
hhst['vphimass_pass_all_exclusivity'].Draw()
hhst['vphimass_pass_all_excl_pass_coplancc2'].SetLineColor(kRed)
hhst['vphimass_pass_all_excl_pass_coplancc2'].Draw('same')
hhsl['vphimass_pass_all_excl_pass_coplancc2'].SetLineColor(kBlue)
hhsl['vphimass_pass_all_excl_pass_coplancc2'].Draw('same')
c0.cd(4)
hhst['vphimass_pass_all_exclusivity'].Draw()
hhst['vphimass_pass_all_exclusivity_pass_all_coplanarity'].SetLineColor(kRed)
hhst['vphimass_pass_all_exclusivity_pass_all_coplanarity'].Draw('same')
hhsl['vphimass_pass_all_exclusivity_pass_all_coplanarity'].SetLineColor(kBlue)
hhsl['vphimass_pass_all_exclusivity_pass_all_coplanarity'].Draw('same')
c0.SaveAs('hist_coplan_all_exclusivity_excl.png')


