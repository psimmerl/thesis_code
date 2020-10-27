from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile
from ROOT import TGraphErrors
from ROOT import TLorentzVector
from ROOT import gROOT, gBenchmark, gStyle, gPad, kRed
from collections import OrderedDict
from array import array
import numpy as np
import math
import sys, os


gStyle.SetOptStat(00000)
gStyle.SetPalette(55);

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj

pr_b_thry = TF1("pr_b_thry","x/sqrt(x*x + 0.938*0.938)",0.10,8.0)
kp_b_thry = TF1("kp_b_thry","x/sqrt(x*x + 0.493*0.493)",0.10,8.0)
km_b_thry = TF1("km_b_thry","x/sqrt(x*x + 0.497*0.493)",0.10,8.0)

c1 = TCanvas("c1","c1",1200,450)
c1.Divide(3,1)
c1.cd(1)
hhs['pro_beta_vs_p_pass_all_exclusivity'].SetTitle("Final Proton Beta vs P; p (GeV); #beta")
hhs['pro_beta_vs_p_pass_all_exclusivity'].Draw('colz')
pr_b_thry.SetLineColor(kRed)
pr_b_thry.SetLineStyle(2)
pr_b_thry.Draw("same")

c1.cd(2)
hhs['kp_beta_vs_p_pass_all_exclusivity'].SetTitle("Final K^{+} Beta vs P; p (GeV); #beta")
hhs['kp_beta_vs_p_pass_all_exclusivity'].Draw('colz')
kp_b_thry.SetLineColor(kRed)
kp_b_thry.SetLineStyle(2)
kp_b_thry.Draw("same")

c1.cd(3)
hhs['km_beta_vs_p_pass_all_exclusivity'].SetTitle("Final K^{-} Beta vs P; p (GeV); #beta")
hhs['km_beta_vs_p_pass_all_exclusivity'].Draw('colz')
km_b_thry.SetLineColor(kRed)
km_b_thry.SetLineStyle(2)
km_b_thry.Draw("same")
c1.SaveAs("beta_analysis_all_hadrons_all_comb.png")
