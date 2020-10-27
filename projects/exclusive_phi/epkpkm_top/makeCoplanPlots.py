from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

gStyle.SetOptStat(000000)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


hpro_a=hhs['cplpro_pass_all_exclusivity']
hpro_b=hhs['cplpro_raw']
hkp_a=hhs['cplkp_pass_all_exclusivity']
hkp_b=hhs['cplkp_raw']
hkm_a=hhs['cplkm_pass_all_exclusivity']
hkm_b=hhs['cplkm_raw']

c0 = TCanvas("c0","c0",1200,450)
c0.Divide(3,1)
c0.cd(1)
hpro_b.SetTitle('Angle Between Measured and Calculated Proton with eK^{+}K^{-}X; #theta (deg); counts')
hpro_b.SetLineColor(kBlack)
hpro_a.SetLineColor(kRed)
hpro_b.SetAxisRange(0.0,20.0,'X')
hpro_b.Draw()
hpro_a.Draw('same')
hpro_b.SetMinimum(0)
c0.cd(2)
hkp_b.SetTitle('Angle Between Measured and Calculated K^{+} with epK^{-}X; #theta (deg); counts')
hkp_b.SetLineColor(kBlack)
hkp_a.SetLineColor(kRed)
hkp_b.SetAxisRange(0.0,20.0,'X')
hkp_b.Draw()
hkp_a.Draw('same')
c0.cd(3)
hkm_b.SetTitle('Angle Between Measured and Calculated K^{-} with epK^{+}X; #theta (deg); counts')
hkm_b.SetLineColor(kBlack)
hkm_a.SetLineColor(kRed)
hkm_b.SetAxisRange(0.0,20.0,'X')
hkm_b.Draw()
hkm_a.Draw('same')
c0.SaveAs('hist_coplan_all_exclusivity.png')

