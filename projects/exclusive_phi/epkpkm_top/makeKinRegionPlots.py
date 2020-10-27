from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

gStyle.SetOptStat(0000)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


c0 = TCanvas("c0","c0",1200,600)
c0.Divide(2,1)
c0.cd(1)
hhs['xq2_pass_all_exclusivity'].SetTitle('Q2 vs xb; xb; Q2 (GeV^{2})')
hhs['xq2_pass_all_exclusivity'].Draw('colz')
c0.cd(2)
hhs['txb_pass_all_exclusivity'].SetTitle('-t vs xb; -t (GeV^{2}); xb')
hhs['txb_pass_all_exclusivity'].Draw('colz')
c0.SaveAs('q2x_all_excl.png')

c1 = TCanvas("c1","c1",1000,1000)
c1.cd(1)
hhs['txb_pass_all_exclusivity'].SetTitle('-t vs xb; -t (GeV^{2}); xb')
hhs['txb_pass_all_exclusivity'].Draw('colz')
c1.SaveAs('tx_all_excl.png')



