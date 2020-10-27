from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

gStyle.SetOptStat(00000)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj



c0 = TCanvas("c0","c0",1200,450)
lraw=hhs['lambdamass_raw']
lfin=hhs['lambdamass_pass_all_exclusivity']
#plot raw and final lamda spectrum on top of each other
lraw.SetLineColor(kBlack)
lfin.SetLineColor(kRed)
lraw.Draw()
lfin.Draw('same')
c0.SaveAs('hist_lambdamass_before_after_inbending_epkpkm.png')


#c0.Update()
#p1=c3b.GetPad(1)
#ymax=gPad.GetUymax()


c1 = TCanvas('c1','c1',1800,900)
c1.Divide(2,1)
c1.cd(1)
hhs['lambdaphimass_raw'].SetTitle('Invariant Mass of K^{+}K^{-} vs Invariant Mass of PrK^{-}; Invariant Mass of PrK^{-} (GeV); Invariant Mass of  K^{+}K^{-} (GeV)')
hhs['lambdaphimass_raw'].Draw('colz')
c1.cd(2)
hhs['lambdaphimass_pass_all_exclusivity'].SetTitle('Invariant Mass of K^{+}K^{-} vs Invariant Mass of PrK^{-}; Invariant Mass of PrK^{-} (GeV); Invariant Mass of  K^{+}K^{-} (GeV)')
hhs['lambdaphimass_pass_all_exclusivity'].Draw('colz')
#lmin = TLine(1.49, 0.8, 1.49, 1.5)
#lmax = TLine(1.587,0.8, 1.587, 1.5)
#lmin.SetLineWidth(3)
#lmin.SetLineColor(kRed)
#lmin.Draw('same')
#lmax.SetLineWidth(3)
#lmax.SetLineColor(kRed)
#lmax.Draw('same')
c1.SaveAs('hist_lambdaphimass_before_after_inbending_epkpkm.png')

c2 = TCanvas('c2','c2',1000,1000)
c2.Divide(1,1)
c2.cd(1)
#hhs['lambdaphimass_pass_all_exclusivity'].Rebin2D(2,2)
hhs['lambdaphimass_pass_all_exclusivity'].Draw('colz')
c2.SaveAs('hist_lambdaphimass_pass_all_exclusivity_inbending_epkpkm.png')

