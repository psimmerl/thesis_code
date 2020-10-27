from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TMath, TFormula
from ROOT import TLegend
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad, kBlack, kRed, kBlue
from array import array
import math
import sys


gStyle.SetOptStat(0000)


fin = TFile(sys.argv[1])

phim = fin.Get('vphimass_pass_all_exclusivity')

gaus = "[0]*TMath::Gaus(x,[1],[2])"
hvyside = "(x>[4])*2/TMath::Pi()*TMath::ATan( (x-[4])/[5])"
bck = "[3]*TMath::Exp(-TMath::Power(TMath::Abs([6]*(x-[4])),[7]))"
form = gaus +"+"+ hvyside + "*" + bck

#print(" my function to fit phi is " + form)
phif = TF1("fitPhi",form, 0.98, 1.3)
print(phif)
fitparm = [300.2, 1.02, 0.022, 50, 0.987, 0.0045, 2.19, 3]
for pp in range(len(fitparm)):
    print("par number %d, value of %f" %(pp, fitparm[pp]) )
    phif.SetParameter(pp,fitparm[pp])
    phif.SetParLimits(pp, 0.5*fitparm[pp], 1.5*fitparm[pp])

phif.SetNpx(10000)    
phif.FixParameter(4,fitparm[4])
phim.Fit("fitPhi","BRN")

bck = TF1("bck",form,0.98,1.3)
bck.SetNpx(10000)
for pp in range(len(fitparm)):
    print("par number %d, value of %f" %(pp, fitparm[pp]) )
    bck.SetParameter(pp,phif.GetParameter(pp))

bck.FixParameter(0,0)


c1 = TCanvas("c1","c1",900,500)
c1.cd(1)
phim.SetTitle("Invariant Mass of Charged Kaons from epK^{+}K^{-};IM K^{+}K^{-} (GeV);counts")
phim.GetXaxis().CenterTitle()
phim.GetYaxis().CenterTitle()
phim.Draw("PE")
bck.SetLineColor(4)
phif.Draw("same")
bck.Draw("same")
c1.SaveAs("phiPeak.png")


phiraw = fin.Get('vphimass_raw')
phipall= fin.Get('vphimass_pass_all_exclusivity')
phipalllmbd = fin.Get('vphimass_pass_all_exclusivity_pass_lmdacut')


############################################
############################################
#compare phi mass 
c2 = TCanvas('c2','c2',900,900)
c2.cd(1)
phiraw.SetTitle('Invariant Mass of Charged Kaons from ep->epK^{+}K^{-}; I.M. K^{+}K^{-} (GeV);counts')
phiraw.SetLineColor(kBlack)
phipall.SetLineColor(kRed)
#phipalllmbd.SetLineColor(kRed)
phiraw.Draw()
phipall.Draw('same')
l2 = TLegend(0.6,0.9, 0.9, 0.9)
l2.AddEntry(phiraw,'All Events')
l2.AddEntry(phipall,'Final Events')
l2.Draw('same')
#phipalllmbd.Draw('same')
c2.SaveAs('hist_phi_mass_peak_comparison.png')
