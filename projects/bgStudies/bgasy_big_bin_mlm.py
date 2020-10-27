from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TLatex
from ROOT import TGraphErrors, TBox
from ROOT import gROOT, gBenchmark, gStyle, gPad, kRed, kMagenta
from array import array
import math
import sys


bg_center = array('d')
bg_asy = array('d')
bg_center_err = array('d')
bg_asy_err = array('d')

deltamass = 0.07
asy = array('d')
asy.append(0.0624)
asy.append(0.0227)
asy.append(-0.0146)
asy.append(0.00947)
asy.append(0.1065)
asy.append(0.02037)
asyerr = array('d')
asyerr.append(0.0331)
asyerr.append(0.03711)
asyerr.append(0.03855)
asyerr.append(0.03831)
asyerr.append(0.0420)
asyerr.append(0.0460)


bins = 6
center_start = 1.027
for ii in range(0, bins):
    temp_center = (center_start +  (center_start + deltamass))/2.0
    print(temp_center)
    bg_center.append(temp_center)
    center_start += deltamass
    bg_center_err.append(0)
    bg_asy.append(asy[ii])
    bg_asy_err.append(asyerr[ii])
    

lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)

fit_bg_bsa = TF1('fit_bg_bsa','[0]*x + [1]',1.065, 1.415)

g_asy = TGraphErrors(len(bg_center), bg_center, bg_asy, bg_center_err, bg_asy_err)
c1 = TCanvas('c1','c1',900,900)
c1.Divide(1,1)
c1.cd(1)
g_asy.SetTitle('')
g_asy.SetMarkerSize(6)
g_asy.SetMarkerColor(kMagenta)
g_asy.Draw('AP')
lab.DrawLatex(0.1, 0.925, 'Background Asymmetry Contribution')
lab.DrawLatex(0.45, 0.015, 'Invariant Mass K^{ +}K^{ -} (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.035, 0.5, 'Asymmetry')
lab.SetTextAngle(0)

g_asy.Fit('fit_bg_bsa')
fit_bg_bsa.Draw('same')

c1.SaveAs('background_asy_mlm.pdf')
