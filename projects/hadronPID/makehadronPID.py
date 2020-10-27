from ROOT import TFile, TH1F, TH2F, TF1
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, TLatex
import sys,os

gStyle.SetOptStat(00000)
gStyle.SetPalette(55)

beam=10.604

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def setup_global_options():
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetPalette(57)

setup_global_options()

ff = TFile(sys.argv[1])

hhs= load_histos(ff)



lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)


pro = TF1("pro","x/TMath::Sqrt(x*x + 0.938*0.938)", 0.2, 6)
pro.SetLineColor(kRed)
pro.SetLineWidth(4)
kaon = TF1("kaon","x/TMath::Sqrt(x*x + 0.493*0.493)", 0.2, 6)
kaon.SetLineColor(kRed)
kaon.SetLineWidth(4)
pion = TF1("pion","x/TMath::Sqrt(x*x + 0.139*0.139)", 0.2, 6)
pion.SetLineColor(kRed)
pion.SetLineWidth(4)

hposbet = hhs["p_beta_pos_FTOF"]
hnegbet = hhs["p_beta_neg_FTOF"]

cn = TCanvas('cn','cn',900,900)
cn.cd(1)
gPad.SetLogz()
hposbet.SetTitle("")
hposbet.GetXaxis().SetRangeUser(0,5)
hposbet.GetYaxis().SetRangeUser(0.85,1.1)
hposbet.Draw('colz')
pro.Draw("same")
kaon.Draw("same")
pion.Draw("same")
lab.DrawLatex(0.1, 0.925,"#beta vs Momentum")
lab.DrawLatex(0.5, 0.015,"Momentum (GeV)")
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, "#beta")
lab.SetTextAngle(0)

cn.SaveAs('h_bp_posV1.pdf')


cp = TCanvas('cp','cp',900,900)
cp.cd(1)
gPad.SetLogz()
hnegbet.SetTitle("")
hnegbet.GetXaxis().SetRangeUser(0,5)
hnegbet.GetYaxis().SetRangeUser(0.85,1.1)
hnegbet.Draw('colz')
kaon.Draw("same")
pion.Draw("same")
lab.DrawLatex(0.1, 0.925,"#beta vs Momentum")
lab.DrawLatex(0.5, 0.015,"Momentum (GeV)")
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, "#beta")
lab.SetTextAngle(0)
cp.SaveAs('h_bp_negV1.pdf')


