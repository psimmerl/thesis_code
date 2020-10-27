from ROOT import TFile, TH1F, TH2F, TF1
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed
import sys,os

gStyle.SetOptStat(00000)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


pip_b_thry = TF1("pr_b_thry","x/sqrt(x*x + 0.135*0.135)",0.10,8.0)
pr_b_thry = TF1("pr_b_thry","x/sqrt(x*x + 0.938*0.938)",0.10,8.0)
kp_b_thry = TF1("kp_b_thry","x/sqrt(x*x + 0.493*0.493)",0.10,8.0)
km_b_thry = TF1("km_b_thry","x/sqrt(x*x + 0.493*0.493)",0.10,8.0)

c1 = TCanvas("c1","c1",1200,1200)
c1.cd(1)
gPad.SetLogz()
hhs['ftof_betap_pos'].SetTitle("Beta vs P For All Positives; p (GeV); #beta")
hhs['ftof_betap_pos'].Draw('colz')
hhs['ftof_betap_pos'].GetYaxis().SetRangeUser(0.92,1.04)
hhs['ftof_betap_pos'].Draw('colz')

pip_b_thry.SetLineColor(kRed)
pip_b_thry.SetLineStyle(2)
pip_b_thry.SetLineWidth(4)
pip_b_thry.Draw('same')
pr_b_thry.SetLineColor(kRed)
pr_b_thry.SetLineStyle(2)
pr_b_thry.SetLineWidth(4)
pr_b_thry.Draw("same")
kp_b_thry.SetLineColor(kRed)
kp_b_thry.SetLineStyle(2)
kp_b_thry.SetLineWidth(4)
kp_b_thry.Draw("same")
c1.SaveAs('h_betap_all_positives.png')


c2 = TCanvas("c2","c2",1200,1200)
c2.cd(1)
gPad.SetLogz()
hhs['ftof_betap_neg'].SetTitle("Beta vs P For All Negatives; p (GeV); #beta")
hhs['ftof_betap_neg'].Draw('colz')
hhs['ftof_betap_neg'].GetYaxis().SetRangeUser(0.92,1.04)
hhs['ftof_betap_neg'].Draw('colz')
km_b_thry.SetLineColor(kRed)
km_b_thry.SetLineStyle(2)
km_b_thry.Draw("same")
c2.SaveAs('h_betap_all_negatives.png')
