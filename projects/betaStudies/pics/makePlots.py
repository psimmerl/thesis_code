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


pr_b_thry = TF1("pr_b_thry","x/sqrt(x*x + 0.938*0.938)",0.10,8.0)
kp_b_thry = TF1("kp_b_thry","x/sqrt(x*x + 0.493*0.493)",0.10,8.0)
km_b_thry = TF1("km_b_thry","x/sqrt(x*x + 0.497*0.493)",0.10,8.0)

c1 = TCanvas("c1","c1",1200,450)
c1.Divide(3,1)
c1.cd(1)
hhs['pr_betap'].SetTitle("Proton Beta vs P For All Combinations; p (GeV); #beta")
hhs['pr_betap'].Draw('colz')
pr_b_thry.SetLineColor(kRed)
pr_b_thry.SetLineStyle(2)
pr_b_thry.Draw("same")

c1.cd(2)
hhs['kp_betap'].SetTitle("K^{+} Beta vs P For All Combinations; p (GeV); #beta")
hhs['kp_betap'].Draw('colz')
kp_b_thry.SetLineColor(kRed)
kp_b_thry.SetLineStyle(2)
kp_b_thry.Draw("same")

c1.cd(3)
hhs['km_betap'].SetTitle("K^{-} Beta vs P For All Combinations; p (GeV); #beta")
hhs['km_betap'].Draw('colz')
km_b_thry.SetLineColor(kRed)
km_b_thry.SetLineStyle(2)
km_b_thry.Draw("same")
c1.SaveAs("beta_analysis_all_hadrons_all_comb.png")


c2 = TCanvas("c2","c2",1200,450)
c2.Divide(3,1)
c2.cd(1)
hhs['pr_betap_final'].SetTitle("Proton Beta vs P For All Combinations with ALL Excl. Cuts; p (GeV); #beta")
hhs['pr_betap_final'].Draw('colz')
pr_b_thry.SetLineColor(kRed)
pr_b_thry.SetLineStyle(2)
pr_b_thry.Draw("same")

c2.cd(2)
hhs['kp_betap_final'].SetTitle("K^{+} Beta vs P For All Combinations  with ALL Excl. Cuts; p (GeV); #beta")
hhs['kp_betap_final'].Draw('colz')
kp_b_thry.SetLineColor(kRed)
kp_b_thry.SetLineStyle(2)
kp_b_thry.Draw("same")

c2.cd(3)
hhs['km_betap_final'].SetTitle("K^{-} Beta vs P For All Combinations with ALL Excl Cuts; p (GeV); #beta")
hhs['km_betap_final'].Draw('colz')
km_b_thry.SetLineColor(kRed)
km_b_thry.SetLineStyle(2)
km_b_thry.Draw("same")
c2.SaveAs("beta_analysis_all_hadrons_all_comb_allcuts.png")




    
