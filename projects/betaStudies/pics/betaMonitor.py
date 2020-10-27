from ROOT import TFile, TH1F, TH2F, TF1
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed
import sys,os

gStyle.SetOptStat(00000)

hhs = {}
runs = set()
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    if '_r' in obj.GetName():
        hhs[obj.GetName()] = obj
        #print(obj.GetName()[-4:])
        runs.add(int(obj.GetName()[-4:]))

pr_b_thry = TF1("pr_b_thry","x/sqrt(x*x + 0.938*0.938)",0.10,8.0)
kp_b_thry = TF1("kp_b_thry","x/sqrt(x*x + 0.493*0.493)",0.10,8.0)
pip_b_thry = TF1("pip_b_thry","x/sqrt(x*x + 0.139*0.139)",0.10,8.0)

pr_b_thry.SetLineColor(kRed)
kp_b_thry.SetLineColor(kRed)
pip_b_thry.SetLineColor(kRed)



#c0 = TCanvas('c0','c0',500*6,500*len(runs))
#c0.Divide(6,len(runs))


for rr in runs:
    print(' on run %s ' %(rr))
    c0 = TCanvas('c0'+str(rr),'c0'+str(rr),1000,1500)
    c0.Divide(2,3)
    cc=1

    for ss in range(1,7):
        c0.cd(cc)
        hhs['ftof_betap_neg_s'+str(ss)+'_r'+str(rr)].Draw('colz')
        #pr_b_thry.Draw('same')
        #kp_b_thry.Draw('same')
        #pip_b_thry.Draw('same')
    
        cc+=1


    c0.SaveAs('hist_betap_neg_per_sector_run'+str(rr)+'.pdf')
