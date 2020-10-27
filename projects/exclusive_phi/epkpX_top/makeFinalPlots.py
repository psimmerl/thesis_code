from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue
from ROOT import TGraph
import sys,os
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

#gStyle.SetOptStat(1111)

hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


c0 = TCanvas("c0","c0",1200,450)
c0.cd(1)
hhs['epkpkmXe_pass_c123'].SetTitle('epK^{+}K^{-}X Energy With All Cuts;Missing Energy (GeV);counts')
hhs['epkpkmXe_raw'].SetLineColor(kBlack)
hhs['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs['epkpkmXe_fail_c123'].SetLineColor(kBlue)
hhs['epkpkmXe_raw'].Draw()
hhs['epkpkmXe_fail_c123'].Draw('same')
hhs['epkpkmXe_pass_c123'].Draw('same')
c0.SaveAs('hist_epkpkmXe_c123_passfail.png')

c1 = TCanvas("c1","c1",1200,450)
c1.cd(1)
hhs['ekpkmX_pass_c013'].SetTitle('eK^{+}K^{-}X MM2 Cuts;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_raw'].SetLineColor(kBlack)
hhs['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs['ekpkmX_fail_c013'].SetLineColor(kBlue)
hhs['ekpkmX_raw'].Draw()
hhs['ekpkmX_fail_c013'].Draw('same')
hhs['ekpkmX_pass_c013'].Draw('same')
hhs['ekpkmX_raw'].SetMinimum(0)
c1.SaveAs('hist_ekpkmXe_c013_passfail.png')


c2 = TCanvas("c2","c2",1200,450)
c2.cd(1)
hhs['epkpX_pass_c023'].SetTitle('epK^{+}X MM2 Cuts;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_raw'].SetLineColor(kBlack)
hhs['epkpX_pass_c023'].SetLineColor(kRed)
hhs['epkpX_fail_c023'].SetLineColor(kBlue)
hhs['epkpX_raw'].Draw()
hhs['epkpX_pass_c023'].Draw('same')
hhs['epkpX_fail_c023'].Draw('same')
hhs['epkpX_pass_c023'].SetMinimum(0)
c2.SaveAs('hist_epkpX_c023_passfail.png')

c3 = TCanvas("c3","c3",1200,450)
c3.cd(1)
hhs['epkmX_pass_c012'].SetTitle('epK^{-}X MM2 Cuts;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_raw'].SetLineColor(kBlack)
hhs['epkmX_pass_c012'].SetLineColor(kRed)
hhs['epkmX_fail_c012'].SetLineColor(kBlue)
hhs['epkmX_raw'].Draw()
hhs['epkmX_pass_c012'].Draw('same')
hhs['epkmX_fail_c012'].Draw('same')
hhs['epkmX_pass_c012'].SetMinimum(0)
c3.SaveAs('hist_epkmX_c012_passfail.png')

c3a=TCanvas('c3a','c3a',1000,1000)
c3a.Divide(2,2)
c3a.cd(1)
hhs['epkpkmXe_pass_c123'].SetTitle('epK^{+}K^{-}X E Raw and All but ME cut;Missing Energy (GeV);counts')
hhs['epkpkmXe_raw'].SetLineColor(kBlack)
hhs['epkpkmXe_pass_c123'].SetLineColor(kRed)
hhs['epkpkmXe_fail_c123'].SetLineColor(kBlue)
hhs['epkpkmXe_raw'].Draw()
hhs['epkpkmXe_pass_c123'].Draw('same')
hhs['epkpkmXe_fail_c123'].Draw('same')
c3a.cd(2)
hhs['ekpkmX_pass_c013'].SetTitle('eK^{+}K^{-}X Raw and All but MM2 of eK^{+}K^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['ekpkmX_raw'].SetLineColor(kBlack)
hhs['ekpkmX_pass_c013'].SetLineColor(kRed)
hhs['ekpkmX_fail_c013'].SetLineColor(kBlue)
hhs['ekpkmX_raw'].Draw()
hhs['ekpkmX_pass_c013'].Draw('same')
hhs['ekpkmX_fail_c013'].Draw('same')
hhs['ekpkmX_raw'].SetMinimum(0)
c3a.cd(3)
hhs['epkpX_pass_c023'].SetTitle('epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkpX_raw'].SetLineColor(kBlack)
hhs['epkpX_pass_c023'].SetLineColor(kRed)
hhs['epkpX_fail_c023'].SetLineColor(kBlue)
hhs['epkpX_raw'].Draw()
hhs['epkpX_pass_c023'].Draw('same')
hhs['epkpX_fail_c023'].Draw('same')
hhs['epkpX_fail_c023'].SetMinimum(0)
c3a.cd(4)
hhs['epkmX_raw'].SetTitle('epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['epkmX_raw'].SetLineColor(kBlack)
hhs['epkmX_pass_c012'].SetLineColor(kRed)
hhs['epkmX_fail_c012'].SetLineColor(kBlue)
hhs['epkmX_raw'].Draw()
hhs['epkmX_pass_c012'].Draw('same')
hhs['epkmX_fail_c012'].Draw('same')
hhs['epkmX_fail_c012'].SetMinimum(0)
c3a.SaveAs('hist_all_before_after_passfail.png')
