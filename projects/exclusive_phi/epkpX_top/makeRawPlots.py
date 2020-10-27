from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack
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
c0.Divide(2,3)
c0.cd(1)
hhs['LeX_raw'].Draw()
c0.cd(2)
hhs['LepX_raw'].SetAxisRange(-1., 5.,'X');
hhs['LepX_raw'].Draw()
c0.cd(3)
hhs['LekpX_raw'].SetAxisRange(-1., 5.,'X');
hhs['LekpX_raw'].Draw()
c0.cd(4)
hhs['LekpkmX_raw'].SetAxisRange(-1., 5.,'X');
hhs['LekpkmX_raw'].Draw()
c0.cd(5)
hhs['LepkpX_raw'].SetAxisRange(-1., 2.5,'X');
hhs['LepkpX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['LepkpX_raw'].Draw()
hhs['LepkpX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c0.cd(6)
hhs['LepkpkmXe_raw'].Draw()
c0.SaveAs('hist_epkpXtop_Lraw_eXs.png')


c1 = TCanvas('c1','c1',1000,1500)
c1.Divide(2,2)
c1.cd(1)
hhs['Lvphimass_raw'].Draw()
c1.cd(2)
hhs['lambdamass_raw'].Draw()
c1.cd(3)
hhs['lambdaphimass_raw'].Draw('colz')
c1.SaveAs('hist_phi_lambda_mass_raw_epkpXtop.png')

c2 = TCanvas('c2','c2',1000,1000)
c2.cd(1)
hhs['kp_as_pip_raw'].SetTitle('epK^{+}X with Kaon as pip plus MM; epK^{+}X MM; Pion Seperation')
hhs['kp_as_pip_raw'].Draw('colz')
c2.SaveAs('hist_kp_as_pip_raw_epkpXtop.png')


c3 = TCanvas('c3','c3',1000,1000)
c3.Divide(2,1)
c3.cd(1)
hhs['cplkp_raw'].Draw()
c3.cd(2)
#hhs['cplpr_raw'].Draw()
c3.SaveAs('hist_coplan_angle_kppr_epkpXtop.png')


c4 = TCanvas("c4","c4",1200,450)
c4.Divide(2,3)
c4.cd(1)
hhs['LeX_raw'].Draw()
hhs['LeX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['LeX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c4.cd(2)
hhs['LepX_raw'].SetAxisRange(-1., 5.,'X');
hhs['LepX_raw'].Draw()
hhs['LepX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['LepX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c4.cd(3)
hhs['LekpX_raw'].SetAxisRange(-1., 5.,'X');
hhs['LekpX_raw'].Draw()
hhs['LekpX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['LekpX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c4.cd(4)
hhs['LekpkmX_raw'].SetAxisRange(-1., 5.,'X');
hhs['LekpkmX_raw'].Draw()
hhs['LekpkmX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['LekpkmX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c4.cd(5)
hhs['LepkpX_raw'].SetAxisRange(-1., 2.5,'X');
hhs['LepkpX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['LepkpX_raw'].Draw()
hhs['LepkpX_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c4.cd(6)
hhs['LepkpkmXe_raw'].Draw()
c4.SaveAs('hist_epkpXtop_beforeafter_eXs.png')

c5 = TCanvas('c5','c5',1000,500)
c5.Divide(2,1)
c5.cd(1)
hhs['Lvphimass_raw'].Draw()
hhs['Lvphimass_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['Lvphimass_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c5.cd(2)
hhs['lambdamass_raw'].Draw()
hhs['lambdamass_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].SetLineColor(kRed)
hhs['lambdamass_pass_modified_exclusive_cuts_epkpX_cutepkpX_cut'].Draw('same')
c5.SaveAs('hist_phi_lambda_mass_beforeafter_epkpXtop.png')
