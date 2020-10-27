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



c0 =TCanvas('c0','c0',900,900)
c0.cd(1)
hhs['kp_mntm_phimass_pass_all_exclusivity'].Rebin2D(5,5)
hhs['kp_mntm_phimass_pass_all_exclusivity'].SetAxisRange(0,4.75,'X')
hhs['kp_mntm_phimass_pass_all_exclusivity'].SetAxisRange(0,3,'Y')
hhs['kp_mntm_phimass_pass_all_exclusivity'].Draw('colz')
c0.SaveAs('hist_kp_mntm_phi_mass_pass_all_epkpkm_outbending.png')

c1 =TCanvas('c1','c1',900,900)
c1.cd(1)
hhs['km_mntm_phimass_pass_all_exclusivity'].Rebin2D(5,5)
hhs['km_mntm_phimass_pass_all_exclusivity'].SetAxisRange(0,4.75,'X')
hhs['km_mntm_phimass_pass_all_exclusivity'].SetAxisRange(0,3,'Y')
hhs['km_mntm_phimass_pass_all_exclusivity'].Draw('colz')
c1.SaveAs('hist_km_mntm_phi_mass_pass_all_epkpkm_outbending.png')

c2 =TCanvas('c2','c2',900,900)
c2.cd(1)
hhs['phi_mntm_phimass_pass_all_exclusivity'].Rebin2D(5,5)
hhs['phi_mntm_phimass_pass_all_exclusivity'].SetAxisRange(0,4.75,'X')
hhs['phi_mntm_phimass_pass_all_exclusivity'].SetAxisRange(0,3,'Y')
hhs['phi_mntm_phimass_pass_all_exclusivity'].Draw('colz')
c2.SaveAs('hist_phi_mntm_phi_mass_pass_all_epkpkm_outbending.png')

