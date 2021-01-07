import argparse
from ROOT import TFile, TH1F, TH2F, TF1, TLine, TLatex, TLegend, TBox
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kBlack, kBlue, kMagenta, kGreen
from ROOT import TGraph
import sys,os, math
from array import array

#cuts used are stored in this map indexed 0, 1, 2, 3
#exclusivity_cut_results = (np.absolute(epkpkmX.E())<0.6, epkpXmm2 < 0.6, ekpkmXmm2 < 1.2, epkmXmm2 < 0.6)
#pass cuts designated with pass_c012

ap = argparse.ArgumentParser()#Paul addition
ap.add_argument(
    '-i',
    '--input_file',
    required=True
)
ap.add_argument(
    '-o',
    '--output_prefix',
    required=True
)
args = ap.parse_args()

input_rootfile = args.input_file 
out_name = args.output_prefix 
ff = TFile(input_rootfile)#sys.argv[1])
#out_name='pidtype1_inb_V2'

gStyle.SetOptStat(0000)
lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)

hhs = {}
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj

c3a=TCanvas('c3a','c3a',1000,1000)
c3a.Divide(2,2)
c3a.cd(1)
hhs['missing_energy_combssize1'].SetTitle('')
hhs['missing_energy_combssize1'].SetLineColor(kBlack)
hhs['missing_energy_pass_all_but_missing_energy'].SetLineColor(kRed)
hhs['missing_energy_fail_at_least_missing_energy'].SetLineColor(kBlue)
hhs['missing_energy_combssize1'].Draw()
hhs['missing_energy_pass_all_but_missing_energy'].Draw('same')
hhs['missing_energy_fail_at_least_missing_energy'].Draw('same')
lab.DrawLatex(0.1, 0.925, 'Missing Energy: ep #rightarrow epK^{ +}K^{ -}X')
lab.DrawLatex(0.5, 0.015, 'Missing Energy (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'counts')
lab.SetTextAngle(0)

c3a.cd(2)
hhs['mm_ekpkmX_combssize1'].SetTitle('')
hhs['mm_ekpkmX_combssize1'].SetLineColor(kBlack)
hhs['mm_ekpkmX_pass_all_but_missing_proton'].SetLineColor(kRed)
hhs['mm_ekpkmX_fail_at_least_missing_proton'].SetLineColor(kBlue)
hhs['mm_ekpkmX_combssize1'].Draw()
hhs['mm_ekpkmX_pass_all_but_missing_proton'].Draw('same')
hhs['mm_ekpkmX_fail_at_least_missing_proton'].Draw('same')
hhs['mm_ekpkmX_combssize1'].SetMinimum(0)
lab.DrawLatex(0.1, 0.925, 'Missing Mass^{2}: ep #rightarrow eK^{ +}K^{ -}X')
lab.DrawLatex(0.5, 0.015, 'Proton Mass^{2} (GeV^{ 2})')
lab.SetTextAngle(90)
lab.DrawLatex(0.036, 0.5, 'counts')
lab.SetTextAngle(0)

c3a.cd(3)
#hhs['epkpX_raw'].SetTitle('epK^{+}X Raw and All but MM2 of epK^{+}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['mm_epkpX_combssize1'].SetTitle('')
hhs['mm_epkpX_combssize1'].SetLineColor(kBlack)
hhs['mm_epkpX_pass_all_but_missing_kaonM'].SetLineColor(kRed)
hhs['mm_epkpX_fail_at_least_missing_kaonM'].SetLineColor(kBlue)
hhs['mm_epkpX_combssize1'].Draw()
hhs['mm_epkpX_pass_all_but_missing_kaonM'].Draw('same')
hhs['mm_epkpX_fail_at_least_missing_kaonM'].Draw('same')
hhs['mm_epkpX_fail_at_least_missing_kaonM'].SetMinimum(0)
lab.DrawLatex(0.1, 0.925, 'Missing Mass^{2}: ep #rightarrow epK^{ +}X')
lab.DrawLatex(0.5, 0.015, 'K^{ -} Mass^{2} (GeV^{ 2})')
lab.SetTextAngle(90)
lab.DrawLatex(0.036, 0.5, 'counts')
lab.SetTextAngle(0)
c3a.cd(4)
#hhs['epkmX_raw'].SetTitle('epK^{-}X Raw and All but MM2 of epK^{-}X cut;Missing Mass^{2} (GeV^{2});counts')
hhs['mm_epkmX_combssize1'].SetTitle('')
hhs['mm_epkmX_combssize1'].SetLineColor(kBlack)
hhs['mm_epkmX_pass_all_but_missing_kaonP'].SetLineColor(kRed)
hhs['mm_epkmX_fail_at_least_missing_kaonP'].SetLineColor(kBlue)
hhs['mm_epkmX_combssize1'].Draw()
hhs['mm_epkmX_pass_all_but_missing_kaonP'].Draw('same')
hhs['mm_epkmX_fail_at_least_missing_kaonP'].Draw('same')
hhs['mm_epkmX_fail_at_least_missing_kaonP'].SetMinimum(0)
lab.DrawLatex(0.1, 0.925, 'Missing Mass^{2}: ep #rightarrow epK^{ -}X')
lab.DrawLatex(0.5, 0.015, 'K^{ +} Mass^{2} (GeV^{ 2})')
lab.SetTextAngle(90)
lab.DrawLatex(0.036, 0.5, 'counts')
lab.SetTextAngle(0)
c3a.SaveAs('hist_all_before_after_passfail_{}.pdf'.format(out_name))

################################################################################################################################################
## draw the missing e results with cuts on all other variables
## 
c3b=TCanvas('c3b','c3b',1000,1000)
c3b.Divide(2,2)
c3b.cd(1)
hhs['missing_energy_combssize1'].Draw()
hhs['missing_energy_pass_all_but_missing_energy'].Draw('same')
hhs['missing_energy_fail_at_least_missing_energy'].Draw('same')
lab.DrawLatex(0.1, 0.925, 'Missing Energy: ep #rightarrow epK^{ +}K^{ -}X')
lab.DrawLatex(0.5, 0.015, 'Energy (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.038, 0.5, 'counts')
lab.SetTextAngle(0)
leg  = TLegend(0.12, 0.68, 0.48, 0.88)
leg.AddEntry(hhs['missing_energy_combssize1'],'No Cuts','l')
leg.AddEntry(hhs['missing_energy_pass_all_but_missing_energy'],'Pass Cuts','l')
leg.AddEntry(hhs['missing_energy_fail_at_least_missing_energy'],'Failed Cuts','l')
leg.SetLineWidth(0)
leg.Draw('same')

c3b.cd(2)
hhs['mm_ekpkmX_combssize1'].Draw()
hhs['mm_ekpkmX_combssize1'].SetMinimum(0)
fit_miss_p = TF1("fit_miss_p","gaus(0)",0.83, 0.95); 
fit_miss_p.SetParLimits(3,-10,10)
fit_miss_p.SetParLimits(4,-10,10)
fit_miss_p.SetParLimits(5,-10,10)
hhs['mm_ekpkmX_combssize1'].Fit("fit_miss_p","R")
fit_miss_p.SetLineWidth(2)
lab.DrawLatex(0.65, 0.8, '#mu:{:.4f}'.format(math.sqrt((fit_miss_p.GetParameter(1)))))
lab.DrawLatex(0.65, 0.75, '#sigma:{:.4f}'.format(math.sqrt(fit_miss_p.GetParameter(2))))
pro_limits  = open('excl_phi_mm2_pro_limits_{}.txt'.format(out_name),'w')
#pro_limits.write(str(fit_miss_p.GetParameter(1)) + ' ' + str(fit_miss_p.GetParameter(2)) )
pro_limits.close()
fit_miss_p.Draw('same')
#draw the cut lines
miss_p_min = fit_miss_p.GetParameter(1) - 5*fit_miss_p.GetParameter(2)
miss_p_max = fit_miss_p.GetParameter(1) + 5*fit_miss_p.GetParameter(2)
c3b.Update()
ymax_miss_p = gPad.GetUymax() 
l_miss_p_min = TLine(miss_p_min, 0.0, miss_p_min, ymax_miss_p)
l_miss_p_max = TLine(miss_p_max, 0.0, miss_p_max, ymax_miss_p)
l_miss_p_min.SetLineWidth(2)
l_miss_p_max.SetLineWidth(2)
l_miss_p_min.SetLineColor(kRed)
l_miss_p_max.SetLineColor(kRed)
l_miss_p_min.Draw('same')
l_miss_p_max.Draw('same')

lab.DrawLatex(0.1, 0.925, 'Missing Mass^{2}: ep #rightarrow eK^{ +}K^{ -}X')
lab.DrawLatex(0.5, 0.015, 'Proton Mass^{2} (GeV^{ 2})')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'counts')
lab.SetTextAngle(0)


c3b.cd(3)
hhs['mm_epkpX_combssize1'].Draw()
fit_miss_km = TF1("fit_miss_km","gaus(0)", 0.18, 0.3);
fit_miss_km.SetLineWidth(2)
hhs['mm_epkpX_combssize1'].Fit("fit_miss_km","R")
fit_miss_km.Draw('same')
lab.DrawLatex(0.70, 0.8, '#mu:{:.4f}'.format(math.sqrt(fit_miss_km.GetParameter(1))))
lab.DrawLatex(0.70, 0.75, '#sigma:{:.4f}'.format(math.sqrt(fit_miss_km.GetParameter(2))))
km_limits  = open('excl_phi_mm2_km_limits_{}.txt'.format(out_name),'w')
#km_limits.write(str(fit_miss_km.GetParameter(1)) + ' ' + str(fit_miss_km.GetParameter(2)) )
km_limits.close()
#missing km mass
miss_km_min = fit_miss_km.GetParameter(1) - 5*fit_miss_km.GetParameter(2)
miss_km_max = fit_miss_km.GetParameter(1) + 5*fit_miss_km.GetParameter(2)
c3b.Update()
ymax_miss_km = gPad.GetUymax() 
l_miss_km_min = TLine(miss_km_min, 0.0, miss_km_min, ymax_miss_km)
l_miss_km_max = TLine(miss_km_max, 0.0, miss_km_max, ymax_miss_km)
l_miss_km_min.SetLineWidth(2)
l_miss_km_max.SetLineWidth(2)
l_miss_km_min.SetLineColor(kRed)
l_miss_km_max.SetLineColor(kRed)
l_miss_km_min.Draw('same')
l_miss_km_max.Draw('same')


lab.DrawLatex(0.1, 0.925, 'Missing Mass^{2}: ep #rightarrow epK^{ +}X')
lab.DrawLatex(0.5, 0.015, 'K^{ -} Mass^{2} (GeV^{ 2})')
lab.SetTextAngle(90)
lab.DrawLatex(0.036, 0.5, 'counts')
lab.SetTextAngle(0)

c3b.cd(4)
hhs['mm_epkmX_combssize1'].Draw()
fit_miss_kp = TF1("fit_miss_kp","gaus(0)", 0.18, 0.3);
hhs['mm_epkmX_combssize1'].Fit("fit_miss_kp","R")
fit_miss_kp.SetLineWidth(2)
fit_miss_kp.Draw('same')
lab.DrawLatex(0.70, 0.8, '#mu:{:.4f}'.format(math.sqrt(fit_miss_kp.GetParameter(1))))
lab.DrawLatex(0.70, 0.75, '#sigma:{:.4f}'.format(math.sqrt(fit_miss_kp.GetParameter(2))))
kp_limits  = open('excl_phi_mm2_kp_limits_{}.txt'.format(out_name),'w')
#kp_limits.write(str(fit_miss_kp.GetParameter(1)) + ' ' + str(fit_miss_kp.GetParameter(2)) )
kp_limits.close()
miss_kp_min = fit_miss_kp.GetParameter(1) - 5*fit_miss_kp.GetParameter(2)
miss_kp_max = fit_miss_kp.GetParameter(1) + 5*fit_miss_kp.GetParameter(2)
c3b.Update()
ymax_miss_kp = gPad.GetUymax() 
l_miss_kp_min = TLine(miss_kp_min, 0.0, miss_kp_min, ymax_miss_kp)
l_miss_kp_max = TLine(miss_kp_max, 0.0, miss_kp_max, ymax_miss_kp)
l_miss_kp_min.SetLineWidth(2)
l_miss_kp_max.SetLineWidth(2)
l_miss_kp_min.SetLineColor(kRed)
l_miss_kp_max.SetLineColor(kRed)
l_miss_kp_min.Draw('same')
l_miss_kp_max.Draw('same')


lab.DrawLatex(0.1, 0.925, 'Missing Mass^{2}: ep #rightarrow epK^{ -}X')
lab.DrawLatex(0.5, 0.015, 'K^{ +} Mass^{2} (GeV^{ 2})')
lab.SetTextAngle(90)
lab.DrawLatex(0.036, 0.5, 'counts')
lab.SetTextAngle(0)
c3b.SaveAs('hist_all_before_after_passfail_missing_energy_{}.pdf'.format(out_name))


fit_miss_e = TF1("fit_miss_e","gaus(0)", -0.1, 0.1)
hhs['missing_energy_combssize1'].Fit('fit_miss_e','R')
me_limits  = open('excl_phi_me_limits_{}.txt'.format(out_name),'w')
#me_limits.write(str(fit_miss_e.GetParameter(1)) + ' ' + str(fit_miss_e.GetParameter(2)) )
me_limits.close()



'''

#final comparison of phi plots
c4a = TCanvas('phimass0','phimass0',900,900)
c4a.Divide(1,1)
c4a.cd(1)
hhs['vphimass_raw'].SetTitle('')
hhs['vphimass_raw'].SetLineColor(kGreen)
hhs['vphimass_raw'].SetFillColor(kGreen)
hhs['vphimass_raw'].SetFillStyle(3001)
hhs['vphimass_raw'].Draw()
c4a.Update()
ymax = gPad.GetUymax()
l_phi_mass = TLine(1.019, 0.0, 1.019, ymax)
l_phi_mass.SetLineColor(kRed)
l_phi_mass.SetLineWidth(3)
l_phi_mass.Draw('same')
lab.DrawLatex(0.45, 0.725, 'M_{#phi} = 1.019 GeV')
lab.DrawLatex(0.4, 0.015, 'Invariant Mass K^{+}K^{-} (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'Counts')
lab.SetTextAngle(0)
c4a.SaveAs("phi_mass_raw.pdf")

c4 = TCanvas('phimass','phimass',900,900)
c4.Divide(1,1)
c4.cd(1)
#gPad.SetLogy()
hhs['vphimass_raw'].SetTitle('')
hhs['vphimass_raw'].SetLineColor(634)
hhs['vphimass_raw'].SetFillColor(634)
hhs['vphimass_raw'].SetFillStyle(3001)
hhs['vphimass_pass_all_exclusivity'].SetLineColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetFillColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetFillStyle(3001)
hhs['vphimass_raw'].Draw()
hhs['vphimass_pass_all_exclusivity'].Draw('same')
hhs['vphimass_pass_all_exclusivity'].SetMinimum(0.00001)

leg = TLegend(0.7, 0.7, 0.88, 0.88)
leg.SetBorderSize(0)
leg.AddEntry(hhs['vphimass_raw'],"No Cuts","f")
leg.AddEntry(hhs['vphimass_pass_all_exclusivity'],"All Cuts","f")
leg.Draw()
lab.DrawLatex(0.4, 0.015, 'Invariant Mass K^{+}K^{-} (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'Counts')
lab.SetTextAngle(0)
c4.SaveAs('phi_mass_before_after.pdf')

c5 = TCanvas('c5','c5',900,900)
c5.Divide(1,1)
c5.cd(1)
hhs['vphimass_pass_all_exclusivity'].SetTitle('')
hhs['vphimass_pass_all_exclusivity'].SetLineColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetFillColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetFillStyle(3001)
hhs['vphimass_pass_all_exclusivity'].Draw()
c5.Update()
ymax = gPad.GetUymax()
phi_mass_fit=1.020
phi_mass_resolution = 4.25752e-03  
max_phi_mass_signal_limit = phi_mass_fit + 2.5*phi_mass_resolution
min_phi_mass_signal_limit = phi_mass_fit - 2.5*phi_mass_resolution

select_sig = TBox(0.982, 0.0, 1.05, ymax)#max_phi_mass_signal_limit)
select_sig.SetFillColorAlpha(kGreen,0.25)
select_sig.Draw('same')
lab.DrawLatex(0.4, 0.015, 'Invariant Mass K^{+}K^{-} (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'Counts')
lab.SetTextAngle(0)
c5.SaveAs('selected_phi_mass_signal.pdf')


c6 = TCanvas('c6','c6',900,900)
c6.Divide(1,1)
c6.cd(1)
hhs['vphimass_pass_all_exclusivity'].SetTitle('')
hhs['vphimass_pass_all_exclusivity'].SetLineColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetFillColor(kBlack)
hhs['vphimass_pass_all_exclusivity'].SetFillStyle(3001)
hhs['vphimass_pass_all_exclusivity'].Draw()
c6.Update()
ymax = gPad.GetUymax()

select_bck = TBox(1.1, 0.0, 1.17, ymax)#max_phi_mass_signal_limit)
select_bck.SetFillColorAlpha(kMagenta,0.25)
select_bck.Draw('same')
lab.DrawLatex(0.4, 0.015, 'Invariant Mass K^{+}K^{-} (GeV)')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'Counts')
lab.SetTextAngle(0)
lab.DrawLatex(0.6, 0.5, '1.1 - 1.17 GeV')
c6.SaveAs('selected_phi_mass_bck.pdf')



'''
