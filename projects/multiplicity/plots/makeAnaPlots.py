from ROOT import TH1F, TFile, TCanvas
from ROOT import gStyle, gPad
import sys


fin = TFile(sys.argv[1],'open')

h_all_ev_comb = fin.Get('all_combination_all_events')
h_all_ev_comb_cd = fin.Get('CD_combination_all_events')
h_all_ev_comb_fd = fin.Get('FD_combination_all_events')
h_all_final_comb = fin.Get('all_combination_final_events')
h_all_final_cd = fin.Get('CD_combination_final_events')
h_all_final_fd = fin.Get('FD_combination_final_events')

h_npr_fd = fin.Get('proton_counts_all_events_FD')
h_npr_cd = fin.Get('proton_counts_all_events_CD')

h_nkp_fd = fin.Get('kp_counts_all_events_FD')
h_nkp_cd = fin.Get('kp_counts_all_events_CD')

h_nkm_fd = fin.Get('km_counts_all_events_FD')
h_nkm_cd = fin.Get('km_counts_all_events_CD')

c0 = TCanvas('c0','c0',1200,600)
c0.Divide(2,1)
c0.cd(1)
h_npr_fd.SetTitle('Proton Multiplicity Per Event In FD;frequency;counts')
h_npr_fd.Draw()
c0.cd(2)
h_npr_cd.SetTitle('Proton Multiplicity Per Event In CD;frequency;counts')
h_npr_cd.Draw()
c0.SaveAs('hadron_pr_multi.png')

c0a = TCanvas('c0a','c0a',1200,600)
c0a.Divide(2,1)
c0a.cd(1)
h_nkp_fd.SetTitle('K^{+} Multiplicity Per Event In FD;frequency;counts')
h_nkp_fd.Draw()
c0a.cd(2)
h_nkp_cd.SetTitle('K^{+} Multiplicity Per Event In CD;frequency;counts')
h_nkp_cd.Draw()
c0a.SaveAs('hadron_kp_multi.png')

c0b = TCanvas('c0b','c0b',1200,600)
c0b.Divide(2,1)
c0b.cd(1)
h_nkm_fd.SetTitle('K^{-} Multiplicity Per Event In FD;frequency;counts')
h_npr_fd.Draw()
c0b.cd(2)
h_nkm_cd.SetTitle('K^{-} Multiplicity Per Event In CD;frequency;counts')
h_nkm_cd.Draw()
c0b.SaveAs('hadron_km_multi.png')

gStyle.SetOptStat(0000)

c1=TCanvas("c1","c1",900,900)
c1.cd(1)
gPad.SetLogy()
h_all_ev_comb.SetTitle("Total Combinations of epK^{+}K^{-};number of combinations;counts")
h_all_ev_comb.Draw()
c1.SaveAs("h_all_ev_comb.png")

c2=TCanvas("c2","c2",900,900)
c2.cd(1)
gPad.SetLogy()
h_all_ev_comb_cd.SetTitle("Total Combinations of epK^{+}K^{-} in CD;number of combinations;counts")
h_all_ev_comb_cd.Draw()
c2.SaveAs("h_all_ev_comb_cd.png")

c3=TCanvas("c3","c3",900,900)
c3.cd(1)
gPad.SetLogy()
h_all_ev_comb_fd.SetTitle("Total Combinations of epK^{+}K^{-} in FD;number of combinations;counts")
h_all_ev_comb_fd.Draw()
c3.SaveAs("h_all_ev_comb_fd.png")


c4=TCanvas("c4","c4",900,900)
c4.cd(1)
gPad.SetLogy()
h_all_final_comb.SetTitle("Total Combinations After Excl. Cuts of epK^{+}K^{-};number of combinations;counts")
h_all_final_comb.Draw()
c4.SaveAs("h_all_final_comb.png")

c5=TCanvas("c5","c5",900,900)
c5.cd(1)
gPad.SetLogy()
h_all_final_cd.SetTitle("Total Combinations After Excl. Cuts of epK^{+}K^{-} in CD;number of combinations;counts")
h_all_final_cd.Draw()
c5.SaveAs("h_all_final_comb_cd.png")

c6=TCanvas("c6","c6",900,900)
c6.cd(1)
gPad.SetLogy()
h_all_final_fd.SetTitle("Total Combinations After Excl. Cuts of epK^{+}K^{-} in FD;number of combinations;counts")
h_all_final_fd.Draw()
c6.SaveAs("h_all_final_comb_fd.png")

