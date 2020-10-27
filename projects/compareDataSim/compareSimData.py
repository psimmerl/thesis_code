from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue
import numpy as np

from array import array
import math
import sys,os

fd = TFile(sys.argv[1])
fs = TFile(sys.argv[2])
datatype = sys.argv[3]

base=os.path.basename(fd.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_DataVsSim_kinematics'+datatype+'.pdf'
can = TCanvas('can','can',1200,600)
can.Print('{}['.format(mon_out_file_name))
can.Divide(2,1)

gStyle.SetOptStat(00000)

#########################################################
#########################################################
# code to compare reconstructed DATA to reconstructed SIM
#########################################################
#########################################################

# compare electron kinematics - p, theta,
select_phi='_select_phi'
hd_el_p = fd.Get('el_p_pass_all_exclusivity'+select_phi)
hd_el_theta = fd.Get('el_theta_pass_all_exclusivity'+select_phi)

hs_el_p = fs.Get('el_p')
hs_el_theta = fs.Get('el_theta')

can.cd(1)
hd_el_p.SetLineColor(kRed)
hs_el_p.SetLineColor(kBlue)
hs_el_p.Scale( hd_el_p.GetMaximum() / hs_el_p.GetMaximum())
hd_el_p.SetTitle("Electron Mntm; p (GeV); counts")
hd_el_p.Draw()
hs_el_p.Draw("same+hist")
lel1=TLegend(0.7, 0.6, 0.9, 0.75)
lel1.AddEntry(hd_el_p,'DATA')
lel1.AddEntry(hs_el_p,'SIM')
lel1.Draw('same')

can.cd(2)
hd_el_theta.SetLineColor(kRed)
hs_el_theta.SetLineColor(kBlue)
hs_el_theta.Scale( hd_el_theta.GetMaximum() / hs_el_theta.GetMaximum())
hd_el_theta.SetTitle("Electron #theta; #theta (deg); counts")
hd_el_theta.Draw()
hs_el_theta.Draw("same+hist")
lel2=TLegend(0.7, 0.6, 0.9, 0.75)
lel2.AddEntry(hd_el_theta,'DATA')
lel2.AddEntry(hs_el_theta,'SIM')
lel2.Draw('same')

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,2)

######################################### compare kaonN kinematics - p, theta,
hd_el_ptheta = fd.Get('el_ptheta_pass_all_exclusivity'+select_phi)
hd_el_thetaphi = fd.Get('el_phitheta_pass_all_exclusivity'+select_phi)

hs_el_ptheta = fs.Get('el_ptheta')
hs_el_thetaphi = fs.Get('el_thetaphi')

can.cd(1)
hd_el_ptheta.SetTitle("DATA El. Mntm vs #theta; p (GeV); #theta (deg)")
hd_el_ptheta.Draw("colz")
can.cd(2)
hd_el_thetaphi.SetTitle("DATA El. #phi vs #theta; #phi (deg); #theta (deg)")
hd_el_thetaphi.Draw("colz")
can.cd(3)
hs_el_ptheta.SetTitle("SIM. El. Mntm vs #theta; p (GeV); #theta (deg)")
hs_el_ptheta.Draw("colz")
can.cd(4)
hs_el_thetaphi.SetTitle("SIM. El. #phi vs #theta; #phi (deg); #theta (deg)")
hs_el_thetaphi.Draw("colz")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)

########################## compare protons kinematics - p, theta,
hd_pr_p = fd.Get('pr_p_pass_all_exclusivity'+select_phi)
hd_pr_theta = fd.Get('pr_theta_pass_all_exclusivity'+select_phi)

hs_pr_p = fs.Get('pr_p')
hs_pr_theta = fs.Get('pr_theta')

can.cd(1)
hd_pr_p.SetLineColor(kRed)
hs_pr_p.SetLineColor(kBlue)
hs_pr_p.Scale( hd_pr_p.GetMaximum() / hs_pr_p.GetMaximum())
hd_pr_p.SetTitle("Proton Mntm; p (GeV); counts")
hd_pr_p.Draw()
hs_pr_p.Draw("same+hist")
lpr1=TLegend(0.7, 0.6, 0.9, 0.75)
lpr1.AddEntry(hd_pr_p,'DATA')
lpr1.AddEntry(hs_pr_p,'SIM')
lpr1.Draw('same')
can.cd(2)
hd_pr_theta.SetLineColor(kRed)
hs_pr_theta.SetLineColor(kBlue)
hs_pr_theta.Scale( hd_pr_theta.GetMaximum() / hs_pr_theta.GetMaximum())
hd_pr_theta.SetTitle("Proton #theta; #theta (deg); counts")
hd_pr_theta.Draw()
hs_pr_theta.Draw("same+hist")
lpr2=TLegend(0.7, 0.6, 0.9, 0.75)
lpr2.AddEntry(hd_pr_theta,'DATA')
lpr2.AddEntry(hs_pr_theta,'SIM')
lpr2.Draw('same')
can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,2)

######################################### compare Proton kinematics - ptheta, thetaphi
hd_pr_ptheta = fd.Get('pr_ptheta_pass_all_exclusivity'+select_phi)
hd_pr_thetaphi = fd.Get('pr_phitheta_pass_all_exclusivity'+select_phi)

hs_pr_ptheta = fs.Get('pr_ptheta')
hs_pr_thetaphi = fs.Get('pr_thetaphi')

can.cd(1)
hd_pr_ptheta.SetTitle("DATA Pr Mntm vs #theta; p (GeV); #theta (deg)")
hd_pr_ptheta.Draw("colz")
can.cd(2)
hd_pr_thetaphi.SetTitle("DATA Pr #phi vs #theta; #phi (deg); #theta (deg)")
hd_pr_thetaphi.Draw("colz")
can.cd(3)
hs_pr_ptheta.SetTitle("SIM. Pr Mntm vs #theta; p (GeV); #theta (deg)")
hs_pr_ptheta.Draw("colz")
can.cd(4)
hs_pr_thetaphi.SetTitle("SIM. Pr #phi vs #theta; #phi (deg); #theta (deg)")
hs_pr_thetaphi.Draw("colz")
can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)


###########################3# compare kaonP kinematics - p, theta,
hd_kp_p = fd.Get('kp_p_pass_all_exclusivity'+select_phi)
hd_kp_theta = fd.Get('kp_theta_pass_all_exclusivity'+select_phi)

hs_kp_p = fs.Get('kp_p')
hs_kp_theta = fs.Get('kp_theta')

can.cd(1)
hd_kp_p.SetLineColor(kRed)
hs_kp_p.SetLineColor(kBlue)
hs_kp_p.Scale( hd_kp_p.GetMaximum() / hs_kp_p.GetMaximum())
hd_kp_p.SetTitle(" K^{+} Mntm; p (GeV); counts")
hd_kp_p.Draw()
hs_kp_p.Draw("same+hist")
lkp1=TLegend(0.7, 0.6, 0.9, 0.75)
lkp1.AddEntry(hd_kp_p,'DATA')
lkp1.AddEntry(hs_kp_p,'SIM')
lkp1.Draw('same')
can.cd(2)
hd_kp_theta.SetLineColor(kRed)
hs_kp_theta.SetLineColor(kBlue)
hs_kp_theta.Scale( hd_kp_theta.GetMaximum()/ hs_kp_theta.GetMaximum())
hd_kp_theta.SetTitle("K^{+} #theta; #theta (deg); counts")
hd_kp_theta.Draw()
hs_kp_theta.Draw("same+hist")
lkp2=TLegend(0.7, 0.6, 0.9, 0.75)
lkp2.AddEntry(hd_kp_theta,'DATA')
lkp2.AddEntry(hs_kp_theta,'SIM')
lkp2.Draw('same')
can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,2)


######################################### compare kaonN kinematics - p, theta,
hd_kp_ptheta = fd.Get('kp_ptheta_pass_all_exclusivity'+select_phi)
hd_kp_thetaphi = fd.Get('kp_phitheta_pass_all_exclusivity'+select_phi)

hs_kp_ptheta = fs.Get('kp_ptheta')
hs_kp_thetaphi = fs.Get('kp_thetaphi')

can.cd(1)
hd_kp_ptheta.SetTitle("DATA K^{+} Mntm vs #theta; p (GeV); #theta (deg)")
hd_kp_ptheta.Draw("colz")
can.cd(2)
hd_kp_thetaphi.SetTitle("DATA K^{+} #phi vs #theta; #phi (deg); #theta (deg)")
hd_kp_thetaphi.Draw("colz")
can.cd(3)
hs_kp_ptheta.SetTitle("SIM. K^{+} Mntm vs #theta; p (GeV); #theta (deg)")
hs_kp_ptheta.Draw("colz")
can.cd(4)
hs_kp_thetaphi.SetTitle("SIM. K^{+} #phi vs #theta; #phi (deg); #theta (deg)")
hs_kp_thetaphi.Draw("colz")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)


######################################### compare kaonN kinematics - p, theta,
hd_km_p = fd.Get('km_p_pass_all_exclusivity'+select_phi)
hd_km_theta = fd.Get('km_theta_pass_all_exclusivity'+select_phi)

hs_km_p = fs.Get('km_p')
hs_km_theta = fs.Get('km_theta')

can.cd(1)
hd_km_p.SetLineColor(kRed)
hs_km_p.SetLineColor(kBlue)
hs_km_p.Scale( hd_km_p.GetMaximum() / hs_km_p.GetMaximum())
hd_km_p.SetTitle(" K^{-} Mntm; p (GeV); counts")
hd_km_p.Draw()
hs_km_p.Draw("same+hist")
lkm1=TLegend(0.7, 0.6, 0.9, 0.75)
lkm1.AddEntry(hd_km_p,'DATA')
lkm1.AddEntry(hs_km_p,'SIM')
lkm1.Draw('same')
can.cd(2)
hd_km_theta.SetLineColor(kRed)
hs_km_theta.SetLineColor(kBlue)
hs_km_theta.Scale( hd_km_theta.GetMaximum()/hs_km_theta.GetMaximum())
hd_km_theta.SetTitle("K^{-} #theta; #theta (deg); counts")
hd_km_theta.Draw()
hs_km_theta.Draw("same+hist")
lkm2=TLegend(0.7, 0.6, 0.9, 0.75)
lkm2.AddEntry(hd_km_p,'DATA')
lkm2.AddEntry(hs_km_p,'SIM')
lkm2.Draw('same')
can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,2)

######################################### compare kaonN kinematics - p, theta,
hd_km_ptheta = fd.Get('km_ptheta_pass_all_exclusivity'+select_phi)
hd_km_thetaphi = fd.Get('km_phitheta_pass_all_exclusivity'+select_phi)

hs_km_ptheta = fs.Get('km_ptheta')
hs_km_thetaphi = fs.Get('km_thetaphi')

can.cd(1)
hd_km_ptheta.SetTitle("DATA K^{-} Mntm vs #theta; p (GeV); #theta (deg)")
hd_km_ptheta.Draw("colz")
can.cd(2)
hd_km_thetaphi.SetTitle("DATA K^{-} #phi vs #theta; #phi (deg); #theta (deg)")
hd_km_thetaphi.Draw("colz")
can.cd(3)
hs_km_ptheta.SetTitle("SIM. K^{-} Mntm vs #theta; p (GeV); #theta (deg)")
hs_km_ptheta.Draw("colz")
can.cd(4)
hs_km_thetaphi.SetTitle("SIM K^{-} #phi vs #theta; #phi (deg); #theta (deg)")
hs_km_thetaphi.Draw("colz")
can.Print(mon_out_file_name)
can.Clear()
can.Divide(3,2)

##############################################################################
## compare missing mass between reconstructed data and reconstructed simulation

hd_eX = fd.Get('LeX_pass_all_exclusivity'+select_phi)
hd_epX = fd.Get('epX_pass_all_exclusivity'+select_phi)
hd_epkpX = fd.Get('epkpX_pass_all_exclusivity'+select_phi)
hd_epkmX = fd.Get('epkmX_pass_all_exclusivity'+select_phi)
hd_ekpkmX = fd.Get('ekpkmX_pass_all_exclusivity'+select_phi)
hd_epkpkmX = fd.Get('epkpkmX_pass_all_exclusivity'+select_phi)
hd_pt = fd.Get('pt_pass_all_exclusivity'+select_phi)

hs_eX = fs.Get('LeX')
hs_epX = fs.Get('epX')
hs_epkpX = fs.Get('epkpX')
hs_epkmX = fs.Get('epkmX')
hs_ekpkmX = fs.Get('ekpkmX')
hs_epkpkmX = fs.Get('epkpkmX')
hs_pt = fs.Get('pt')

can.cd(1)
hd_eX.SetLineColor(kRed)
hs_eX.SetLineColor(kBlue)
hd_eX.Scale(hs_eX.GetMaximum() / hd_eX.GetMaximum())
hd_eX.SetTitle("ep #rightarrow eX MM2;MM2 (GeV^2);counts")
hd_eX.Draw("hist")
hs_eX.Draw("same+hist")
can.cd(2)
hd_epX.SetLineColor(kRed)
hs_epX.SetLineColor(kBlue)
hd_epX.Scale(hs_epX.GetMaximum() / hd_epX.GetMaximum())
hd_epX.SetTitle("ep #rightarrow epX MM2;MM2 (GeV^2);counts")
hd_epX.Draw("hist")
hs_epX.Draw("same+hist")
can.cd(3)
hd_epkpX.SetLineColor(kRed)
hs_epkpX.SetLineColor(kBlue)
hd_epkpX.Scale(hs_epkpX.GetMaximum() / hd_epkpX.GetMaximum())
print('number of data entries %d, number of sim enries %d ' % (hd_epkpX.GetMaximum(), hs_epkpX.GetMaximum()))
hd_epkpX.SetTitle("ep #rightarrow epK^{+}X MM2;MM2 (GeV^2);counts")
hs_epkpX.Draw("hist")
hd_epkpX.Draw("same+hist")
can.cd(4)
hd_epkmX.SetLineColor(kRed)
hs_epkmX.SetLineColor(kBlue)
hd_epkmX.Scale(hs_epkmX.GetMaximum() / hd_epkmX.GetMaximum())
hd_epkmX.SetTitle("ep #rightarrow epK^{-}X MM2;MM2 (GeV^2);counts")
hd_epkmX.Draw("hist")
hs_epkmX.Draw("same+hist")
can.cd(5)
hd_ekpkmX.SetLineColor(kRed)
hs_ekpkmX.SetLineColor(kBlue)
print('number of data entries %d, number of sim enries %d ' % (hd_ekpkmX.GetMaximum(), hs_ekpkmX.GetMaximum()))
hd_ekpkmX.Scale(hs_ekpkmX.GetMaximum() / hd_ekpkmX.GetMaximum())
hd_ekpkmX.SetTitle("ep #rightarrow eK^{+}K^{-}X MM2;MM2 (GeV^2);counts")
hd_ekpkmX.Draw("hist")
hs_ekpkmX.Draw("same+hist")
can.cd(6)
hd_epkpkmX.SetLineColor(kRed)
hs_epkpkmX.SetLineColor(kBlue)
hd_epkpkmX.Scale(hs_epkpkmX.GetMaximum() / hd_epkpkmX.GetMaximum())
hd_epkpkmX.SetTitle("ep #rightarrow epK^{+}K^{-}X MM2;MM2 (GeV^2);counts")
hd_epkpkmX.Draw("hist")
hs_epkpkmX.Draw("same+hist")


can.Print(mon_out_file_name)
can.Clear()
can.SetCanvasSize(600,600)
can.Divide(1,1)


hd_q2x = fd.Get('xq2_pass_all_exclusivity'+select_phi)
hd_txb = fd.Get('txb_pass_all_exclusivity'+select_phi)
hd_q2t = fd.Get('q2t_pass_all_exclusivity'+select_phi)

hs_q2x = fs.Get('xq2')
hs_txb = fs.Get('txb')
hs_q2t = fs.Get('q2t')

hs_q = fs.Get('q2')
hs_t = fs.Get('t')
can.cd(1)
hd_q=hd_q2x.ProjectionY()
hd_q.SetLineColor(kRed)
hs_q.SetLineColor(kBlue)
hd_q.Scale(hs_q.GetMaximum()/hd_q.GetMaximum())
hd_q.SetTitle('Q2 Data (red) vs Sim (blue);Q^{2} (GeV^2);counts')
hd_q.Draw('hist')
hs_q.Draw('same+hist')
lq1=TLegend(0.7, 0.6, 0.9, 0.75)
lq1.AddEntry(hd_q,'DATA')
lq1.AddEntry(hs_q,'SIM')
lq1.Draw('same')
can.Print(mon_out_file_name)


can.Clear()
can.SetCanvasSize(600,600)
can.Divide(1,1)
can.cd(1)
hd_t=hd_txb.ProjectionX()
hd_t.SetLineColor(kRed)
hs_t.SetLineColor(kBlue)
hd_t.Scale(hs_t.GetMaximum()/hd_t.GetMaximum())
hd_t.SetTitle('-t Data (red) vs Sim (blue);-t (GeV^2);counts')
hd_t.Draw('hist')
hs_t.Draw('same+hist')
lt1=TLegend(0.7, 0.6, 0.9, 0.75)
lt1.AddEntry(hd_t,'DATA')
lt1.AddEntry(hs_t,'SIM')
lt1.Draw('same')

can.Print(mon_out_file_name)

can.Clear()
can.SetCanvasSize(1800,1200)
can.Divide(3,2)
can.cd(1)
hd_q2x.Rebin2D(2,2)
hd_q2x.SetTitle('DATA: Q2 vs xb ; xb; Q2 (GeV^2)')
hd_q2x.Draw('colz')
can.cd(2)
hd_txb.Rebin2D(2,2)
hd_txb.SetTitle('DATA: xb vs -t ; -t (GeV^2); xb')
hd_txb.Draw('colz')
can.cd(3)
hd_q2t.SetTitle('DATA: -t vs Q2 ; Q2 (GeV^2); -t (GeV^2)')
hd_q2t.Rebin2D(6,6)
hd_q2t.Draw('colz')
can.cd(4)
hs_q2x.Rebin2D(2,2)
hs_q2x.SetTitle('SIM: Q2 vs xb ; xb; Q2 (GeV^2)')
hs_q2x.Draw('colz')
can.cd(5)
hs_txb.Rebin2D(2,2)
hs_txb.SetTitle('SIM: xb vs -t ; -t (GeV^2); xb')
hs_txb.Draw('colz')
can.cd(6)
hs_q2t.Rebin2D(2,2)
hs_q2t.SetTitle('SIM: -t vs Q2 ; Q2 (GeV^2); -t (GeV^2)')
hs_q2t.Draw('colz')


can.Print(mon_out_file_name)

hd_prb = fd.Get('pro_beta_vs_p_raw')
hd_kpb = fd.Get('kp_beta_vs_p_raw')
hd_kmb = fd.Get('km_beta_vs_p_raw')
hs_prb = fs.Get('pro_beta_vs_p_raw')
hs_kpb = fs.Get('kp_beta_vs_p_raw')
hs_kmb = fs.Get('km_beta_vs_p_raw')

can.Clear()
can.SetCanvasSize(1800,1200)
can.Divide(3,2)
can.cd(1)
hd_prb.SetTitle('DATA proton #beta vs p; p (GeV); #beta')
hd_prb.Draw('colz')
can.cd(2)
hd_kpb.SetTitle('DATA K^{+} #beta vs p; p (GeV); #beta')
hd_kpb.Draw('colz')
can.cd(3)
hd_kmb.SetTitle('DATA K^{-} #beta vs p; p (GeV); #beta')
hd_kmb.Draw('colz')
can.cd(4)
hs_prb.SetTitle('SIM proton #beta vs p; p (GeV); #beta')
hs_prb.Draw('colz')
can.cd(5)
hs_kpb.SetTitle('SIM K^{+} #beta vs p; p (GeV); #beta')
hs_kpb.Draw('colz')
can.cd(6)
hs_kmb.SetTitle('SIM K^{-} #beta vs p; p (GeV); #beta')
hs_kmb.Draw('colz')


can.Print(mon_out_file_name)



'''
can.cd(1)
hd_pt.SetLineColor(kRed)
hs_pt.SetLineColor(kBlue)
hd_pt.Scale(hs_pt.GetMaximum() / hd_pt.GetMaximum())
hd_pt.SetTitle("Missing Transverse Mntm;p_{T} (GeV);counts")
hd_pt.Draw("hist")
hs_pt.Draw("same+hist")

can.Print(mon_out_file_name)
can.Clear()
can.SetCanvasSize(600,600)
can.Divide(1,1)

##############################################################################
## compare phi meson mass between data and simulation

hd_vphi = fd.Get('vphimass_pass_all_exclusivity')
hs_vphi = fs.Get('vphimass')

can.cd(1)
hd_vphi.SetLineColor(kRed)
hs_vphi.SetLineColor(kBlue)
hd_vphi.RebinX(3)
hd_vphi.SetTitle("Invariant Masss of K^{+}K^{-} from ep #rightarrow ep#phi;mass (GeV);counts")
hd_vphi.Scale(hs_vphi.GetMaximum() / hd_vphi.GetMaximum() )
hd_vphi.Draw("hist")
hs_vphi.Draw("same+hist")
can.Print(mon_out_file_name)
can.Clear()
can.Divide(1,1)

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)
'''

can.cd


can.Print('{}]'.format(mon_out_file_name))


##############################################################################
## compare cos_kp
'''
hd_coskp = fd.Get('cos_kp_evcut')
hs_coskp = fs.Get('cos_kp_evcut')

ccoskp = TCanvas("ccoskp","ccoskp",900,900)
ccoskp.cd(1)
hd_coskp.SetLineColor(kRed)
hs_coskp.SetLineColor(kBlue)
hd_coskp.SetTitle("cos(#theta_{K^{+}});cos(#theta_{K^{+}});counts")
hd_coskp.Scale(hs_coskp.GetMaximum() / hd_coskp.GetMaximum() )
hd_coskp.Draw("hist")
hs_coskp.Draw("same+hist")
ccoskp.SaveAs("pics/comp_coskp_"+base+".pdf")
'''















    
    
    
    
