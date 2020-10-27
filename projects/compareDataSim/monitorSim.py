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
fs = TFile(sys.argv[1])
datatype = sys.argv[2]
tag=sys.argv[3]
base=os.path.basename(fd.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_mc_events_'+tag+'.pdf'
can = TCanvas('can','can',1200,3200)
can.Print('{}['.format(mon_out_file_name))
can.Divide(3,8)
####################################################
####################################################
 # code to compare reconstructed MC to generated MC
####################################################
####################################################



#print out all p, theta, phi, ptheta,thetaphi, pphi 
# for each
l_el_kin = [fs.Get('el_p'),fs.Get('el_theta'),fs.Get('el_phi'), fs.Get('el_ptheta'),fs.Get('el_thetaphi'),fs.Get('el_pphi') ]
l_pr_kin = [fs.Get('pr_p'),fs.Get('pr_theta'),fs.Get('pr_phi'), fs.Get('pr_ptheta'),fs.Get('pr_thetaphi'),fs.Get('pr_pphi') ]
l_kp_kin = [fs.Get('kp_p'),fs.Get('kp_theta'),fs.Get('kp_phi'), fs.Get('kp_ptheta'),fs.Get('kp_thetaphi'),fs.Get('kp_pphi') ]
l_km_kin = [fs.Get('km_p'),fs.Get('km_theta'),fs.Get('km_phi'), fs.Get('km_ptheta'),fs.Get('km_thetaphi'),fs.Get('km_pphi') ]
l_part_kin = [l_el_kin, l_pr_kin, l_kp_kin, l_km_kin]
can_counter=1
for pp in l_part_kin:
    for hh in pp:
        can.cd(can_counter)
        can_counter+=1
        #if( hh.InheritsFrom('TH1') ) :
        #    hh.Draw()
        #if( hh.InheritsFrom('TH2') ):
        hh.Draw('colz')
        
can.Print(mon_out_file_name)
can.Clear()
can.SetCanvasSize(1200,600)
can.Divide(2,1)

# compare electron kinematics - p, theta,
hd_el_p = fs.Get('el_p')
hd_el_theta = fs.Get('el_theta')

hs_el_p = fs.Get('gen_el_p')
hs_el_theta = fs.Get('gen_el_theta')

#celp = TCanvas('celp','celp',900,900)
can.cd(1)
hd_el_p.SetLineColor(kRed)
hs_el_p.SetLineColor(kBlue)
#hs_el_p.Scale( hd_el_p.GetMaximum() / hs_el_p.GetMaximum())
hd_el_p.SetTitle("Electron Mntm; p (GeV); counts")
hd_el_p.Draw()
hs_el_p.Draw("same+hist")
leg=TLegend(0.7,0.7,0.9,0.9)
leg.AddEntry(hd_el_p,'data')
leg.AddEntry(hs_el_p,'sim')

can.cd(2)
hd_el_theta.SetLineColor(kRed)
hs_el_theta.SetLineColor(kBlue)
#hs_el_theta.Scale( hd_el_theta.GetMaximum() / hs_el_theta.GetMaximum())
hd_el_theta.SetTitle("Electron #theta; #theta (deg); counts")
hd_el_theta.Draw()
hs_el_theta.Draw("same+hist")
leg=TLegend(0.7,0.7,0.9,0.9)
leg.AddEntry(hd_el_theta,'data')
leg.AddEntry(hs_el_theta,'sim')

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)
########################## compare protons kinematics - p, theta,
hd_pr_p = fd.Get('pr_p')
hd_pr_theta = fd.Get('pr_theta')

hs_pr_p = fs.Get('gen_pr_p')
hs_pr_theta = fs.Get('gen_pr_theta')

can.cd(1)
hd_pr_p.SetLineColor(kRed)
hs_pr_p.SetLineColor(kBlue)
#hs_pr_p.Scale( hs_pr_p.GetMaximum() / hd_pr_p.GetMaximum())
hd_pr_p.SetTitle("Proton Mntm; p (GeV); counts")
hd_pr_p.Draw()
hs_pr_p.Draw("same+hist")
can.cd(2)
hd_pr_theta.SetLineColor(kRed)
hs_pr_theta.SetLineColor(kBlue)
#hs_pr_theta.Scale( hs_pr_theta.GetMaximum() / hd_pr_theta.GetMaximum())
hd_pr_theta.SetTitle("Proton #theta; #theta (deg); counts")
hd_pr_theta.Draw()
hs_pr_theta.Draw("same+hist")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)

###########################3# compare kaonP kinematics - p, theta,
hd_kp_p = fd.Get('kp_p')
hd_kp_theta = fd.Get('kp_theta')

hs_kp_p = fs.Get('gen_kp_p')
hs_kp_theta = fs.Get('gen_kp_theta')

can.cd(1)
hd_kp_p.SetLineColor(kRed)
hs_kp_p.SetLineColor(kBlue)
#hs_kp_p.Scale( hd_kp_p.GetMaximum() / hs_kp_p.GetMaximum())
hd_kp_p.SetTitle(" K^{+} Mntm; p (GeV); counts")
hd_kp_p.Draw()
hs_kp_p.Draw("same+hist")
can.cd(2)
hd_kp_theta.SetLineColor(kRed)
hs_kp_theta.SetLineColor(kBlue)
#hs_kp_theta.Scale( hd_kp_theta.GetMaximum()/ hs_kp_theta.GetMaximum())
hd_kp_theta.SetTitle("K^{+} #theta; #theta (deg); counts")
hd_kp_theta.Draw()
hs_kp_theta.Draw("same+hist")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(2,1)

######################################### compare kaonN kinematics - p, theta,
hd_km_p = fd.Get('km_p')
hd_km_theta = fd.Get('km_theta')

hs_km_p = fs.Get('gen_km_p')
hs_km_theta = fs.Get('gen_km_theta')

can.cd(1)
hd_km_p.SetLineColor(kRed)
hs_km_p.SetLineColor(kBlue)
#hs_km_p.Scale( hd_kp_p.GetMaximum() / hs_kp_p.GetMaximum())
hd_km_p.SetTitle(" K^{-} Mntm; p (GeV); counts")
hd_km_p.Draw()
hs_km_p.Draw("same+hist")
can.cd(2)
hd_km_theta.SetLineColor(kRed)
hs_km_theta.SetLineColor(kBlue)
#hs_km_theta.Scale( hd_km_theta.GetMaximum()/ hs_km_theta.GetMaximum())
hd_km_theta.SetTitle("K^{-} #theta; #theta (deg); counts")
hd_km_theta.Draw()
hs_km_theta.Draw("same+hist")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(3,2)

##############################################################################
## compare missing mass between reconstructed  and generated

hd_eX = fd.Get('LeX')
hd_epX = fd.Get('epX')
hd_epkpX = fd.Get('epkpX')
hd_epkmX = fd.Get('epkmX')
hd_ekpkmX = fd.Get('ekpkmX')
hd_epkpkmX = fd.Get('epkpkmXe')
hd_pt = fd.Get('pt')

hs_eX = fs.Get('gen_LeX')
hs_epX = fs.Get('gen_epX')
hs_epkpX = fs.Get('gen_epkpX')
hs_epkmX = fs.Get('gen_epkmX')
hs_ekpkmX = fs.Get('gen_ekpkmX')
hs_epkpkmX = fs.Get('gen_epkpkmXe')
hs_pt = fs.Get('gen_pt')

can.cd(1)
hd_eX.SetLineColor(kRed)
hs_eX.SetLineColor(kBlue)
hd_eX.Scale(hs_eX.GetMaximum() / hd_eX.GetMaximum())
hd_eX.SetTitle("ep #rightarrow eX MM2;MM2 (GeV^2);counts")
hd_eX.SetAxisRange(3, 11, "X")
hd_eX.Draw("hist")
hs_eX.Draw("same+hist")

can.cd(2)
hd_epX.SetLineColor(kRed)
hs_epX.SetLineColor(kBlue)
hd_epX.Scale(hs_epX.GetMaximum() / hd_epX.GetMaximum())
hd_epX.SetTitle("ep #rightarrow epX MM2;MM2 (GeV^2);counts")
hd_epX.SetAxisRange(0.6, 1.5, "X")
hd_epX.Draw("hist")
hs_epX.Draw("same+hist")

can.cd(3)
hd_epkpX.SetLineColor(kRed)
hs_epkpX.SetLineColor(kBlue)
hd_epkpX.Scale(hs_epkpX.GetMaximum() / hd_epkpX.GetMaximum())
hd_epkpX.SetTitle("ep #rightarrow epK^{+}X MM2;MM2 (GeV^2);counts")
hd_epkpX.SetAxisRange(0.0, 0.5, "X")
hd_epkpX.Draw("hist")
hs_epkpX.Draw("same+hist")

can.cd(4)
hd_epkmX.SetLineColor(kRed)
hs_epkmX.SetLineColor(kBlue)
hd_epkmX.Scale(hs_epkmX.GetMaximum() / hd_epkmX.GetMaximum())
hd_epkmX.SetTitle("ep #rightarrow epK^{-}X MM2;MM2 (GeV^2);counts")
hd_epkmX.SetAxisRange(0.0, 0.5, "X")
hd_epkmX.Draw("hist")
hs_epkmX.Draw("same+hist")

can.cd(5)
hd_ekpkmX.SetLineColor(kRed)
hs_ekpkmX.SetLineColor(kBlue)
#print('number of data entries %d, number of sim enries %d ' % (hd_ekpkmX.GetMaximum(), hs_ekpkmX.GetMaximum()))
hd_ekpkmX.Scale(hs_ekpkmX.GetMaximum() / hd_ekpkmX.GetMaximum())
hd_ekpkmX.SetTitle("ep #rightarrow eK^{+}K^{-}X MM2;MM2 (GeV^2);counts")
hd_ekpkmX.SetAxisRange(0.6, 1.2, "X")
hd_ekpkmX.Draw("hist")
hs_ekpkmX.Draw("same+hist")

can.cd(6)
hd_epkpkmX.SetLineColor(kRed)
hs_epkpkmX.SetLineColor(kBlue)
hd_epkpkmX.Scale(hs_epkpkmX.GetMaximum() / hd_epkpkmX.GetMaximum())
hd_epkpkmX.SetTitle("ep #rightarrow epK^{+}K^{-}X ; Missing Energy (GeV);counts")
hd_epkpkmX.SetAxisRange(-0.2, 0.2, "X")
hd_epkpkmX.Draw("hist")
hs_epkpkmX.Draw("same+hist")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(1,1)
can.SetCanvasSize(600,600)

#can.cd(1)
#hd_pt.SetLineColor(kRed)
#hs_pt.SetLineColor(kBlue)
#hd_pt.Scale(hs_pt.GetMaximum() / hd_pt.GetMaximum())
#hd_pt.SetTitle("Missing Transverse Mntm;p_{T} (GeV);counts")
#hd_pt.Draw("hist")
#hs_pt.Draw("same+hist")
#can.Print(mon_out_file_name)
#can.Clear()
#can.Divide(1,1)
#can.SetCanvasSize(600,600)

##############################################################################
## compare phi meson mass between data and simulation
'''hd_vphi = fd.Get('vphimass_evcut')
hs_vphi = fs.Get('vphimass')

can.cd(1)
hd_vphi.SetLineColor(kRed)
hs_vphi.SetLineColor(kBlue)
hd_vphi.SetTitle("Invariant Masss of K^{+}K^{-} from ep #rightarrow ep#phi;mass (GeV);counts")
#hd_vphi.Scale(hs_vphi.GetMaximum() / hd_vphi.GetMaximum() )
hd_vphi.Draw("hist")
hs_vphi.Draw("same+hist")

can.Print(mon_out_file_name)
can.Clear()
can.Divide(1,1)
can.SetCanvasSize(600,600)

##############################################################################
## compare cos_kp
hd_coskp = fd.Get('cos_kp_evcut')
hs_coskp = fs.Get('cos_kp_evcut')

can.cd(1)
hd_coskp.SetLineColor(kRed)
hs_coskp.SetLineColor(kBlue)
hd_coskp.SetTitle("cos(#theta_{K^{+}});cos(#theta_{K^{+}});counts")
#hd_coskp.Scale(hs_coskp.GetMaximum() / hd_coskp.GetMaximum() )
hd_coskp.Draw("hist")
hs_coskp.Draw("same+hist")
'''
can.Print('{}]'.format(mon_out_file_name))
