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
datatype = sys.argv[2]
tag=sys.argv[3]
base=os.path.basename(fd.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_philamba_events_'+tag+'.pdf'
can = TCanvas('can','can',1200,1200)
can.Print('{}['.format(mon_out_file_name))
#can.Divide(3,8)
####################################################
####################################################
 # code to compare reconstructed MC to generated MC
####################################################
####################################################

str_pass_all = 'pass_all_exclusivity'
str_remove_lambda1 = 'pass_all_exclusivity_phi_rmv_lambda'
str_remove_lambda2 = 'pass_all_exclusivity_rmv_both_lambdas'

hd_lambda_pass_all = fd.Get('lambdamass_'+str_pass_all)
hd_philambda_pass_all = fd.Get('lambdaphimass_'+str_pass_all)
hd_phi_pass_all = fd.Get('vphimass_'+str_pass_all)
hd_phi_pass_all_rmv_lambda1 = fd.Get('vphimass_'+str_remove_lambda1)
hd_phi_pass_all_rmv_lambda12 = fd.Get('vphimass_'+str_remove_lambda2)

#draw dalitz plot first
#can.Clear()
can.Divide(1,1)
hd_philambda_pass_all.SetTitle("Phi vs Lambda - 2D; inv. mass p K^{-} (GeV); inv. mass K^{+}K^{-}")
hd_philambda_pass_all.Draw("colz")
can.Print(mon_out_file_name)

#show 1D plots of lambda (prKm)
can.Clear()
can.Divide(1,1)
hd_lambda_pass_all.SetTitle("Inv. Mass Pr K^{-}; inv. mass Pr K^{-}; counts")
hd_lambda_pass_all.Draw()
can.Print(mon_out_file_name)

#show phi after remove lambdas
can.Clear()
can.Divide(1,1)
hd_phi_pass_all.SetTitle("Inv. Mass K^{+}K^{-}; inv. mass K^{+}K^{-}; counts")
hd_phi_pass_all.Draw()
hd_phi_pass_all_rmv_lambda1.SetLineColor(kRed)
hd_phi_pass_all_rmv_lambda1.Draw("same")
hd_phi_pass_all_rmv_lambda12.SetLineColor(kBlue)
hd_phi_pass_all_rmv_lambda12.Draw("same")
can.Print(mon_out_file_name)


can.Print('{}]'.format(mon_out_file_name))

