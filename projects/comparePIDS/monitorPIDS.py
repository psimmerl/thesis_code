from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue
import numpy as np

from array import array
import math
import sys,os
'''
fd = [TFile(sys.argv[1]),TFile(sys.argv[2]),TFile(sys.argv[3])]
datatype = sys.argv[4]
tag=sys.argv[5]
base=os.path.basename(fd[0].GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_philamba_events_pidmethods'+tag+'.pdf'
can = TCanvas('can','can',1200,1200)
can.Print('{}['.format(mon_out_file_name))

#fd[0] - EB only
#fd[1] - all cuts but new cut
#fd[2] - all cuts with new cut


str_pass_all = 'pass_all_exclusivity'
str_remove_lambda1 = 'pass_all_exclusivity_phi_rmv_lambda'
str_remove_lambda2 = 'pass_all_exclusivity_rmv_both_lambdas'
 

can.Divide(1,1)
can.cd(1)
hd_phi_pass_all1 = fd[0].Get('vphimass_'+str_pass_all)
hd_phi_pass_all2 = fd[1].Get('vphimass_'+str_pass_all)
hd_phi_pass_all3 = fd[2].Get('vphimass_'+str_pass_all)
hd_phi_pass_all1.SetTitle("Inv. Mass K^{+}K^{-}; inv. mass K^{+}K^{-}; counts")

hd_phi_pass_all1.Draw()
hd_phi_pass_all2.SetLineColor(kRed)
hd_phi_pass_all3.SetLineColor(kBlue)
hd_phi_pass_all2.Draw("same")
hd_phi_pass_all3.Draw("same")

can.Print(mon_out_file_name)

can.Clear()
can.Divide(1,1)
can.cd(1)

hd_lambda_pass_all1 = fd[0].Get('lambdamass_'+str_pass_all)
hd_lambda_pass_all2 = fd[1].Get('lambdamass_'+str_pass_all)
hd_lambda_pass_all3 = fd[2].Get('lambdamass_'+str_pass_all)
hd_lambda_pass_all1.SetTitle("Inv. Mass Pr K^{-}; inv. mass Pr K^{-}; counts")

hd_lambda_pass_all1.Draw()
hd_lambda_pass_all2.SetLineColor(kRed)
hd_lambda_pass_all3.SetLineColor(kBlue)
hd_lambda_pass_all2.Draw("same")
hd_lambda_pass_all3.Draw("same")

can.Print(mon_out_file_name)
'''


can.Print('{}]'.format(mon_out_file_name))
