from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLatex, TLine
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kBlack, kViolet, kTeal

from ROOT import Fit
from ROOT import Math

import numpy as np
from collections import OrderedDict 

from array import array
import math
import sys,os

mon_out_file_name='monitor_elastic_lumi_rates_per_sector.pdf'
can = TCanvas('can','can',1200,1200)
can.Print('{}['.format(mon_out_file_name))
can.Divide(2,1)

lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)

gStyle.SetOptStat(00000)

dump=[]
map_sector_rates = {}

runs = [5418, 5419, 5306, 5443, 5444, 5453]
### 5443 , 5444, 5306
special_runs = [5443, 5444, 5306]
special_run_current=[5, 20, 45]
special_rates=[]



my_colors=[kGreen, kRed, kViolet, kBlue, kTeal, kBlack]
for rr in runs:
    sect_rates = []
    fin=open('mon_incl_elastic_r{}_per_sector.txt'.format(rr))
    for ll in fin:
        sect_rates.append( float( ll.split(' ')[1] ) )

    
    map_sector_rates[rr] = sect_rates

#########################
## special bs
for rr in special_runs:
    sect_rates = []
    fin=open('mon_incl_elastic_r{}_per_sector.txt'.format(rr))
    for ll in fin:
        sect_rates.append( float( ll.split(' ')[1] ) )
        
    special_rates.append(sect_rates[0])


print special_rates


print map_sector_rates

legend = TLegend(0.1,0.7,0.23,0.9)
legend.SetHeader("Legend","C") #option "C" allows to center the header
mg_rates_per_sector = TMultiGraph()
jj=0
for key, items in map_sector_rates.items():
    
    x_val = array('d')
    y_val = array('d')
    x_err = array('d')
    y_err = array('d')
    ss=1

    sector1_rates = 0
    for ii in items:
        x_val.append(ss)
        y_val.append(ii)
        x_err.append(0)
        y_err.append(0)
        if ss == 1 :
            sector1_rates = ii

        ss+=1
    g_rates = TGraphErrors(len(x_val), x_val, y_val, x_err, y_err)
    g_rates.SetTitle('R{}'.format(key) )
    g_rates.SetMarkerStyle(21)
    g_rates.SetMarkerColor(my_colors[jj])
    dump.append(g_rates)
 
    legend.AddEntry(g_rates)
    mg_rates_per_sector.Add(g_rates)
    dump.append(mg_rates_per_sector)
    

    jj+=1


print od

mg_rates_one_sector = TMultiGraph()

x_run_val = array('d')
y_run_val = array('d')
x_zer = array('d')
y_zer = array('d')
for kk in range(0, len(special_rates)):        
    x_run_val.append(special_run_current[kk])
    y_run_val.append(special_rates[kk])
    x_zer.append(0.0)
    y_zer.append(0.0)


gPad.SetLogy()                 
g_run_rates = TGraphErrors(len(x_run_val), x_run_val, y_run_val, x_zer, y_zer)
g_run_rates.SetTitle("")
g_run_rates.SetMarkerStyle(21)
g_run_rates.SetMarkerColor(my_colors[0])
g_run_rates.Draw('AP')
lab.DrawLatex(0.1, 0.925, 'Sector 1 Signal to FC vs Run')
lab.DrawLatex(0.45, 0.015, 'Run')
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, 'Signal / Faraday Cup')
lab.SetTextAngle(0)
can.Print(mon_out_file_name)

mg_rates_per_sector.SetTitle("")
title="Signal/Faraday Cup Per Sector for ep->eX"
xtitle="Sector"
ytitle="Signal / Faraday Cup"
mg_rates_per_sector.Draw("AP")
lab.DrawLatex(0.1, 0.925, title)
lab.DrawLatex(0.45, 0.015, xtitle)
lab.SetTextAngle(90)
lab.DrawLatex(0.04, 0.5, ytitle)
lab.SetTextAngle(0)

legend.Draw("same");

can.Print(mon_out_file_name)

can.Print('{}]'.format(mon_out_file_name))



