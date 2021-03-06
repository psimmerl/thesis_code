#!/usr/bin/env python 
import argparse 
import numpy as np
import math
# Trick for docker install to run headless
# and never use plt.show() 
# dont use matplotlib while running this code on the batch farm
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt 

################################################################
## note - this is used for analyzing the RAD elastic events
## that includes both ISR and FSR. First half is ISR CTOF

from array import array 
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors, TMultiGraph, TLegend,
                  TLine, kRed, kBlack, kBlue, kWhite, kGreen, kViolet, kPink, kAzure, kOrange, 
                  kTRUE, gROOT)

from ROOT import( RooRealVar, RooArgSet, RooArgList, RooDataHist,
                  RooGaussian, RooFit, RooAddPdf, RooPolynomial, RooFormulaVar, RooLinkedList)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
global_dump = []

gROOT.SetBatch(kTRUE);

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def setup_global_options():
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gStyle.SetPalette(57)



def add_text_page(can, label, text, save_name):
    """ Write some text onto a page. """
    can.Clear()
    can.cd(1)
    label.DrawLatex(0.1, 0.5, text)
    can.Print(save_name)

def plot_asy_methods( can, file_in, file_names, save_name, label, title, xtitle, ytitle, bin_infos=None, field=None ):
    dumpster=[]

    # txt file with bin center and asymetrry value with errors
    all_phi_trento_results=[]
    phi_asy_meas_result = []
    cc=0
    for nn in file_names:
        fin = open(file_in.format(nn),'r')

        x_cnt = array('d')
        y_val = array('d')
        x_err = array('d')
        y_err = array('d')

        ## load parameters for method 1 (asy vs kpkm mass )
        for ll in fin:
            rr = ll.split(' ')
            print rr
            bc = float(rr[0])
            asy = float(rr[1])
            asy_err = float(rr[2])
            x_cnt.append(bc)
            y_val.append(asy)
            x_err.append(0)
            y_err.append(asy_err)


        # asy graph for each name instane
        asy_graph = TGraphErrors(len(x_cnt), x_cnt, y_val, x_err, y_err)
        asy_graph.SetTitle('')
        asy_graph.SetMarkerStyle(22)
        asy_graph.SetMarkerSize(2)
        asy_graph.SetLineColor(bin_infos['color'][cc])
        asy_graph.SetMarkerColor(bin_infos['color'][cc])
        dumpster.append(asy_graph)

        ## fit the trento dist for
        fit_sin = TF1("fit_sin_{}".format(cc),"[0]*sin(x * TMath::Pi()/180.0)",-180,180)
        fit_sin.SetParameter(0,1)
        fit_sin.SetLineColor(bin_infos['color'][cc])
        fit_sin.SetLineStyle(2)
        asy_graph.Fit("fit_sin_{}".format(cc),"R")
        dumpster.append(fit_sin)

        phi_asy_meas_result.append(fit_sin)
        all_phi_trento_results.append(asy_graph)
        cc+=1

    mg_asy_comp = TMultiGraph()
    ## legend
    leg = TLegend(0.7, 0.7, 0.9, 0.9)
    cc=0
    for gg in all_phi_trento_results:
        ## draw both graphs togethers to compare m1 and m2
        mg_asy_comp.Add(gg)
        leg.AddEntry(gg,bin_infos['name'][cc],'p')
        cc+=1
    dumpster.append(mg_asy_comp)

    mg_asy_comp.GetHistogram().SetMaximum(0.5)
    mg_asy_comp.GetHistogram().SetMinimum(-0.5)
    mg_asy_comp.GetXaxis().SetLimits(-180.0, 180.0)
    dumpster.append(mg_asy_comp)
    
    # draw grid
    h = gPad.DrawFrame(1.4, -0.5, 0.95, 0.5)
    h.GetXaxis().SetNdivisions(-9)
    h.GetYaxis().SetNdivisions(-9)
    gPad.SetGrid()
    dumpster.append(h)

    mg_asy_comp.Draw('AP')    
    dumpster.append(mg_asy_comp)

    leg.Draw('same')
    dumpster.append(leg)

    ## draw asy fit information between two methods
    cc=0
    for asy_fit in phi_asy_meas_result:
        label.DrawLatex(0.13, 0.85-0.062*cc,'#color[{}]'.format(bin_infos['color'][cc]) + '{'+'Asy^{m}_{'+bin_infos['name'][cc]+'}'+'{0:.4f} #pm {1:.4f}'.format(asy_fit.GetParameter(0), asy_fit.GetParError(0)) +'}')        
        cc+=1
    if title:
        label.DrawLatex(0.1, 0.925, title)        
    if xtitle:
        label.DrawLatex(0.65, 0.025, xtitle)            
    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.035, 0.5, ytitle)
        label.SetTextAngle(0)
    dumpster.append(label)

    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.13, 0.125, field)#125, 0.847, field)
        dumpster.append(flabel)
    

    can.Print(output_pdfname)
    

    
    


    
if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    output_pdfname = args.output_prefix + '.pdf'
    setup_global_options() 

    can = TCanvas('can', 'can', 1000, 600)
    can.SetBatch(kTRUE)


    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    can.Clear()
    can.Update()

    # start looking at total integrated asy
    #asy_per_mass = open('asy_vs_imkpkm_pidtype1_additionalcuts_binVar3.txt','r')
    
    
    bin_info = {'min_sig_bin':0,
                'max_sig_bin':2,
                'color': [880-5, 632, 860-6, 416-1],
                'name':['#pm 1.5#sigma','#pm 2#sigma', '#pm 2.5#sigma', '#pm 3#sigma'],
                
            }
    field_setting = 'inbNoutb'
    

    fin_name='asy_vs_trento_method1_imkpkm_bin1_pidtype1_additionalcuts_rmvres12_binVar11{}_'+field_setting+'.txt'
    plot_asy_methods( can, fin_name, ['a','b','c','d'], output_pdfname, lab, 'BSA Field Comparison at Different #sigma - Left Sideband',
                      '#phi_{trento} (deg)', 'Asy', bin_infos=bin_info, field=field_setting)


    fin_name='asy_vs_trento_method1_imkpkm_bin2_pidtype1_additionalcuts_rmvres12_binVar11{}_'+field_setting+'.txt'
    plot_asy_methods( can, fin_name, ['a','b','c','d'], output_pdfname, lab, 'BSA Field Comparison at Different #sigma - Peak', 
                      '#phi_{trento} (deg)', 'Asy', bin_infos=bin_info, field=field_setting)

    fin_name='asy_vs_trento_method1_imkpkm_bin3_pidtype1_additionalcuts_rmvres12_binVar11{}_'+field_setting+'.txt'
    plot_asy_methods( can, fin_name, ['a','b','c','d'], output_pdfname, lab, 'BSA Field Comparison at Different #sigma - Right Sideband', 
                      '#phi_{trento} (deg)', 'Asy', bin_infos=bin_info, field=field_setting)


    '''
    # compare asy vs -t, but integrated over Q2, xb
    for tt in range(1,len(t_bins)):
        #plot_asy( can, asy_per_mass, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice', bin_info)
        fin_m1 = open('asy_vs_trento_method1_imkpkm_bin2_pidtype1_additionalcuts_binVar3_tbin{0}.txt'.format(tt), 'r')
        fin_m2 = open('asy_method2_integrated_additionalcuts_binVar3_tbin{0}.txt'.format(tt),'r')
        plot_asy_methods( can, fin_m1, fin_m2, output_pdfname, lab, 'BSA Method Comparison {0:.4f}<-t<{1:.4f}'.format(t_bins[tt-1], t_bins[tt]), '#phi_{trento} (deg)', 'Asy')

    '''
    '''

    
    # now plot averages asy as a function of -t 
    t_bin_info={ 'name':'-t',
                 'kinematic_bin' : t_bins,
                 'min_bin': 1,
                 'max_bin': len(t_bins),
                 'xlim':[0.0, 7],
                 'min_sig_bin': 0,
                 'max_sig_bin': 2,
                 'bin_centers': t_bin_centers
             }
    plot_binned_asy( can, 'asy_vs_imkpkm_tbin{}_additionalcuts_binVar3.txt', output_pdfname, lab, 'BSA vs -t via Slices in I.M. of K^{ +}K^{ -}' , t_bin_info, xtitle='-t (GeV)', ytitle='A^{sin(#phi)}_{LU}')
    
    # check over Q2, but integrated over xb, -t
    for qq in range(1,len(q2_bins)):
        asy_per_mass_q2bin = open('asy_vs_imkpkm_q2bin{}_additionalcuts_binVar3.txt'.format(qq),'r')
        plot_asy( can, asy_per_mass_q2bin, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice,'+' {0:.4f}< Q2 <{1:.4f}'.format(q2_bins[qq-1],q2_bins[qq]), bin_info)

    # now plot averages asy as a function of Q2
    q2_bin_info={ 'name':'Q^{2}',
                  'kinematic_bin' : q2_bins,
                  'min_bin': 1,
                  'max_bin': len(q2_bins),
                  'xlim':[0.0, 10],
                  'min_sig_bin': 0,
                  'max_sig_bin': 2,
                  'bin_centers': q2_bin_centers
             }
    plot_binned_asy( can, 'asy_vs_imkpkm_q2bin{}_additionalcuts_binVar3.txt', output_pdfname, lab, 'BSA vs Q^{2} via Slices in I.M. of K^{ +}K^{ -}' , q2_bin_info, xtitle='Q^{2} (GeV)',ytitle='A^{sin(#phi)}_{LU}')


    # check over xb, but integrated over Q2, -t
    
    for xx in range(1,len(xb_bins)):
        asy_per_mass_xbbin = open('asy_vs_imkpkm_xbbin{}_additionalcuts_binVar3.txt'.format(xx),'r')
        plot_asy( can, asy_per_mass_xbbin, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice,'+' {0:.4f}< Xb <{1:.4f}'.format(xb_bins[xx-1],xb_bins[xx]), bin_info)
        
    
    # now plot averages asy as a function of Xb
    xb_bin_info={ 'name':'X_{b}',
                  'kinematic_bin' : xb_bins,
                  'min_bin': 1,
                  'max_bin': len(xb_bins),
                  'xlim':[0.0, 1],
                  'min_sig_bin': 0,
                  'max_sig_bin': 2,
                  'bin_centers': xb_bin_centers
             }
    plot_binned_asy( can, 'asy_vs_imkpkm_xbbin{}_additionalcuts_binVar3.txt', output_pdfname, lab, 'BSA vs X_{b} via Slices in I.M. of K^{ +}K^{ -}' , xb_bin_info, xtitle='X_{b}',ytitle='A^{sin(#phi)}_{LU}')
    '''
    can.Print('{}]'.format(output_pdfname))
