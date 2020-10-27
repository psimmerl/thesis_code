#!/usr/bin/env python 
import argparse 
import numpy as np
import math

from array import array 
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors, TLegend, TMath,
                  TLine, kRed, kBlack, kBlue, kWhite, kGreen, kTRUE, gROOT)

from ROOT import( RooRealVar, RooArgSet, RooArgList, RooDataHist,
                  RooGaussian, RooFit, RooAddPdf, RooPolynomial, RooFormulaVar, RooLinkedList)


#####################################################################
## this macro is for examining the phi mass distribution in 
## bins of phi_trento to compare with the method in the mon-phi.py
## macro which plots asy as a function of the invariant mass of kpkm


default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
global_dump = []

gROOT.SetBatch(kTRUE)

def readInBinningInfo(fname):
    fin = open(fname,'r')
    temp = []
    for ll in fin:
        print ll
        temp.append(float(ll))
    return temp    


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

def plot_slices( canvas, histos, title_formatter, save_name, label,
                      title=None, xtitle=None, ytitle=None , plot_details=None, rebinX_option=None, field=None):

    canvas.Clear()
    canvas.Divide(3,4)
    

    root_is_dumb = []
    page_max = 0
    signal_sum=0
    upper_bound=1.1#125

    h_phi = TH1F("h_phi","h_phi", plot_details['hist_bins'], plot_details['hist_array'])#10, -180, 180)                        
    
    for bb in range(1, plot_details['max_bins']):
        canvas.cd(bb)
        #if bb == 5 or bb == 6 : continue
        trento_upper = h_phi.GetBinCenter(bb) + h_phi.GetBinWidth(bb)/2.0# plot_details['bin_center'][bb-1] + plot_details['bin_width']/2.0
        trento_lower = h_phi.GetBinCenter(bb) - h_phi.GetBinWidth(bb)/2.0#plot_details['bin_center'][bb-1] - plot_details['bin_width']/2.0
        if isinstance(histos.get(title_formatter.format(bb), default_histo), TH1F):
            h_temp = histos.get(title_formatter.format(bb), default_histo)            
            
            #if rebinX_option:
            #    h_temp.RebinX(rebinX_option[bb-1])

            print('--> Number of Entries: {}'.format(h_temp.GetEntries()))
            
            h_temp.SetLineColor(kBlack)            
            h_temp.Draw('same')
            root_is_dumb.append(h_temp)

            h_temp.SetMaximum( h_temp.GetMaximum()*1.12 )
            canvas.Update()

            h_temp.GetXaxis().SetRangeUser(0.97, upper_bound)
            
            mass_line_shape = "[0]*TMath::Gaus(x,[1],[2]) + [3] * ( x> 2*0.4937) * TMath::Power(x - 2*0.4937 , [4]) * TMath::Exp(-[5] * (x - 2*0.4937) )"
            mass_fit = TF1("mass_fit",mass_line_shape,2*0.4937,upper_bound)
            mass_fit.SetLineWidth(3)
            mass_fit.SetNpx(10000)
            exp_fit_par=7.87329
            if bb == 5 or bb == 6:
                exp_fit_par = 6
            mass_fit.SetParameters( h_temp.GetMaximum() , 1.0198 , 0.008 , 50 , 3.39580e-01 , exp_fit_par )
            mass_fit.SetParLimits(0, 0, h_temp.GetMaximum()*1.2)
            #mass_fit.SetParLimits(1, 1.018, 1.025)
            #mass_fit.SetParLimits(2, h_temp.GetBinWidth(1)/2.0, h_temp.GetRMS()*1.1)
            if plot_details['fix_fit']:        
                print(' --> fixing parameters 1 and 2 ')
                mass_fit.FixParameter(1, 1.0198)
                mass_fit.FixParameter(2, 0.00433)
            else:
                mass_fit.SetParLimits(1, 1.018, 1.025)
                mass_fit.SetParLimits(2, h_temp.GetBinWidth(1)/2.0, h_temp.GetRMS()*1.2)
                
            mass_fit.SetParLimits(4, 0, 7)
            mass_fit.SetParLimits(5, 0, 7)
            h_temp.Fit("mass_fit","R")
                        
            can.Update()

            show_back = TF1("show_back",mass_line_shape, 2*0.4937,upper_bound)
            for ip in range(0,6):
                show_back.SetParameter( ip , mass_fit.GetParameter(ip) )
            show_back.SetParameter( 0 , 0 )
            show_back.SetLineWidth(2)
            show_back.SetNpx(100)
            show_back.SetLineColor(kBlue)
            show_back.SetLineStyle(2)            
            show_back.Draw("same")

            fit_norm = mass_fit.GetParameter(0)
            fit_norm_err = mass_fit.GetParError(0)
            fit_mean = mass_fit.GetParameter(1)
            fit_mean_err = mass_fit.GetParError(1)
            fit_sigma = mass_fit.GetParameter(2)
            fit_sigma_err = mass_fit.GetParError(2)
            fit_chi2 = mass_fit.GetChisquare()
            fit_ndf = mass_fit.GetNDF()

            bin_width = h_temp.GetBinWidth(1)

            norm_const = np.sqrt(2*TMath.Pi())  / bin_width
            n_signal = fit_norm * fit_sigma * norm_const
            d_nsignal = np.sqrt( pow( fit_norm_err*fit_sigma ,2) + pow( fit_norm*fit_sigma_err ,2)  ) * norm_const
                        
            print('trento center {}'.format(plot_details['bin_center'][bb-1]))
            print('trento lower {}'.format(trento_lower))
            print('trento upper {}'.format(trento_upper))
            label.DrawLatex(0.15, 0.85,'B{0}'.format(bb))
            label.DrawLatex(0.66, 0.925,',{1:3.2f}<#phi<{2:3.2f}'.format(bb, trento_lower, trento_upper ))
            # label.DrawLatex(0.6,0.85,'N={0}'.format(int(n_signal)))
            #label.DrawLatex(0.6,0.80,'#Delta N={0}'.format(int(d_nsignal)))

            temp_lab = TLatex()
            temp_lab.SetNDC()
            temp_lab.SetTextFont(42)
            temp_lab.SetTextSize(0.082)
            temp_lab.SetTextColor(1)
            temp_lab.DrawLatex(0.6,0.83,'N={0}'.format(int(n_signal)))
            temp_lab.DrawLatex(0.6,0.74,'#Delta N={0}'.format(int(d_nsignal)))

            label.DrawLatex(0.6,0.69,'#mu = {0:4.4f}'.format(fit_mean))
            label.DrawLatex(0.6,0.64,'#sigma = {0:5.5f}'.format(fit_sigma))
            label.DrawLatex(0.6,0.59,'#Delta #sigma = {0:1.5f}'.format(fit_sigma_err))
            label.DrawLatex(0.6,0.54,'#chi^{2}'+' / NDF = {0:2.2f}/{1}'.format(fit_chi2,int(fit_ndf)))
            
            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.65, 0.025, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.055, 0.5, ytitle)
                label.SetTextAngle(0)

            if field:
                flabel=TLatex()
                flabel.SetTextColorAlpha(1, 0.376)
                flabel.SetNDC()
                flabel.SetTextFont(42)
                flabel.SetTextSize(0.08)
                flabel.DrawLatex(0.125, 0.8, field)

            h_phi.SetBinContent(bb,n_signal)
            h_phi.SetBinError(bb,d_nsignal)

            root_is_dumb.append(label)
            root_is_dumb.append(show_back)
            root_is_dumb.append(mass_fit)
            
    root_is_dumb.append(h_phi)

    can.Print(save_name)
    can.Clear()
    if "pos" in title_formatter:
        h_phi.SetLineWidth(3)
        h_phi.SetLineColor(kRed)
    else:
        h_phi.SetLineWidth(3)
        h_phi.SetLineColor(kBlue)
    h_phi.Draw("H+E")
    
    if "pos" in title_formatter:
        label.DrawLatex(0.1, 0.925,"Positive Helicity")
    else:
        label.DrawLatex(0.1, 0.925,"Negative Helicity")
    label.DrawLatex(0.65, 0.0275, '#phi_{trento} (deg)')
    label.SetTextAngle(90)
    label.DrawLatex(0.034, 0.5, 'counts')
    label.SetTextAngle(0)
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)

    canvas.Print(save_name)


def plot_slices_systematic( canvas, histos, title_formatter, save_name, label,
                            title=None, xtitle=None, ytitle=None , plot_details=None, rebinX_option=None, field=None):

    canvas.Clear()
    canvas.Divide(3,4)
    

    root_is_dumb = []
    page_max = 0
    signal_sum=0
    upper_bound=1.1#125
    
    h_phi_syst = TH1F("h_phi_syst","h_phi_syst",10, -180, 180)
    
    for bb in range(1, plot_details['max_bins']):
        canvas.cd(bb)
        #if bb == 5 or bb == 6 : continue
        trento_upper = plot_details['bin_center'][bb-1] + plot_details['bin_width']/2.0
        trento_lower = plot_details['bin_center'][bb-1] - plot_details['bin_width']/2.0
        if isinstance(histos.get(title_formatter.format(bb), default_histo), TH1F):
            h_temp = histos.get(title_formatter.format(bb), default_histo)            
            
            #if rebinX_option:
            #    h_temp.RebinX(rebinX_option[bb-1])

            print('--> Number of Entries: {}'.format(h_temp.GetEntries()))
            
            h_temp.SetLineColor(kBlack)            
            h_temp.Draw('same')
            root_is_dumb.append(h_temp)

            h_temp.SetMaximum( h_temp.GetMaximum()*1.12 )
            canvas.Update()

            h_temp.GetXaxis().SetRangeUser(0.97, upper_bound)
            
            mass_line_shape = "[0]*TMath::Gaus(x,[1],[2]) + [3] * ( x> 2*0.4937) * TMath::Power(x - 2*0.4937 , [4]) * TMath::Exp(-[5] * (x - 2*0.4937) )"
            mass_fit = TF1("mass_fit",mass_line_shape,2*0.4937,upper_bound)
            mass_fit.SetLineWidth(3)
            mass_fit.SetNpx(10000)
            exp_fit_par=7.87329
            if bb == 5 or bb == 6:
                exp_fit_par = 6
            mass_fit.SetParameters( h_temp.GetMaximum() , 1.01937e+00 , 0.008 , 50 , 3.39580e-01 , exp_fit_par )
            mass_fit.SetParLimits(0, 0, h_temp.GetMaximum()*1.2)
            #mass_fit.SetParLimits(1, 1.018, 1.025)
            #mass_fit.SetParLimits(2, h_temp.GetBinWidth(1)/2.0, h_temp.GetRMS()*1.1)
            if plot_details['fix_fit']:        
                print(' --> fixing parameters 1 and 2 ')
                mass_fit.FixParameter(1, 1.0198)
                mass_fit.FixParameter(2, 0.00433)
            else:
                mass_fit.SetParLimits(1, 1.018, 1.025)
                mass_fit.SetParLimits(2, h_temp.GetBinWidth(1)/2.0, h_temp.GetRMS()*1.2)

            mass_fit.SetParLimits(4, 0, 7)
            mass_fit.SetParLimits(5, 0, 7)
            h_temp.Fit("mass_fit","R")
                        
            can.Update()

            show_back = TF1("show_back",mass_line_shape, 2*0.4937,upper_bound)
            for ip in range(0,6):
                show_back.SetParameter( ip , mass_fit.GetParameter(ip) )
            show_back.SetParameter( 0 , 0 )
            show_back.SetLineWidth(2)
            show_back.SetNpx(100)
            show_back.SetLineColor(kBlue)
            show_back.SetLineStyle(2)            
            show_back.Draw("same")

            fit_norm = mass_fit.GetParameter(0)
            fit_norm_err = mass_fit.GetParError(0)
            fit_mean = mass_fit.GetParameter(1)
            fit_mean_err = mass_fit.GetParError(1)
            fit_sigma = mass_fit.GetParameter(2)
            fit_sigma_err = mass_fit.GetParError(2)
            fit_chi2 = mass_fit.GetChisquare()
            fit_ndf = mass_fit.GetNDF()

            bin_width = h_temp.GetBinWidth(1)

            norm_const = np.sqrt(2*TMath.Pi())  / bin_width
            n_signal = fit_norm * fit_sigma * norm_const
            d_nsignal = np.sqrt( pow( fit_norm_err*fit_sigma ,2) + pow( fit_norm*fit_sigma_err ,2)  ) * norm_const
                        
            print('trento center {}'.format(plot_details['bin_center'][bb-1]))
            print('trento lower {}'.format(trento_lower))
            print('trento upper {}'.format(trento_upper))
            label.DrawLatex(0.15, 0.85,'B{0}'.format(bb))
            label.DrawLatex(0.66, 0.925,',{1:3.2f}<#phi<{2:3.2f}'.format(bb, trento_lower, trento_upper ))
            #label.DrawLatex(0.6,0.85,'N={0}'.format(int(n_signal)))
            #label.DrawLatex(0.6,0.80,'#Delta N={0}'.format(int(d_nsignal)))

            ## try this other method subtract the bin heigh from bck at the bin center
            min_bin_kpkm_range = h_temp.FindBin(fit_mean - 3*fit_sigma)
            max_bin_kpkm_range = h_temp.FindBin(fit_mean + 3*fit_sigma)
            print( ' trying another way to calculate signal events')
            print( ' look at range {} - {} '.format(fit_mean - 3*fit_sigma, fit_mean + 3*fit_sigma))
            print( ' bin range from  {} to {}'.format(min_bin_kpkm_range, max_bin_kpkm_range))
            n_sig_total = 0 
            l_signal = []
            for bs in range( min_bin_kpkm_range, max_bin_kpkm_range):
                print(' in bin {} - bin center {}'.format(bs, h_temp.GetBinCenter(bs)))                
                bin_center = h_temp.GetBinCenter(bs)
                n_count = h_temp.GetBinContent(bs)
                n_calc_bck = show_back.Eval(h_temp.GetBinCenter(bs))
                print(' n counts {} - n bck {}'.format(n_count, n_calc_bck))
                n_sig = n_count - n_calc_bck
                if( n_count < n_calc_bck ):
                    continue
                # if the bck is greater than number in that bin just continue
                # otherwise add it to the total
                n_sig_total+=n_sig
                
                ## draw green line showing different in bin height and bck eval at bin center
                l_sig_total = TLine(bin_center, n_calc_bck, bin_center, n_count)
                l_sig_total.SetLineColor(kGreen)
                l_sig_total.SetLineWidth(2)
                l_sig_total.Draw('same')
                root_is_dumb.append(l_sig_total)

            print('=> n sig total for this phi trento bin is {}'.format(n_sig_total))
            print('=> error is sqrt of counts {}'.format( np.sqrt(n_sig_total)))
            h_phi_syst.SetBinContent(bb, n_sig_total)
            h_phi_syst.SetBinError(bb, np.sqrt(n_sig_total))
            d_n_sig_total = np.sqrt(n_sig_total)

            temp_lab = TLatex()
            temp_lab.SetNDC()
            temp_lab.SetTextFont(42)
            temp_lab.SetTextSize(0.082)
            temp_lab.SetTextColor(1)
            temp_lab.DrawLatex(0.6,0.83,'N={0}'.format(int(n_sig_total)))
            temp_lab.DrawLatex(0.6,0.74,'#Delta N={0}'.format(int(d_n_sig_total)))

            label.DrawLatex(0.6,0.69,'#mu = {0:4.4f}'.format(fit_mean))
            label.DrawLatex(0.6,0.64,'#sigma = {0:5.5f}'.format(fit_sigma))
            label.DrawLatex(0.6,0.59,'#Delta #sigma = {0:1.5f}'.format(fit_sigma_err))
            label.DrawLatex(0.6,0.54,'#chi^{2}'+' / NDF = {0:2.2f}/{1}'.format(fit_chi2,int(fit_ndf)))
            
            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.65, 0.025, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.055, 0.5, ytitle)
                label.SetTextAngle(0)

            if field:
                flabel=TLatex()
                flabel.SetTextColorAlpha(1, 0.376)
                flabel.SetNDC()
                flabel.SetTextFont(42)
                flabel.SetTextSize(0.08)
                flabel.DrawLatex(0.125, 0.8, field)

            root_is_dumb.append(label)
            root_is_dumb.append(show_back)
            root_is_dumb.append(mass_fit)
            
    # now draw results of counts in each phi trento bin
    root_is_dumb.append(h_phi_syst)

    can.Print(save_name)
    can.Clear()
    if "pos" in title_formatter:
        h_phi_syst.SetLineWidth(3)
        h_phi_syst.SetLineColor(kRed)
    else:
        h_phi_syst.SetLineWidth(3)
        h_phi_syst.SetLineColor(kBlue)
    h_phi_syst.Draw("H+E")
    
    if "pos" in title_formatter:
        label.DrawLatex(0.1, 0.925,"Positive Helicity Systematic")
    else:
        label.DrawLatex(0.1, 0.925,"Negative Helicity Systematic")
    label.DrawLatex(0.65, 0.0275, '#phi_{trento} (deg)')
    label.SetTextAngle(90)
    label.DrawLatex(0.034, 0.5, 'counts')
    label.SetTextAngle(0)
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)

    canvas.Print(save_name)


def plot_page_bins( canvas, histos, title_formatter, save_name, label,
                      title=None, xtitle=None, ytitle=None , plot_details=None, rebinX_option=None, field=None):

    canvas.Clear()

    root_is_dumb = []    
    pol=1
    if 'bp' in plot_details.keys():
        pol= plot_details['bp']

    g_list = []
    ncols=3
    nrows=int( (plot_details['max_bins']-1)/ncols)
    print(' number of cols {} number of rows {}'.format(ncols,nrows))
    canvas.Divide(ncols, nrows)

    for bb in range(1, plot_details['max_bins']):
        canvas.cd(bb)
        if isinstance(histos.get(title_formatter.format('pos',bb), default_histo), TH1F):
            print(' --> mass bin {}'.format(bb) )
            h_temp_pos = histos.get(title_formatter.format('pos',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('pos',bb)).GetTitle())
            h_temp_neg = histos.get(title_formatter.format('neg',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('neg',bb)).GetTitle())
            

            h_temp_pos.SetLineColor(kRed)
            h_temp_neg.SetLineColor(kBlue)

            leg=TLegend(0.5, 0.7, 0.7, 0.9)
            leg.AddEntry(h_temp_pos,'pos','l')
            leg.AddEntry(h_temp_neg,'neg','l')
            
            h_temp_pos.Draw()
            h_temp_neg.Draw('same')
            leg.Draw('same')

            root_is_dumb.append(h_temp_pos)
            root_is_dumb.append(h_temp_neg)
            root_is_dumb.append(leg)

            
            if 'bin_list' in plot_details.keys():
                label.DrawLatex(0.2, 0.925, 'BinRange:{0:.2f}<im<{1:.2f}'.format(plot_details['bin_list'][bb-1], plot_details['bin_list'][bb]))            

            if( h_temp_pos.GetMaximum() >  h_temp_neg.GetMaximum() ): 
                h_temp_pos.SetMaximum( h_temp_pos.GetMaximum() + 50)
            else:                    
                h_temp_neg.SetMaximum( h_temp_neg.GetMaximum() + 50 )
                

            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.5, 0.015, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.055, 0.5, ytitle)
                label.SetTextAngle(0)
        

    canvas.Print(save_name)            
    

def fit_phi_mass( h_temp, bb, set_fix, field=None, bin_info=None):
    dump = []
    ## define fit and fit range
    upper_bound=1.125
    h_temp.GetXaxis().SetRangeUser(0.97, upper_bound)
    mass_line_shape = "[0]*TMath::Gaus(x,[1],[2]) + [3] * ( x> 2*0.4937) * TMath::Power(x - 2*0.4937 , [4]) * TMath::Exp(-[5] * (x - 2*0.4937) )"
    mass_fit = TF1("mass_fit",mass_line_shape,2*0.4937,upper_bound)
    
    # set and fix the parameters 
    exp_fit_par=7.87329
    if bb == 5 or bb == 6:
        exp_fit_par = 6
    mass_fit.SetParameters( h_temp.GetMaximum() , 1.01938e+00 , 0.008 , 50 , 3.39580e-01 , exp_fit_par )
    
    if set_fix:        
        print(' --> fixing parameters 1 and 2 ')
        mass_fit.FixParameter(1, 1.0198)
        mass_fit.FixParameter(2, 0.00433)
    else:
        mass_fit.SetParLimits(1, 1.0178, 1.025)
        mass_fit.SetParLimits(2, h_temp.GetBinWidth(1)/2.0, h_temp.GetRMS()*1.2)
    
    mass_fit.SetParLimits(0,0, h_temp.GetMaximum()*1.2)    
    mass_fit.SetParLimits(4, 0,7)
    mass_fit.SetParLimits(5, 0,7)
    h_temp.Fit("mass_fit","R")
    
        
    show_back = TF1("show_back",mass_line_shape, 2*0.4937,upper_bound)
    for ip in range(0,6):
        show_back.SetParameter( ip , mass_fit.GetParameter(ip) )
    show_back.SetParameter( 0 , 0 )
            
    fit_norm = mass_fit.GetParameter(0)
    fit_norm_err = mass_fit.GetParError(0)
    fit_mean = mass_fit.GetParameter(1)
    fit_mean_err = mass_fit.GetParError(1)
    fit_sigma = mass_fit.GetParameter(2)
    fit_sigma_err = mass_fit.GetParError(2)
    fit_chi2 = mass_fit.GetChisquare()
    fit_ndf = mass_fit.GetNDF()
    
    bin_width = h_temp.GetBinWidth(1)
    
    norm_const = np.sqrt(2*TMath.Pi())  / bin_width
    n_signal = fit_norm * fit_sigma * norm_const
    d_nsignal = np.sqrt( pow( fit_norm_err*fit_sigma ,2) + pow( fit_mean*fit_sigma_err ,2)  ) * norm_const
    dump.append(h_temp)
    dump.append(mass_fit)

    return [n_signal, d_nsignal, fit_chi2, fit_ndf, show_back, fit_mean, fit_sigma]

            
    
def plot_asy( canvas, histos, title_formatter, save_name, label,
              title=None, xtitle=None, ytitle=None , plot_details=None, f_out=None, rebinX_option=None, field=None):

    canvas.Clear()
    canvas.Update()

    root_is_dumb = []    
    pol=1
    if 'bp' in plot_details.keys():
        pol= plot_details['bp']

    canvas.Divide(1,1)

    fout=None
    if f_out:
        fout = open(f_out,'w')


    trento_asy = []
    x_center = array('d')
    y_val = array('d')
    x_err = array('d')
    y_err = array('d')

    y_quality_pos = []
    y_quality_neg = []
    
    for bb in range(1, plot_details['max_bins']):
        canvas.cd(bb)
        #if bb == 5 or bb == 6 : continue
        print('- - - - -------------> bb is {}'.format(bb))
        if isinstance(histos.get(title_formatter.format('pos',bb), default_histo), TH1F) and isinstance(histos.get(title_formatter.format('neg',bb), default_histo), TH1F) :
            print('--> Trento Center {}'.format(plot_details['bin_center'][bb-1]))
            trento_upper = plot_details['bin_center'][bb-1] + plot_details['bin_width']/2.0
            trento_lower = plot_details['bin_center'][bb-1] - plot_details['bin_width']/2.0
            print('--> Trento Range {0}< phi <{1}'.format(trento_lower, trento_upper))
            
            h_temp_pos = histos.get(title_formatter.format('pos',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('pos',bb)).GetTitle())
            h_temp_neg = histos.get(title_formatter.format('neg',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('neg',bb)).GetTitle())                                                    
            #if rebinX_option:
            #    h_temp_pos.RebinX(rebinX_option[0][bb-1])
            #    h_temp_neg.RebinX(rebinX_option[1][bb-1])

            print('--> POS Helicity Number of Entries: {}'.format(h_temp_pos.GetEntries()))
            print('--> NEG Helicity Number of Entries: {}'.format(h_temp_pos.GetEntries()))
            root_is_dumb.append(h_temp_pos)
            root_is_dumb.append(h_temp_neg)
            
            ##########################
                        
            pos_results = fit_phi_mass( h_temp_pos, bb, plot_details['fix_fit'])
            neg_results = fit_phi_mass( h_temp_neg, bb, plot_details['fix_fit'])
    
            ## save quality of fit results for plotting on asymmetry plot per phi trento bin to gauge results
            y_pos_temp = []
            y_pos_temp.append(plot_details['bin_center'][bb-1])
            y_pos_temp.append(pos_results[2:])

            y_neg_temp = []
            y_neg_temp.append(plot_details['bin_center'][bb-1])
            y_neg_temp.append(neg_results[2:])

            y_quality_pos.append(y_pos_temp)
            y_quality_neg.append(y_neg_temp)

            print("  Positive Results " )
            print(pos_results)

            print("  Negative Results " )
            print(neg_results)

            P = pos_results[0]
            N = neg_results[0]

            nom = P-N
            den = P+N

            y=0
            yerr = 0
            print(" number of positive values %d, number of negative values %d" % (P, N))        
            if(den==0):
                y=0
                yerr=0
            else:
                y=nom/den
                yerr = math.sqrt( (1 - y**2)/(P+N))


            print(" bin center {}, asy {}, bin center error {}, asy error {} ".format(plot_details['bin_center'][bb-1], y, 0, yerr))
            x_center.append(plot_details['bin_center'][bb-1])
            y_val.append(y*(1.0/pol))
            x_err.append(0)
            y_err.append(yerr)
            
            ## write asy vs phi_trento out to compare with otehr methods
            if f_out:
                fout.write(str(plot_details['bin_center'][bb-1]) + ' ' + str(y*(1.0/pol)) + ' ' + str(yerr) + ' \n')


        ##########################
        ## now draw the results for each trento phi bin in a graph
        graph = TGraphErrors(len(x_center), x_center, y_val, x_err, y_err)
        graph.SetMarkerStyle(22)
        graph.GetHistogram().SetMaximum(0.65)
        graph.GetHistogram().SetMinimum(-0.65)
        graph.GetXaxis().SetRangeUser(-180.0, 180.0)
        graph.SetMarkerColor(kBlack)
        graph.SetMarkerSize(9)
        root_is_dumb.append(graph)
        
        # draw grid
        h = gPad.DrawFrame(-180, -0.65, -180, 0.65)
        h.GetXaxis().SetNdivisions(-7)
        h.GetYaxis().SetNdivisions(-7)
        h.GetXaxis().SetLabelSize(0.025)
        h.GetYaxis().SetLabelSize(0.025)
        gPad.SetGrid()
        root_is_dumb.append(h)
        
        graph.Draw('AP')
        gPad.RedrawAxis("g")
        
        fit_asy = TF1('fasy{}'.format(bb),"[0]*sin(x*TMath::Pi()/180.0)", -180, 180)
        fit_asy.SetParameter(0,1)
        graph.Fit('fasy{}'.format(bb))
        fit_asy.SetLineColor(kRed)
        fit_asy.SetLineWidth(3)
        fit_asy_chi2 = fit_asy.GetChisquare()
        fit_asy_ndf = fit_asy.GetNDF()
        root_is_dumb.append(fit_asy)
                        
        # draw bin related info
        temp_lab = TLatex()
        temp_lab.SetNDC()
        temp_lab.SetTextFont(42)
        temp_lab.SetTextSize(0.06)
        temp_lab.SetTextColor(1)
        # for poster change back
        temp_lab.DrawLatex(0.55, 0.762,  'Asy_{s}'+': {0:.4f}'.format(fit_asy.GetParameter(0)))
        temp_lab.DrawLatex(0.55, 0.695,  '#Delta Asy_{s}'+': {0:.4f}'.format(fit_asy.GetParError(0)))        
        temp_lab.DrawLatex(0.55, 0.625,  '#chi^{2}'+'/NDF:{0:.1f}/{1}'.format(fit_asy_chi2,fit_asy_ndf))
        
        #label.DrawLatex(0.55, 0.755,  'Asy:{0:.3f}'.format(fit_asy.GetParameter(0)))
        #label.DrawLatex(0.55, 0.695,  'AsyErr:{0:.4f}'.format(fit_asy.GetParError(0)))        
        #label.DrawLatex(0.55, 0.635,  '#chi^{2}'+'/NDF:{0:.1f}/{1}'.format(fit_asy_chi2,fit_asy_ndf))
        
        
        if title:
            label.DrawLatex(0.1, 0.925, title)
            
        if xtitle:
            label.DrawLatex(0.5, 0.0275, xtitle)
            
        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.055, 0.5, ytitle)
            label.SetTextAngle(0)
            
        if field:
            flabel=TLatex()
            flabel.SetTextColorAlpha(1, 0.376)
            flabel.SetNDC()
            flabel.SetTextFont(42)
            flabel.SetTextSize(0.08)
            flabel.DrawLatex(0.125, 0.847, field)


        #### draw the quality makers on the graph
        ## change for poster contest put back
        '''
        label.DrawLatex(0.8,0.85,'#color[2]{POS}')
        label.DrawLatex(0.8,0.79,'#color[4]{NEG}')
        label.SetTextAngle(60)            
        for qq in range(0, len(y_quality_pos)):
            print(y_quality_pos[qq])
            print( y_quality_pos[qq][1][0])
            print( y_quality_pos[qq][1][1])
            label.DrawLatex( (y_quality_pos[qq][0]+180 - 10)/360, 0.2,"{},".format(qq) + "#color[2]{"+ "{0:2.1f} / {1}".format(y_quality_pos[qq][1][0], y_quality_pos[qq][1][1]) + "}")
            label.DrawLatex( (y_quality_neg[qq][0]+180 - 10)/360 + 0.06 - 0.01, 0.2, "{},".format(qq) + "#color[4]{"+ "{0:2.1f} / {1}".format(y_quality_neg[qq][1][0], y_quality_neg[qq][1][1]) + "}")
            root_is_dumb.append(label)
        label.SetTextAngle(0)
        '''

    if f_out:
        fout.close()
    canvas.Print(save_name)


### systematic calc for this method
def plot_asy_systematic( canvas, histos, title_formatter, save_name, label,
                         title=None, xtitle=None, ytitle=None , plot_details=None, f_out=None, rebinX_option=None, field=None):

    canvas.Clear()
    canvas.Update()

    root_is_dumb = []    
    pol=1
    if 'bp' in plot_details.keys():
        pol= plot_details['bp']

    canvas.Divide(1,1)

    fout=None
    if f_out:
        fout = open(f_out,'w')


    trento_asy = []
    x_center = array('d')
    y_val = array('d')
    x_err = array('d')
    y_err = array('d')

    y_quality_pos = []
    y_quality_neg = []
    
    for bb in range(1, plot_details['max_bins']):
        canvas.cd(bb)
        #if bb == 5 or bb == 6 : continue
        if isinstance(histos.get(title_formatter.format('pos',bb), default_histo), TH1F) and isinstance(histos.get(title_formatter.format('neg',bb), default_histo), TH1F) :
            print('--> Trento Center {}'.format(plot_details['bin_center'][bb-1]))
            trento_upper = plot_details['bin_center'][bb-1] + plot_details['bin_width']/2.0
            trento_lower = plot_details['bin_center'][bb-1] - plot_details['bin_width']/2.0
            print('--> Trento Range {0}< phi <{1}'.format(trento_lower, trento_upper))
            
            h_temp_pos = histos.get(title_formatter.format('pos',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('pos',bb)).GetTitle())
            h_temp_neg = histos.get(title_formatter.format('neg',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('neg',bb)).GetTitle())                                                    
            #if rebinX_option:
            #    h_temp_pos.RebinX(rebinX_option[0][bb-1])
            #    h_temp_neg.RebinX(rebinX_option[1][bb-1])

            print('--> POS Helicity Number of Entries: {}'.format(h_temp_pos.GetEntries()))
            print('--> NEG Helicity Number of Entries: {}'.format(h_temp_pos.GetEntries()))
            root_is_dumb.append(h_temp_pos)
            root_is_dumb.append(h_temp_neg)
            
            ##########################
            pos_results = fit_phi_mass( h_temp_pos, bb, plot_details['fix_fit'])
            neg_results = fit_phi_mass( h_temp_neg, bb, plot_details['fix_fit'])
            
            ## calc difference between bin peak and bck eval at bin center +/- 3 sigma
            pos_mean = pos_results[5]
            pos_sig = pos_results[6]

            neg_mean = neg_results[5]
            neg_sig = neg_results[6]

            min_bin_kpkm_pos = h_temp_pos.FindBin(pos_mean - 3*pos_sig)
            max_bin_kpkm_pos = h_temp_pos.FindBin(pos_mean + 3*pos_sig)

            min_bin_kpkm_neg = h_temp_neg.FindBin(neg_mean - 3*neg_sig)
            max_bin_kpkm_neg = h_temp_neg.FindBin(neg_mean + 3*neg_sig)
            # first calc for positive helicty state
            n_sig_total_pos=0
            for bs_pos in range(min_bin_kpkm_pos, max_bin_kpkm_pos):
                print('POS: in bin {} - bin center {}'.format(bs_pos, h_temp_pos.GetBinCenter(bs_pos)))                
                bin_center = h_temp_pos.GetBinCenter(bs_pos)
                n_count = h_temp_pos.GetBinContent(bs_pos)
                n_calc_bck = pos_results[4].Eval(h_temp_pos.GetBinCenter(bs_pos))
                print('POS: n counts {} - n bck {}'.format(n_count, n_calc_bck))
                n_sig = n_count - n_calc_bck
                if( n_count < n_calc_bck ):
                    continue
                # if the bck is greater than number in that bin just continue
                # otherwise add it to the total
                n_sig_total_pos+=n_sig

            # second calc for negative helicty state
            n_sig_total_neg=0
            for bs_neg in range(min_bin_kpkm_neg, max_bin_kpkm_neg):
                print('NEG: in bin {} - bin center {}'.format(bs_neg, h_temp_neg.GetBinCenter(bs_neg)))                
                bin_center = h_temp_neg.GetBinCenter(bs_neg)
                n_count = h_temp_neg.GetBinContent(bs_neg)
                n_calc_bck = neg_results[4].Eval(h_temp_neg.GetBinCenter(bs_neg))
                print('NEG: n counts {} - n bck {}'.format(n_count, n_calc_bck))
                n_sig = n_count - n_calc_bck
                if( n_count < n_calc_bck ):
                    continue
                # if the bck is greater than number in that bin just continue
                # otherwise add it to the total
                n_sig_total_neg+=n_sig
                
            print(' Number of positive {}, and Number of Negative {}'.format(n_sig_total_pos, n_sig_total_neg))
                
            ## save quality of fit results for plotting on asymmetry plot per phi trento bin to gauge results
            y_pos_temp = []
            y_pos_temp.append(plot_details['bin_center'][bb-1])
            y_pos_temp.append(pos_results[2:])

            y_neg_temp = []
            y_neg_temp.append(plot_details['bin_center'][bb-1])
            y_neg_temp.append(neg_results[2:])

            y_quality_pos.append(y_pos_temp)
            y_quality_neg.append(y_neg_temp)

            print("  Positive Results " )
            print(n_sig_total_pos)

            print("  Negative Results " )
            print(n_sig_total_neg)

            P = n_sig_total_pos
            N = n_sig_total_neg

            nom = P-N
            den = P+N

            y=0
            yerr = 0
            print(" number of positive values %d, number of negative values %d" % (P, N))        
            if(den==0):
                y=0
                yerr=0
            else:
                y=nom/den
                yerr = math.sqrt( (1 - y**2)/(P+N))


            print(" bin center {}, asy {}, bin center error {}, asy error {} ".format(plot_details['bin_center'][bb-1], y, 0, yerr))
            x_center.append(plot_details['bin_center'][bb-1])
            y_val.append(y*(1.0/pol))
            x_err.append(0)
            y_err.append(yerr)
            
            ## write asy vs phi_trento out to compare with otehr methods
            if f_out:
                fout.write(str(plot_details['bin_center'][bb-1]) + ' ' + str(y*(1.0/pol)) + ' ' + str(yerr) + ' \n')


        ##########################
        ## now draw the results for each trento phi bin in a graph
        graph = TGraphErrors(len(x_center), x_center, y_val, x_err, y_err)
        graph.SetMarkerStyle(22)
        graph.GetHistogram().SetMaximum(0.65)
        graph.GetHistogram().SetMinimum(-0.65)
        graph.GetXaxis().SetRangeUser(-180.0, 180.0)
        graph.SetMarkerColor(kBlack)
        graph.SetMarkerSize(9)
        root_is_dumb.append(graph)
        
        # draw grid
        h = gPad.DrawFrame(-180, -0.65, -180, 0.65)
        h.GetXaxis().SetNdivisions(-7)
        h.GetYaxis().SetNdivisions(-7)
        h.GetXaxis().SetLabelSize(0.025)
        h.GetYaxis().SetLabelSize(0.025)
        gPad.SetGrid()
        root_is_dumb.append(h)
        
        graph.Draw('AP')
        gPad.RedrawAxis("g")
        
        fit_asy = TF1('fasy{}'.format(bb),"[0]*sin(x*TMath::Pi()/180.0)", -180, 180)
        fit_asy.SetParameter(0,1)
        graph.Fit('fasy{}'.format(bb))
        fit_asy.SetLineColor(kRed)
        fit_asy.SetLineWidth(3)
        fit_asy_chi2 = fit_asy.GetChisquare()
        fit_asy_ndf = fit_asy.GetNDF()
        root_is_dumb.append(fit_asy)
                        
        # draw bin related info
        temp_lab = TLatex()
        temp_lab.SetNDC()
        temp_lab.SetTextFont(42)
        temp_lab.SetTextSize(0.06)
        temp_lab.SetTextColor(1)
        # for poster change back
        temp_lab.DrawLatex(0.55, 0.762,  'Asy_{s}'+': {0:.4f}'.format(fit_asy.GetParameter(0)))
        temp_lab.DrawLatex(0.55, 0.695,  '#Delta Asy_{s}'+': {0:.4f}'.format(fit_asy.GetParError(0)))        
        temp_lab.DrawLatex(0.55, 0.625,  '#chi^{2}'+'/NDF:{0:.1f}/{1}'.format(fit_asy_chi2,fit_asy_ndf))
        
        #label.DrawLatex(0.55, 0.755,  'Asy:{0:.3f}'.format(fit_asy.GetParameter(0)))
        #label.DrawLatex(0.55, 0.695,  'AsyErr:{0:.4f}'.format(fit_asy.GetParError(0)))        
        #label.DrawLatex(0.55, 0.635,  '#chi^{2}'+'/NDF:{0:.1f}/{1}'.format(fit_asy_chi2,fit_asy_ndf))
        
        
        if title:
            label.DrawLatex(0.1, 0.925, title)
            
        if xtitle:
            label.DrawLatex(0.5, 0.0275, xtitle)
            
        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.055, 0.5, ytitle)
            label.SetTextAngle(0)
            
        if field:
            flabel=TLatex()
            flabel.SetTextColorAlpha(1, 0.376)
            flabel.SetNDC()
            flabel.SetTextFont(42)
            flabel.SetTextSize(0.08)
            flabel.DrawLatex(0.125, 0.847, field)


        #### draw the quality makers on the graph
        ## change for poster contest put back
        '''
        label.DrawLatex(0.8,0.85,'#color[2]{POS}')
        label.DrawLatex(0.8,0.79,'#color[4]{NEG}')
        label.SetTextAngle(60)            
        for qq in range(0, len(y_quality_pos)):
            print(y_quality_pos[qq])
            print( y_quality_pos[qq][1][0])
            print( y_quality_pos[qq][1][1])
            label.DrawLatex( (y_quality_pos[qq][0]+180 - 10)/360, 0.2,"{},".format(qq) + "#color[2]{"+ "{0:2.1f} / {1}".format(y_quality_pos[qq][1][0], y_quality_pos[qq][1][1]) + "}")
            label.DrawLatex( (y_quality_neg[qq][0]+180 - 10)/360 + 0.06 - 0.01, 0.2, "{},".format(qq) + "#color[4]{"+ "{0:2.1f} / {1}".format(y_quality_neg[qq][1][0], y_quality_neg[qq][1][1]) + "}")
            root_is_dumb.append(label)
        label.SetTextAngle(0)
        '''

    if f_out:
        fout.close()
    canvas.Print(save_name)





def add_text_page(can, label, text, save_name):
    """ Write some text onto a page. """
    can.Clear()
    can.cd(1)
    label.DrawLatex(0.1, 0.5, text)
    can.Print(save_name)
    
    
if __name__ == '__main__':

    ap = argparse.ArgumentParser()
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
    output_pdfname = args.output_prefix + '.pdf'
    rootfile = TFile(input_rootfile)
    histos = load_histos(rootfile)

    # uncomment to view all histograms in file
    #for k,v in histos.items():
    #    print(k, v)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 5000, 5000)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    bp=0
    field_setting='inbNoutb'

    if field_setting == 'inb':
        bp = 0.8592
    elif field_setting == 'outb':
        bp = 0.8920 # wein angle same as end of inbending data
    elif field_setting == 'inbNoutb':
        bp=0.874

    can.Print('{}['.format(output_pdfname))

    ### get the bin centers for the phi trento distribution
    trento_centers = [histos['trento_pos_mass_bin1'].GetBinCenter(xx) for xx in range(1, len(histos['trento_pos_mass_bin1'])-1)]
    print(trento_centers)

    trento_sliced = { 'max_bins':11, # 13 for 12 bins
                      'colors': True,                      
                      'bin_center': trento_centers,
                      'bin_width':histos['trento_pos_mass_bin2'].GetBinWidth(1),
                      'fix_fit':True,
                      'hist_bins':10,
                      'hist_array': array('d',[-180, -144, -108, -72, -36, 0, 36, 72, 108, 144, 180])
                    }
    
    plot_slices( can, histos, 'trento_bin{}_pos' , output_pdfname, lab,
                      title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced, field=field_setting)
    
    plot_slices( can, histos, 'trento_bin{}_neg' , output_pdfname, lab,
                      title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced, field=field_setting)
    
    plot_asy( can, histos, 'trento_bin{1}_{0}', output_pdfname, lab,
              title='Asy - Method 2, Remove #Lambda(1520),#Lambda(1820)', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=trento_sliced, field=field_setting)
              #f_out='asy_method2_integrated_additionalcuts_binVar10_rmvres12_'+field_setting+'.txt', field=field_setting)

    '''
    
    ### systematic check for the method 2
    ## this method subtracts the bck value from the bin value within +/- 3 sigma of peak.
    print(' >>>>>>>>>>>>>>>> starting systematic to subtract background bin by bin')
    plot_slices_systematic( can, histos, 'trento_bin{}_pos' , output_pdfname, lab,
                            title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced, field=field_setting)
    
    plot_slices_systematic( can, histos, 'trento_bin{}_neg' , output_pdfname, lab,
                            title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced, field=field_setting)

    plot_asy_systematic( can, histos, 'trento_bin{1}_{0}', output_pdfname, lab,
                         title='Asy - Method 2 Systematic a', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=trento_sliced, 
                         f_out='asy_method2_systematic_integrated_additionalcuts_binVar10_rmvres12_'+field_setting+'.txt', field=field_setting)


    ### systematic check of method 2
    ## this method allows the fit parameters to vary, some of the parameters are only allowed to vary within a small window 
    ## otherwise the stability of the fits can vary rapidly
    ## only change the value of the fix_fit key , keep everything else the same
    trento_sliced['fix_fit'] = False
    
    
    plot_slices( can, histos, 'trento_bin{}_pos' , output_pdfname, lab,
                      title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced, field=field_setting)
    
    plot_slices( can, histos, 'trento_bin{}_neg' , output_pdfname, lab,
                      title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced, field=field_setting)
    
    plot_asy( can, histos, 'trento_bin{1}_{0}', output_pdfname, lab,
              title='Asy - Method 2 Systematic b, No Resonances', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=trento_sliced, 
              f_out='asy_method2_systematic2b_additionalcuts_binVar10_rmvres12_'+field_setting+'.txt', field=field_setting)

    
    ## done with main techniques
    ####


    #################################################################
    ### analysis check on bin effects


    ## systematic check of method 2
    ## 12 bins in trento
    ## checks the effect of bin migration on asy result near +/- 90 deg
    ## noticed this from simulation - so have to check in data if it is a big deal,
    ## results of order of half a percent or less.
    print(' ' )
    print(' ' )
    print(' >>>>>>>>>> starting systematic study 2')
    trento_centers_mod2 = [histos['mod2_trento_pos_mass_bin1'].GetBinCenter(xx) for xx in range(1, len(histos['mod2_trento_pos_mass_bin1'])-1)]
    print(trento_centers_mod2)

    trento_sliced_mod2 = { 'max_bins':13, # 13 for 12 bins
                           'colors': True,                      
                           'bin_center': trento_centers_mod2,
                           'bin_width':histos['mod2_trento_pos_mass_bin2'].GetBinWidth(1),
                           'fix_fit':True,
                           'hist_bins':12,
                           'hist_array': array('d', [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
                   }
    
    plot_slices( can, histos, 'trento_bin{}_config2_pos' , output_pdfname, lab,
                      title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced_mod2, field=field_setting)
    
    plot_slices( can, histos, 'trento_bin{}_config2_neg' , output_pdfname, lab,
                      title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced_mod2, field=field_setting)
    
    plot_asy( can, histos, 'trento_bin{1}_config2_{0}', output_pdfname, lab,
              title='Asy - Method 2, No Resonances, 12 #phi_{trento} bins', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=trento_sliced_mod2, 
              f_out='asy_method2_integrated_additionalcuts_binVar10_rmvres12_config2_'+field_setting+'.txt', field=field_setting)

    ########################################################33
    ## next systematic check of phi trento using config3 / mod3
    ## merge the two middle bins from -30 to 30 for 10 trento bins, so [-180, -144, -108, -72, -36, 36, 72, 108, 144]
    print(' ' )
    print(' ' )
    print(' >>>>>>>>>> starting systematic study 3')
    trento_centers_mod3 = [histos['mod3_trento_pos_mass_bin1'].GetBinCenter(xx) for xx in range(1, len(histos['mod3_trento_pos_mass_bin1'])-1)]
    print(trento_centers_mod3)

    trento_sliced_mod3 = { 'max_bins': 10, #10, # 13 for 12 bins
                           'colors': True,                      
                           'bin_center': trento_centers_mod3,
                           'bin_width':histos['mod3_trento_pos_mass_bin2'].GetBinWidth(1),
                           'fix_fit':True,
                           'hist_bins': 9,
                           'hist_array': array('d',[-180, -144, -108, -72, -36, 36, 72, 108, 144, 180])

                       }
    
    plot_slices( can, histos, 'trento_bin{}_config3_pos' , output_pdfname, lab,
                      title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced_mod3, field=field_setting)
    
    plot_slices( can, histos, 'trento_bin{}_config3_neg' , output_pdfname, lab,
                      title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced_mod3, field=field_setting)
    
    #trento_bin9_config3_pos
    plot_asy( can, histos, 'trento_bin{1}_config3_{0}', output_pdfname, lab,
              title='Asy - Method 2, No Resonances, (10 bins) Merge [-30,30] #phi_{trento} bins', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=trento_sliced_mod3, 
              f_out='asy_method2_integrated_additionalcuts_binVar10_rmvres12_config3_'+field_setting+'.txt', field=field_setting)

    ########################################################33
    ## next systematic check of phi trento using config4 / mod4
    ## merge the two middle bins from -30 to 30 for 10 trento bins, so [-180, -150, -120, -90, -60, -30, 30, 60, 90, 120, 150]
    print(' ' )
    print(' ' )
    print(' >>>>>>>>>> starting systematic study 4')

    trento_centers_mod4 = [histos['mod4_trento_pos_mass_bin1'].GetBinCenter(xx) for xx in range(1, len(histos['mod4_trento_pos_mass_bin1'])-1)]
    print(trento_centers_mod4)

    trento_sliced_mod4 = { 'max_bins': 12, # bins + 1
                           'colors': True,                      
                           'bin_center': trento_centers_mod4,
                           'bin_width':histos['mod4_trento_pos_mass_bin2'].GetBinWidth(1),
                           'fix_fit':True,
                           'hist_bins': 11,
                           'hist_array': array('d',[-180, -150, -120, -90, -60, -30, 30, 60, 90, 120, 150, 180])
                       }
    
    plot_slices( can, histos, 'trento_bin{}_config4_pos' , output_pdfname, lab,
                      title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced_mod4, field=field_setting)
    
    plot_slices( can, histos, 'trento_bin{}_config4_neg' , output_pdfname, lab,
                      title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced_mod4, field=field_setting)

    plot_asy( can, histos, 'trento_bin{1}_config4_{0}', output_pdfname, lab,
              title='Asy - Method 2, No Resonances, (12 bins) Merge [-30,30] #phi_{trento} bins', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=trento_sliced_mod4, 
              f_out='asy_method2_integrated_additionalcuts_binVar10_rmvres12_config4_'+field_setting+'.txt', field=field_setting)

    ## 
    ## end anaylsis

    '''






    
    '''
    mass_bins = readInBinningInfo('bin_mass_rangesVar9_'+field_setting+'.txt')
    q2_bins=readInBinningInfo('bin_limits_q2_'+field_setting+'.txt')
    xb_bins=readInBinningInfo('bin_limits_xb_'+field_setting+'.txt')
    t_bins=readInBinningInfo('bin_limits_t_'+field_setting+'.txt')
    
    q2_centers = readInBinningInfo('bin_centers_q2_'+field_setting+'.txt')
    xb_centers = readInBinningInfo('bin_centers_xb_'+field_setting+'.txt')
    t_centers = readInBinningInfo('bin_centers_t_'+field_setting+'.txt')


    t_im_asy_details={'max_bins':11,
                      'colors': True,
                      'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                      'bin_list_fx': [1.0039, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                      1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                      'bin_list_og':[0.95, 1.01000, 1.01466, 1.01933, 1.024, 1.02866, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                     1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                      'bin_list':[0.95, 1.01000, 1.02866, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                  1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                      'bin_center': trento_centers,
                      'bin_width':histos['trento_pos_mass_bin2'].GetBinWidth(1),
                      'rebinX_pos' : [ [1, 1, 1, 1, 2, 2, 1, 1, 1, 1], [1, 1, 1, 2, 2, 2, 1, 2, 1, 1]],
                      'rebinX_neg' : [ [1, 1, 1, 1, 1, 3, 1, 1, 1, 1], [1, 1, 1, 2, 2, 2, 2, 1, 1, 2]],
                      'fix_fit': True
                }

    for tt in range(1,len(t_bins)):
        print('------------------------------------------------------------> creating asymettry plots for t bin {}'.format(tt))
    
        plot_slices( can, histos, 'trento_bin{}_pos_t_bin'+str(tt) , output_pdfname, lab,
                     title='FD:Inv, Mass, POS, {0:.3f}<-t<{1:.3f}'.format(t_bins[tt-1],t_bins[tt]), xtitle='K^{ +}  K^{ -} (GeV)', ytitle='counts', 
                     plot_details=trento_sliced, rebinX_option=t_im_asy_details['rebinX_pos'][tt-1] )

        plot_slices( can, histos, 'trento_bin{}_neg_t_bin'+str(tt) , output_pdfname, lab,
                     title='FD:Inv, Mass, NEG, {0:.3f}<-t<{1:.3f}'.format(t_bins[tt-1],t_bins[tt]), xtitle='K^{ +}  K^{ -} (GeV)', ytitle='counts', 
                     plot_details=trento_sliced, rebinX_option=t_im_asy_details['rebinX_neg'][tt-1] )

        plot_asy( can, histos, 'trento_bin{1}_{0}'+'_t_bin'+str(tt), output_pdfname, lab,
                  title='Asy,{0:.3f}<-t<{1:.3f}'.format(t_bins[tt-1],t_bins[tt]), xtitle='#phi_{trento} (deg)', ytitle='Asy', plot_details = t_im_asy_details,
                  f_out='asy_method2_integrated_additionalcuts_binVar8_{1}_tbin{0}.txt'.format(tt, field_setting), rebinX_option = [t_im_asy_details['rebinX_pos'][tt-1], t_im_asy_details['rebinX_neg'][tt-1]] )
                
    ''' 
    '''
    #######################################
    
    q2_im_asy_details={'max_bins':15,
                       'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                       'bin_list_fx': [1.0039, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                       1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                       'bin_list_og':[0.95, 1.01000, 1.01466, 1.01933, 1.024, 1.02866, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                      1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                       'bin_list':[0.95, 1.01000, 1.02866, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                   1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                       'max': 1100,
                       'bins': q2_bins,
                       'bin_centers':q2_centers
                       
                }

    ## q2 bins
    for qq in range(1,len(q2_bins)):
        print('------------------------------------------------------------> creating asymettry plots for Q2 bin {}'.format(qq))
        plot_slices( can, histos, 'q2_bin'+str(qq)+'_mass_bin{}' , output_pdfname, lab,
                      title='FD:Final, Inv.MassBins. {0:.3f}<Q2<{1:.3f}'.format(q2_bins[qq-1], q2_bins[qq]), xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_sliced)

        plot_asy( can, histos, 'trento_{}_q2_bin'+str(qq)+'_mass_bin{}', output_pdfname, lab,
                  title='Asy,{0:.3f}<Q2<{1:.3f}'.format(q2_bins[qq-1], q2_bins[qq]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=q2_im_asy_details, f_out='asy_vs_imkpkm_q2bin{0}_additionalcuts_binVar6_{1}.txt'.format(qq, field_setting))

    #######################################
    xb_im_asy_details={'max_bins':15,
                       'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                       'bin_list_fx': [1.0039, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                       1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                       'bin_list_og': [0.95, 1.01000, 1.01466, 1.01933, 1.024, 1.02866, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                       1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                       'bin_list':[0.95, 1.01000, 1.02866, 1.0374, 1.0708, 1.1043, 1.1378, 1.1712, 1.2047, 1.2382, 1.2716, 1.3051,
                                   1.3386, 1.3720, 1.4055, 1.4390, 1.4724, 1.5059, 1.5394, 1.5728, 1.6063, 1.6397, 1.6732],
                       'max': 1100,
                       'bins': xb_bins,
                       'bin_centers': xb_centers
                       }
    
    ## xb bins 
    for xx in range(1, len(xb_bins)):
        print('------------------------------------------------------------> creating asymettry plots for Xb bin {}'.format(xx))
        plot_slices( can, histos, 'xb_bin'+str(xx)+'_mass_bin{}' , output_pdfname, lab,
                     title='FD:Final, Inv.MassBins. {0:.3f}<Xb<{1:.3f}'.format(xb_bins[xx-1], xb_bins[xx]), xtitle='K^{ +}  K^{ -}  (GeV)', ytitle='counts' , plot_details=phi_mass_sliced)

        plot_asy( can, histos, 'trento_{}_xb_bin'+str(xx)+'_mass_bin{}', output_pdfname, lab,
                  title='Asy,{0:.3f}<Xb<{1:.3f}'.format(xb_bins[xx-1], xb_bins[xx]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=xb_im_asy_details, f_out='asy_vs_imkpkm_xbbin{0}_additionalcuts_binVar6_{1}.txt'.format(xx,field_setting))
     
    '''

    can.Print('{}]'.format(output_pdfname))
