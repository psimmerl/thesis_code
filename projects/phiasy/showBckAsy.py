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
                  TLine, kRed, kBlack, kBlue, kWhite, kGreen, kMagenta, kTRUE, gROOT)

from ROOT import( RooRealVar, RooArgSet, RooArgList, RooDataHist,
                  RooGaussian, RooFit, RooAddPdf, RooPolynomial, RooFormulaVar, RooLinkedList)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
global_dump = []

gROOT.SetBatch(kTRUE);

def readInBinningInfo(fname):
    fin = open(fname,'r')
    temp = []
    for ll in fin:
        print ll
        temp.append(float(ll))
    return temp    

def readInOverlapBinningInfo(fname):
    ovrlp_list = []
    fin = open(fname,'r')
    for ll in fin:
        temp_ll = ll.split(' ')
        print(temp_ll)
        ovrlp_list.append([float(temp_ll[0]), float(temp_ll[1])])
    return ovrlp_list


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
    
def plot_asy( can,  file_in, save_name, label, title, bin_info, field=None, fbin_out=None, bin_centers=None):

    dumpster=[]

    # txt file with bin center and asymetrry value with errors
    asy_per_mass = file_in
    x_cnt = array('d')
    y_val = array('d')
    x_err = array('d')
    y_err = array('d')

    x_cnt_sig = array('d')
    y_val_sig = array('d')
    x_err_sig = array('d')
    y_err_sig = array('d')

    chi2_values = []
    ndf_values = []
    x_all_values = []
    y_all_values = []

    cc=0
    weighted_avg = 0
    weighted_avg_den = 0
    x_cnt_sig_avg=0

    for ll in asy_per_mass:
        rr = ll.split(' ')
        print rr
        a=float(rr[0]) # mass
        b=float(rr[1])  # asy
        c=float(rr[2])  # asy error
        d=float(rr[3])  # chi2
        e=float(rr[4]) # NDF
        x_all_values.append(a)
        y_all_values.append(b)
        chi2_values.append(d)
        ndf_values.append(e)

        ## signal
        if cc > bin_info['min_sig_bin'] and cc < bin_info['max_sig_bin']: ## parameter to change here 
            x_cnt_sig.append(a)
            y_val_sig.append(b)
            x_err_sig.append(0.0)
            y_err_sig.append(c)
            ## calculate weighted average
            ## assuming sampling from IID , gaus
            
            weighted_avg+= ( b*c)
            weighted_avg_den+= (c)
            x_cnt_sig_avg+=a 

            print('weighted avg den {}'.format(weighted_avg_den))
            print('weighted avg {} '.format(weighted_avg/weighted_avg_den))
                                        
        else:
            x_cnt.append(a)
            y_val.append(b)
            x_err.append(0.0)
            y_err.append(c)
        cc+=1


    x_cnt_weight =  array('d')
    y_weight = array('d')
    x_cnt_weight_er = array('d')
    y_weight_err = array('d')

    # get weighted variance
    temp_weight_err = 0
    weight_err_den = 0


    print('--> weighted average is {}'.format(weighted_avg/weighted_avg_den))
    print('--> bin center is for weighted avg: {}'.format(x_cnt_sig_avg/len(x_cnt_sig)))
          
    for mm in range(0,len(y_err_sig)):
        print(' err {}'.format(y_err_sig[mm]))
        print(' 1/err^2 {}'.format(1.0/(y_err_sig[mm]**2)))
        temp_weight_err+= (1.0/(y_err_sig[mm]**2))

    print('--> sum of 1/sig2 {}'.format(temp_weight_err))
    weight_err = 1.0/np.sqrt(temp_weight_err)
    print('--> error {}'.format(y_err))

    asy_meas = weighted_avg/weighted_avg_den
    asy_meas_err = weight_err
    

    x_cnt_weight.append(x_cnt_sig_avg/len(x_cnt_sig)) ## avg bin center of phi peak
    y_weight.append(asy_meas)                         ## weighted asy
    x_cnt_weight_er.append(0.0)                       ## x error 
    y_weight_err.append(asy_meas_err)                 ## y error
    
    ## averaged measured asy to report
    avg_asy_graph = TGraphErrors(len(x_cnt_weight), x_cnt_weight, y_weight, x_cnt_weight_er , y_weight_err)
    avg_asy_graph.SetTitle('')
    avg_asy_graph.SetMarkerStyle(20)
    avg_asy_graph.SetMarkerSize(1.25)
    avg_asy_graph.SetLineColor(kGreen)
    avg_asy_graph.SetMarkerColor(kGreen)
    dumpster.append(avg_asy_graph)
        
    ## asy outside of signal
    print(' asymmetry outside signal range ')
    print(' x values {}'.format(x_cnt))
    print(' y values {}'.format(y_val))
    asy_per_mass_graph = None
    if 'truncate' in bin_info.keys():
        print(' truncating the graph to only include the left and right data point relevant for analysis')
        #x_cnt=x_cnt[:2]
        #y_val=y_val[:2]
        #x_err=x_err[:2]
        #y_err=y_err[:2]

        print(' truncated x {}'.format(x_cnt))
        print(' truncated x {}'.format(y_val))
        asy_per_mass_graph = TGraphErrors(len(x_cnt[:2]), x_cnt[:2], y_val[:2], x_err[:2], y_err[:2])
    else:
        asy_per_mass_graph = TGraphErrors(len(x_cnt), x_cnt, y_val, x_err, y_err)

    asy_per_mass_graph.SetTitle('')
    asy_per_mass_graph.SetMarkerStyle(21)
    asy_per_mass_graph.SetMarkerSize(1.5)
    asy_per_mass_graph.SetLineColor(kBlack)
    asy_per_mass_graph.SetMarkerColor(kBlack)
    dumpster.append(asy_per_mass_graph)

    ## specify whether to interpolate between the two side bands here with 
    ## bininfo
    background_fit = None 
    asy_bck = y_val[1]
    asy_bck_err = y_err[1]
        
    print(y_val)
    ## get the background fit to calculate the asy
    sidebands_to_draw=[]
    if 'sideband_ranges' in bin_info.keys():
        bb=1
        cc=0
        while bb < 4:
            print(bb)
            bb_mass_low= bin_info['mass_ranges'][bb-1]
            bb_mass_high = bin_info['mass_ranges'][bb]
            print(y_val)
            print(' plotting line for side band at {} from {} to {}'.format(y_val[bb-1], bb_mass_low, bb_mass_high))
            
            background_fit = TLine( bb_mass_low, y_val[cc], bb_mass_high, y_val[cc])            
            background_fit.SetLineColor(kBlue)            
            sidebands_to_draw.append(background_fit)
            dumpster.append(background_fit)
            bb+=2
            cc+=1
        ## change the asy_bck and err value now using both the left and right side band
        print('left asy {} +/- {}'.format(y_val[0],y_err[0]))
        print('right asy {} +/- {}'.format(y_val[1],y_err[1]))        
        asy_bck_err = np.sqrt(1/ (1/y_err[0]**2 + 1/y_err[1]**2) )
        asy_bck = ( y_val[0]/y_err[0]**2 + y_val[1]/y_err[1]**2 ) * asy_bck_err**2
        print(' calculated asy via sideband {} +/- {}'.format(asy_bck, asy_bck_err))
            
    else:
        background_fit = TF1('bck_asy',"[0]", 1.0268, 1.0542)# 1.0268, 1.1400)
        background_fit.SetLineColor(kBlue)
        asy_per_mass_graph.Fit("bck_asy","R")

            
    # asy of defined signal points (averaged and ploted in the avg_asy_graph)
    sig_asy_per_mass_graph = TGraphErrors(len(x_cnt_sig), x_cnt_sig, y_val_sig, x_err_sig, y_err_sig)
    sig_asy_per_mass_graph.SetTitle('')
    sig_asy_per_mass_graph.SetMarkerStyle(20)
    sig_asy_per_mass_graph.SetMarkerSize(2)    
    sig_asy_per_mass_graph.SetLineColor(kRed)
    sig_asy_per_mass_graph.SetMarkerColor(kRed)
    dumpster.append(sig_asy_per_mass_graph)

    ### correct the weighted signal value ( or in the case of one bin, that bins value )
    ## the correct asy we measure if calculated as 
    ## A_s = A_m + B/S * (A_m - A_b)
    B = 1
    S = 1 
    BdivS = 1.0/bin_info['s_to_b_ratio']
    asy_signal = asy_meas + (BdivS) * (asy_meas - asy_bck ) # y_val[1] is the bin next to the right of the signal for sideband subtraction
    print(' - ------ --- ------ > bck asy is {}'.format(asy_bck) )
    asy_signal_err = np.sqrt( asy_meas_err**2 * ( (1+(BdivS))**2 ) + asy_bck_err**2 * ( (BdivS)** 2) ) ## y_err is the error to bsa to the right of signal
    x_cnt_corr = array('d')
    y_val_corr = array('d')
    x_cnt_corr_err = array('d')
    y_cnt_corr_err = array('d')
    x_cnt_corr.append(x_cnt_sig_avg/len(x_cnt_sig) - 0.001) ## artifically shift to view
    y_val_corr.append(asy_signal)
    x_cnt_corr_err.append(0.0)
    y_cnt_corr_err.append( asy_signal_err )
    corrected_asy_graph = TGraphErrors(len(x_cnt_corr), x_cnt_corr, y_val_corr, x_cnt_corr_err, y_cnt_corr_err)
    corrected_asy_graph.SetTitle('')
    corrected_asy_graph.SetMarkerStyle(21)
    corrected_asy_graph.SetMarkerSize(2)
    corrected_asy_graph.SetLineColor(kMagenta)
    corrected_asy_graph.SetMarkerColor(kMagenta)
    dumpster.append(corrected_asy_graph)
    if fbin_out:
        fbin_out.write(str(bin_centers) + ' ' + str(asy_signal) + ' ' + str(asy_signal_err) + '\n')

    ## all plots together!
    mg_asy_per_mass = TMultiGraph()
    mg_asy_per_mass.SetTitle('')
    mg_asy_per_mass.Add(asy_per_mass_graph)
    #mg_asy_per_mass.Add(sig_asy_per_mass_graph) ## here
    mg_asy_per_mass.Add(avg_asy_graph)
    mg_asy_per_mass.Add(corrected_asy_graph)
    dumpster.append(mg_asy_per_mass)

    
    mg_asy_per_mass.GetHistogram().SetMaximum(0.5)
    mg_asy_per_mass.GetHistogram().SetMinimum(-0.5)
    mg_asy_per_mass.GetXaxis().SetLimits(0.95, 1.15)# 1.35) #4) # chaged for poster contest
    dumpster.append(mg_asy_per_mass)
    
    # draw grid
    h = gPad.DrawFrame(1.4, -0.5, 0.95, 0.5)
    h.GetXaxis().SetNdivisions(-9)
    h.GetYaxis().SetNdivisions(-9)
    gPad.SetGrid()
    dumpster.append(h)

    dumpster.append(background_fit)
    
    mg_asy_per_mass.Draw('AP')
    
    dumpster.append(mg_asy_per_mass)
    #background_fit.Draw('same')
    for sb in sidebands_to_draw:
        sb.SetLineColor(kBlue)
        sb.SetLineWidth(2)
        sb.Draw('same')
        dumpster.append(sb)
    
    if 'additional_sidebands' in bin_info.keys():
        sbb_counter=2
        shift_y = 0
        for sbb in bin_info['additional_sidebands']:
            print(sbb)
            minx = sbb[0]
            maxx = sbb[1]
            ## only care about the overlap with right side band
            print(' drawing overlap line from {} to {} at {}'.format(minx, maxx, y_val[sbb_counter]))
            if sbb_counter == 2: 
                shift_y=0.02
                print(' pushing range bar down by 0.02')
            sbb_additional = TLine(minx, y_val[sbb_counter]-shift_y, maxx, y_val[sbb_counter]-shift_y)
            sbb_additional.SetLineColor(kBlue+3)
            sbb_additional.SetLineWidth(2)
            sbb_additional.Draw('same')
            dumpster.append(sbb_additional)
            sbb_counter+=1

    

    #label.DrawLatex(0.1, 0.925, 'BSA Per I.M. of K^{ +}K^{ -} Slice')        
    if title:
        label.DrawLatex(0.1, 0.925, title)        
    label.DrawLatex(0.5, 0.015,'I.M. of K^{ +}K^{ -} (GeV)')            
    label.SetTextAngle(90)
    label.DrawLatex(0.035, 0.5, 'Beam Spin Asy.')
    label.SetTextAngle(0)
    dumpster.append(label)

    ## draw the BSA value of signal on graph
    temp_lab = TLatex()
    temp_lab.SetNDC()
    temp_lab.SetTextFont(42)
    temp_lab.SetTextSize(0.08)
    temp_lab.SetTextColor(1)

    #temp_label.DrawLatex(0.2, 0.25, 'A_{s}'+':{0:.5f}'.format(asy_signal))
    #temp_label.DrawLatex(0.2, 0.20, '#Delta A_{s}'+':{0:.5f}'.format(asy_signal_err))
    temp_lab.DrawLatex(0.23, 0.63, 'A_{s}'+':{0:.5f}'.format(asy_signal))
    temp_lab.DrawLatex(0.23, 0.55, '#Delta A_{s}'+':{0:.5f}'.format(asy_signal_err))

    dumpster.append(temp_lab)
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)

    # dont draw for poster contest 
    '''
    for cc in range(0, len(x_all_values)):
        label.SetTextAngle(60)
        label.DrawLatex((x_all_values[cc]-0.9)/(1.4-0.9), np.abs(y_all_values[cc])/1 + 0.5 , '{0:.2f}/{1}'.format(chi2_values[cc], ndf_values[cc]))
        label.SetTextAngle(0)
        dumpster.append(label)
    '''

    can.Print(output_pdfname)


def plot_binned_asy( can, f_in, save_name, label, title, bin_info, xtitle=None, ytitle=None, field=None):
    dumpster=[]

    x_values =  array('d')
    y_values = array('d')
    x_err = array('d')
    y_err = array('d')

    x_values_corr =  array('d')
    y_values_corr = array('d')
    x_err_corr = array('d')
    y_err_corr = array('d')

    # get weighted variance
    weight_err = 0
    weight_err_den = 0
    print(' starting asy for kinematics ' )
    print(' kinematic bin range from {} to {} '.format(bin_info['min_bin'],bin_info['max_bin']))
    ii=0
    for ff in range(bin_info['min_bin'],bin_info['max_bin']):
        f_asy = open(f_in.format(ff))
        cc = 0
        weighted_avg = 0
        weighted_avg_den = 0
        x_cnt_sig_avg=0
        n_sig_bins=0
        # loop over entries in txt file (mass bin center , asy, asy error)
        # define bin information to assign asymmetry measurement to

        kinematic_bin_min = bin_info['kinematic_bin'][ff-1]
        kinematic_bin_max = bin_info['kinematic_bin'][ff]
        kinematic_bin_avg = bin_info['bin_centers'][ff-1]
        kinematic_bin_err = kinematic_bin_max - kinematic_bin_min

        print(' Getting Asymmetry in range {} to {}'.format(kinematic_bin_min, kinematic_bin_max) )
        print(' Average Value of bin range is {}'.format(kinematic_bin_avg))
        y_results = []
        y_error = []
        bck = []
        bck_err = []
        for ll in f_asy:
            rr = ll.split(' ')
            print rr
            a= float(rr[0]) # mass bin center
            b=float(rr[1]) # asy 
            c=float(rr[2]) # asy err
            d=float(rr[3])  # chi2
            e=float(rr[4]) # NDF

            ## signal - only use those entries in the signal region for this part
            if cc > bin_info['min_sig_bin'] and cc < bin_info['max_sig_bin']: ## parameter to change here 
                ## calculate most prob weighted average
                ## assuming sampling from IID , gaus
            
                weighted_avg+= ( b*c)
                weighted_avg_den+= (c)
                x_cnt_sig_avg+=a 
                n_sig_bins+=1
                print('weighted avg den {}'.format(weighted_avg_den))
                print('weighted avg {} '.format(weighted_avg/weighted_avg_den))
                
                y_results.append(b)
                y_error.append(c)
            else:
                bck.append(b)
                bck_err.append(c)

            cc+=1

            
        print('--> weighted average is {}'.format(weighted_avg/weighted_avg_den))
        print('--> weighted mass bin center is {}'.format(x_cnt_sig_avg/(n_sig_bins)))
        
        print('--> error values of signal {}'.format(y_error))
        weight_err=0
        for mm in range(0,n_sig_bins):
            weight_err+=(1.0/y_error[mm]**2)
            print(' --> error is {}'.format(y_error[mm]))
            print(' --> error calculation: 1/sig = {}'.format(1.0/y_error[mm]))
            print(' --> error calculation: 1/sig^2 = {}'.format(1.0/y_error[mm]**2))
        
        print('---> final 1/error^2={}'.format(weight_err))
        print('---> final error = {}'.format(1.0/np.sqrt(weight_err)))
        
        asy_meas = weighted_avg/weighted_avg_den
        asy_meas_err = 1.0/np.sqrt(weight_err)
        
        
        ## now fill for each kinematic bin
        x_values.append(kinematic_bin_avg)
        y_values.append(asy_meas)
        x_err.append(0.0)        
        y_err.append(asy_meas_err)

        ## now fill for each kinematic bin with the sideband subtracted from the asy
        x_values_corr.append(kinematic_bin_avg - 0.05)
        y_values_corr.append(asy_signal)
        x_err_corr.append(0.0)        
        y_err_corr.append(asy_signal_err)
    

    print(x_values)
    print(y_values)
    print(x_err)
    print(y_err)
    
    kin_asy_graph = TGraphErrors(len(x_values), x_values, y_values, x_err , y_err)
    kin_asy_graph.SetTitle('')
    kin_asy_graph.SetMarkerStyle(29)
    kin_asy_graph.SetMarkerSize(1)
    kin_asy_graph.SetLineColor(kRed)
    kin_asy_graph.SetMarkerColor(kRed)

    kin_asy_corr_graph = TGraphErrors(len(x_values_corr), x_values_corr, y_values_corr, x_err_corr , y_err_corr)
    kin_asy_corr_graph.SetTitle('')
    kin_asy_corr_graph.SetMarkerStyle(28)
    kin_asy_corr_graph.SetMarkerSize(1)
    kin_asy_corr_graph.SetLineColor(kMagenta)
    kin_asy_corr_graph.SetMarkerColor(kMagenta)
    
    mg_kin_asy_graph = TMultiGraph()
    mg_kin_asy_graph.Add(kin_asy_graph)
    mg_kin_asy_graph.Add(kin_asy_corr_graph)

    mg_kin_asy_graph.GetHistogram().SetMaximum(0.5)
    mg_kin_asy_graph.GetHistogram().SetMinimum(-0.5)
    mg_kin_asy_graph.GetXaxis().SetLimits(bin_info['xlim'][0], bin_info['xlim'][1])
        
    dumpster.append(kin_asy_graph)
    dumpster.append(kin_asy_corr_graph)
    dumpster.append(mg_kin_asy_graph)

                      
    # draw grid
    h = gPad.DrawFrame(1.4, -0.5, 0.95, 0.5)
    h.GetXaxis().SetNdivisions(-9)
    h.GetYaxis().SetNdivisions(-9)
    #h.GetXaxis().SetLabelSize(0.025)
    #h.GetYaxis().SetLabelSize(0.025)
    gPad.SetGrid()
    dumpster.append(h)

    mg_kin_asy_graph.Draw('AP')

    #### print out numbers next to data points
    ## print out asy corrected for sideband contribution
    for ii in range(0,len(x_values)):
        print(' draw info fraction here {}'.format(x_values_corr[ii]/bin_info['xlim'][1]))
        label.DrawLatex(x_values[ii]/bin_info['xlim'][1], 0.64,'<{}>:{:3.3f}'.format(bin_info['name'],x_values[ii]))
        label.DrawLatex(x_values[ii]/bin_info['xlim'][1], 0.59,'Asy:{:.4f}'.format(y_values_corr[ii]))
        label.DrawLatex(x_values[ii]/bin_info['xlim'][1], 0.54,'Err:{:.4f}'.format(y_err_corr[ii]))


    if title:
        label.DrawLatex(0.1, 0.925, title)        
    if xtitle:
        label.DrawLatex(0.65, 0.015, xtitle)            
    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.035, 0.5, ytitle)
        label.SetTextAngle(0)
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)


    dumpster.append(label)

    can.Print(output_pdfname)


def plot_binned_asyV2( can, f_in, save_name, label, title, bin_info, xtitle=None, ytitle=None, field=None):
    dumpster=[]

    x_values =  array('d')
    y_values = array('d')
    x_err = array('d')
    y_err = array('d')

    x_values_corr =  array('d')
    y_values_corr = array('d')
    x_err_corr = array('d')
    y_err_corr = array('d')

    # get weighted variance
    weight_err = 0
    weight_err_den = 0
    print(' starting asy for kinematics ' )
    print(' kinematic bin range from {} to {} '.format(bin_info['min_bin'],bin_info['max_bin']))
    ii=0
    print(f_in)
    #fin = open(f_in,'r')
    
    print('opening file ')
    for ll in f_in:
        bi = ll.split(' ')
        print(bi)
        bin_center = float(bi[0])
        bin_y_val = float(bi[1])
        bin_y_err = float(bi[2])
        
        #corrected means taking into account left and right bins and the s/b ratio
        x_values_corr.append(bin_center)
        y_values_corr.append(bin_y_val)
        y_err_corr.append(bin_y_err)
        x_err_corr.append(0.0)
        print(bin_center)
    print(x_values_corr)
    kin_asy_corr_graph = TGraphErrors(len(x_values_corr), x_values_corr, y_values_corr, x_err_corr , y_err_corr)
    kin_asy_corr_graph.SetTitle('')
    kin_asy_corr_graph.SetMarkerStyle(28)
    kin_asy_corr_graph.SetMarkerSize(1)
    kin_asy_corr_graph.SetLineColor(kMagenta)
    kin_asy_corr_graph.SetMarkerColor(kMagenta)
    
    mg_kin_asy_graph = TMultiGraph()
    #mg_kin_asy_graph.Add(kin_asy_graph)
    mg_kin_asy_graph.Add(kin_asy_corr_graph)

    mg_kin_asy_graph.GetHistogram().SetMaximum(0.5)
    mg_kin_asy_graph.GetHistogram().SetMinimum(-0.5)
    mg_kin_asy_graph.GetXaxis().SetLimits(bin_info['xlim'][0], bin_info['xlim'][1])
        
    #dumpster.append(kin_asy_graph)
    dumpster.append(kin_asy_corr_graph)
    dumpster.append(mg_kin_asy_graph)

                      
    # draw grid
    h = gPad.DrawFrame(1.4, -0.5, 0.95, 0.5)
    h.GetXaxis().SetNdivisions(-9)
    h.GetYaxis().SetNdivisions(-9)
    #h.GetXaxis().SetLabelSize(0.025)
    #h.GetYaxis().SetLabelSize(0.025)
    gPad.SetGrid()
    dumpster.append(h)

    mg_kin_asy_graph.Draw('AP')

    #### print out numbers next to data points
    ## print out asy corrected for sideband contribution
    for ii in range(0,len(x_values_corr)):
        if 'plot_title_center_shift' in bin_info.keys():
            print(' draw info fraction here {}'.format(x_values_corr[ii]/bin_info['xlim'][1]))
            x_shift = bin_info['plot_title_center_shift'][ii]
            label.DrawLatex(x_shift + 0.05 + x_values_corr[ii]/bin_info['xlim'][1], 0.64,'<{}>:{:3.3f}'.format(bin_info['name'],x_values_corr[ii]))
            label.DrawLatex(x_shift + 0.05 + x_values_corr[ii]/bin_info['xlim'][1], 0.59,'Asy:{:.4f}'.format(y_values_corr[ii]))
            label.DrawLatex(x_shift + 0.05 + x_values_corr[ii]/bin_info['xlim'][1], 0.54,'Err:{:.4f}'.format(y_err_corr[ii]))
        else:
            print(' draw info fraction here {}'.format(x_values_corr[ii]/bin_info['xlim'][1]))
            label.DrawLatex(0.05 + x_values_corr[ii]/bin_info['xlim'][1], 0.64,'<{}>:{:3.3f}'.format(bin_info['name'],x_values_corr[ii]))
            label.DrawLatex(0.05 + x_values_corr[ii]/bin_info['xlim'][1], 0.59,'Asy:{:.4f}'.format(y_values_corr[ii]))
            label.DrawLatex(0.05 + x_values_corr[ii]/bin_info['xlim'][1], 0.54,'Err:{:.4f}'.format(y_err_corr[ii]))


    if title:
        label.DrawLatex(0.1, 0.925, title)        
    if xtitle:
        label.DrawLatex(0.65, 0.015, xtitle)            
    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.035, 0.5, ytitle)
        label.SetTextAngle(0)
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)


    dumpster.append(label)

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
    field_setting='inbNoutb'
    binvar='Var11'
    #asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_binVar11_inbNoutb.txt
    asy_per_mass = open('asy_vs_imkpkm_pidtype4a_additionalcuts_rmvres12_bin'+binvar+'_'+field_setting+'V2.txt','r')

    ## integrated s/b inb is 0.7959 for range of 1.0126 - 1.0268
    ## integrated s/b out is 1.1237

    ## newer results used are
    ## integrated s/b inb is 0.92640 for range of 1.0120 - 1.0274
    ## integrated s/b out is 1.20689
    sb_ratio = 0
    # list of s/b values for outb
    # 1.5, 2, 2.5, 3
    #sb_ratio_list = [1.3518, 1.1197, 0.9298, 0.7856]
    
    #list of s/b values for inb
    # 1.5, 2, 2.5, 3
    #sb_ratio_list = [1.0030, 0.8296, 0.6899, 0.5829] ## inbending 

    ## list of s/b values for inb n out
    # 1.5, 2, 2.5, 3
    #sb_ratio_list = [1.2180, 1.0088, 0.8377, 0.7077] ## inb + outb this is for without removing resonances
    sb_ratio_list = [1.4479, 1.1992, 0.99571, 0.8411]
    
    # 0.8809 rmv res1
    # 0.8251 rmv res2

    ## outbending remove res s/b ratio
    ## rmv res1 1.21481
    ## rmv res2 1.14334
    ## rmv res12 1.3538

    ## inb and outb remove res s/b ratio
    ## rmv res1 1.08765
    ## rmv res2 1.01143
    ## rmv res12 1.16486
    
    # pidtype 3 is 1.0287
    # pid type 4 is 1.0392
    # pid type 4a  is 1.0464 -> correct cal SF parameterization 

    if( field_setting == 'inb' ):
        sb_ratio = 0.89078 # 0.8296 with res12 #sb_ratio_list[3]
    elif( field_setting == 'outb'):
        sb_ratio =  1.3538 # new nominal for outb #nominal sb ratio 1.1197 with res12
    elif( field_setting == 'inbNoutb'):
        sb_ratio =  1.0464 #1.0287 # 1.06215 #<- this is for eb pid only # nominal sb ratio with res12 0.97962

    mass_bins = readInBinningInfo('bin_mass_ranges'+binvar+'_'+field_setting+'.txt')
    overlap_mass_bins = readInOverlapBinningInfo('bin_mass_slidingVar4_'+field_setting+'.txt')

    bin_info = {'min_sig_bin':0,
                'max_sig_bin':2,
                's_to_b_ratio': sb_ratio,
                'sideband_ranges': [0, 1, 3, 4],
                'mass_ranges' : mass_bins,
                'truncate' : True
            }

    
    q2_bins=readInBinningInfo('bin_limits_q2_'+field_setting+'.txt')
    xb_bins=readInBinningInfo('bin_limits_xb_'+field_setting+'.txt')
    t_bins=readInBinningInfo('bin_limits_t_'+field_setting+'.txt')
    w_bins=readInBinningInfo('bin_limits_w_'+field_setting+'.txt')
    
    q2_centers = readInBinningInfo('bin_centers_q2_'+field_setting+'.txt')
    xb_centers = readInBinningInfo('bin_centers_xb_'+field_setting+'.txt')
    t_centers = readInBinningInfo('bin_centers_t_'+field_setting+'.txt')
    w_centers = readInBinningInfo('bin_centers_w_'+field_setting+'.txt')

    # signal to ratio info from the mon-phimass.py file
    q2_bins_sb_ratio = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/exclusive_phi/epkpkm_top/q2_bin_mass_sb_ratio_'+field_setting+'_final.txt')
    t_bins_sb_ratio = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/exclusive_phi/epkpkm_top/t_bin_mass_sb_ratio_'+field_setting+'_final.txt')
    xb_bins_sb_ratio = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/exclusive_phi/epkpkm_top/xb_bin_mass_sb_ratio_'+field_setting+'_final.txt')
    w_bins_sb_ratio = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/exclusive_phi/epkpkm_top/w_bin_mass_sb_ratio_'+field_setting+'_final.txt')
    #w_bins_sb_ratio_lowq2 = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/exclusive_phi/epkpkm_top/w_bin_mass_lowq2_sb_ratio_'+field_setting+'_final.txt')

    
    ## get asy
    plot_asy( can, asy_per_mass, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice, Remove #Lambda(1520), #Lambda(1820)', bin_info, field=field_setting)
    '''
    
    # get asy from overlap
    
    asy_per_mass_sliding =open('asy_vs_imkpkm_pidtype1_additionalcuts_sliding_Var4_'+field_setting+'.txt','r')

    bin_info_slide = {'min_sig_bin':0,
                      'max_sig_bin':2,
                      's_to_b_ratio': sb_ratio,
                      'sideband_ranges': [0,1, 3,4],
                      'mass_ranges' : mass_bins,            
                      'additional_sidebands' : [[overlap_mass_bins[3][0], overlap_mass_bins[3][1]]]
            }

    
    plot_asy( can, asy_per_mass_sliding, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice Overlap', bin_info_slide, field=field_setting)
    
    
    
    
    
    
    # check over -t, but integrated over Q2, xb
    # V2.txt to differentiate from results first time around which did use correct files
    t_bin_out = open('asy_vs_tbins_additionalcuts_bin{0}_{1}V2.txt'.format(binvar, field_setting),'w')
    for tt in range(1,len(t_bins)):
        asy_per_mass_tbin = open('asy_vs_imkpkm_tbin{0}_additionalcuts_bin{2}_{1}V2.txt'.format(tt,field_setting,binvar),'r')
        plot_asy( can, asy_per_mass_tbin, output_pdfname, lab, 'BSA Per I.M. of K^{ +} K^{ -} Slice'+', {0:.4f}<-t<{1:.4f}'.format(t_bins[tt-1],t_bins[tt]), 
                  bin_info, field=field_setting, fbin_out = t_bin_out, bin_centers=t_centers[tt-1] )
    t_bin_out.close()
        
    
    # now plot averages asy as a function of -t 
    t_bin_info={ 'name':'-t',
                 'kinematic_bin' : t_bins,
                 'min_bin': 1,
                 'max_bin': len(t_bins),
                 'xlim':[0.0, 7],
                 'min_sig_bin': 0,
                 'max_sig_bin': 2,
                 'bin_centers': t_centers,
                 's_to_b_ratio': t_bins_sb_ratio
             }

    t_binned_asy = open('asy_vs_tbins_additionalcuts_bin'+binvar+'_'+field_setting+'V2.txt','r')

    plot_binned_asyV2( can, t_binned_asy, output_pdfname, lab, 'BSA vs -t via Slices in I.M. of K^{ +}K^{ -}' , t_bin_info, 
                     xtitle='-t (GeV)', ytitle='A^{sin(#phi)}_{LU}', field=field_setting)
    
    
    
    
    # check over Q2, but integrated over xb, -t
    q2_bin_out = open('asy_vs_q2bins_additionalcuts_bin{0}_{1}V2.txt'.format(binvar, field_setting),'w')
    for qq in range(1,len(q2_bins)):
        asy_per_mass_q2bin = open('asy_vs_imkpkm_q2bin{0}_additionalcuts_bin{2}_{1}V2.txt'.format(qq, field_setting,binvar),'r')
        plot_asy( can, asy_per_mass_q2bin, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice,'+' {0:.4f}< Q2 <{1:.4f}'.format(q2_bins[qq-1],q2_bins[qq]), 
                  bin_info, field=field_setting, fbin_out = q2_bin_out, bin_centers = q2_centers[qq-1])
    q2_bin_out.close()

    # now plot averages asy as a function of Q2
    q2_bin_info={ 'name':'Q^{2}',
                  'kinematic_bin' : q2_bins,
                  'min_bin': 1,
                  'max_bin': len(q2_bins),
                  'xlim':[0.0, 10],
                  'min_sig_bin': 0,
                  'max_sig_bin': 2,
                  'bin_centers': q2_centers,
                  's_to_b_ratio': q2_bins_sb_ratio
             }
    #plot_binned_asy( can, 'asy_vs_imkpkm_q2bin{}_additionalcuts_'+binvar+'_'+field_setting+'.txt', output_pdfname, lab, 'BSA vs Q^{2} via Slices in I.M. of K^{ +}K^{ -}' , q2_bin_info, xtitle='Q^{2} (GeV)',ytitle='A^{sin(#phi)}_{LU}', field=field_setting)
    q2_binned_asy = open('asy_vs_q2bins_additionalcuts_bin'+binvar+'_'+field_setting+'V2.txt','r')

    plot_binned_asyV2( can, q2_binned_asy, output_pdfname, lab, 'BSA vs Q^{2} via Slices in I.M. of K^{ +}K^{ -}' , q2_bin_info, 
                       xtitle='Q^{2} (GeV^{2})', ytitle='A^{sin(#phi)}_{LU}', field=field_setting)



    # check over xb, but integrated over Q2, -t
    xb_bin_out = open('asy_vs_xbbins_additionalcuts_bin{0}_{1}V2.txt'.format(binvar, field_setting),'w')    
    for xx in range(1,len(xb_bins)):
        asy_per_mass_xbbin = open('asy_vs_imkpkm_xbbin{0}_additionalcuts_bin{2}_{1}V2.txt'.format(xx, field_setting,binvar),'r')
        plot_asy( can, asy_per_mass_xbbin, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice,'+' {0:.4f}< Xb <{1:.4f}'.format(xb_bins[xx-1],xb_bins[xx]),
                  bin_info, field=field_setting, fbin_out = xb_bin_out, bin_centers = xb_centers[xx-1])
    xb_bin_out.close()
    
    # now plot averages asy as a function of Xb
    xb_bin_info={ 'name':'X_{b}',
                  'kinematic_bin' : xb_bins,
                  'min_bin': 1,
                  'max_bin': len(xb_bins),
                  'xlim':[0.0, 1],
                  'min_sig_bin': 0,
                  'max_sig_bin': 2,
                  'bin_centers': xb_centers,
                  's_to_b_ratio': [1,1]
             }
    #plot_binned_asy( can, 'asy_vs_imkpkm_xbbin{}_additionalcuts_'+binvar+'_'+field_setting+'.txt', output_pdfname, lab, 'BSA vs X_{b} via Slices in I.M. of K^{ +}K^{ -}' , xb_bin_info, xtitle='X_{b}',ytitle='A^{sin(#phi)}_{LU}', field=field_setting)
    xb_binned_asy = open('asy_vs_xbbins_additionalcuts_bin'+binvar+'_'+field_setting+'V2.txt','r')
    plot_binned_asyV2( can, xb_binned_asy, output_pdfname, lab, 'BSA vs X_{b} via Slices in I.M. of K^{ +}K^{ -}' , xb_bin_info, 
                       xtitle='X_{b}', ytitle='A^{sin(#phi)}_{LU}', field = field_setting)
    
    
    ######################################
    ## multi dim binning low xb, q2 
    bin_info_xbq2 = {'min_sig_bin':0,
                     'max_sig_bin':2,
                     's_to_b_ratio': 1.4655,
                     'sideband_ranges': [0, 1, 3, 4],
                     'mass_ranges' : mass_bins,
                     'truncate' : True
            }

    asy_per_mass_xbbin1_q2bin1 = open('asy_vs_imkpkm_pidtype1_additionalcuts_xb_bin1_q2_bin1_binVar10_inbNoutbV2.txt','r')
    plot_asy( can, asy_per_mass_xbbin1_q2bin1, output_pdfname, lab, 'BSA Per I.M. of K^{ +} K^{ -} Slice, low X_{b} vs Q^{2}',
              bin_info_xbq2, field=field_setting)

    ######################################
    ## multi dim binning high xb, q2 
    bin_info_xbq2_high = {'min_sig_bin':0,
                          'max_sig_bin':2,
                          's_to_b_ratio': 0.64658,
                          'sideband_ranges': [0, 1, 3, 4],
                          'mass_ranges' : mass_bins,
                          'truncate' : True
                      }

    asy_per_mass_xbbin2_q2bin2 = open('asy_vs_imkpkm_pidtype1_additionalcuts_xb_bin2_q2_bin2_binVar10_inbNoutb.txt','r')
    plot_asy( can, asy_per_mass_xbbin2_q2bin2, output_pdfname, lab, 'BSA Per I.M. of K^{ +} K^{ -} Slice, high X_{b} vs Q^{2}',
              bin_info_xbq2_high, field=field_setting)
    
    

    #################################################3
    # check over W, but integrated over Q2, xb, -t 
    ## v2 for txt file denotes that I this is different than what was shown previously to 08/06/2020 because the binning changed a little bit. 
    w_bin_out = open('asy_vs_wbins_additionalcuts_bin{0}_{1}V2.txt'.format(binvar, field_setting),'w')
    for ww in range(1,len(w_bins)):
        asy_per_mass_wbin = open('asy_vs_imkpkm_wbin{0}_additionalcuts_bin{2}_{1}V2.txt'.format(ww,field_setting,binvar),'r')
        plot_asy( can, asy_per_mass_wbin, output_pdfname, lab, 'BSA Per I.M. of K^{ +} K^{ -} Slice'+', {0:.4f}<W<{1:.4f}'.format(w_bins[ww-1], w_bins[ww]), 
                  bin_info, field=field_setting, fbin_out = w_bin_out, bin_centers=w_centers[ww-1] )
    w_bin_out.close()
            
    # now plot averages asy as a function of W
    w_bin_info={ 'name':'W',
                 'kinematic_bin' : w_bins,
                 'min_bin': 1,
                 'max_bin': len(w_bins),
                 'xlim':[0.0, 5],
                 'min_sig_bin': 0,
                 'max_sig_bin': 2,
                 'bin_centers': w_centers,
                 's_to_b_ratio': w_bins_sb_ratio,
                 'plot_title_center_shift': [-0.15, -0.05]
             }


    w_binned_asy = open('asy_vs_wbins_additionalcuts_bin'+binvar+'_'+field_setting+'V2.txt','r')

    plot_binned_asyV2( can, w_binned_asy, output_pdfname, lab, 'BSA vs W via Slices in I.M. of K^{ +}K^{ -}' , w_bin_info, 
                     xtitle='W (GeV)', ytitle='A^{sin(#phi)}_{LU}', field=field_setting)

    
    ##
    ###########
    ## W integrated overall variables but cut on high q2 , so q2 < 2.75
    w_bin_out_lowq2 = open('asy_vs_wbins_lowq2_additionalcuts_bin{0}_{1}.txt'.format(binvar, field_setting),'w')
    for ww in range(1,len(w_bins)):
        asy_per_mass_wbin_lowq2 = open('asy_vs_imkpkm_wbin{0}_lowq2_additionalcuts_bin{2}_{1}.txt'.format(ww,field_setting,binvar),'r')
        plot_asy( can, asy_per_mass_wbin_lowq2, output_pdfname, lab, 'BSA Per I.M. of K^{ +} K^{ -} Slice, Q^{2}<2.75'+', {0:.4f}<W<{1:.4f}'.format(w_bins[ww-1], w_bins[ww]), 
                  bin_info, field=field_setting, fbin_out = w_bin_out_lowq2, bin_centers = w_centers[ww-1] )
    w_bin_out_lowq2.close()
            
    # now plot averages asy as a function of W
    w_bin_info_lowq2={ 'name':'W',
                       'kinematic_bin' : w_bins,
                       'min_bin': 1,
                       'max_bin': len(w_bins),
                       'xlim':[0.0, 5],
                       'min_sig_bin': 0,
                       'max_sig_bin': 2,
                       'bin_centers': w_centers,
                       's_to_b_ratio': w_bins_sb_ratio_lowq2,
                       'plot_title_center_shift': [-0.03, 0.03]
             }


    w_binned_asy_lowq2 = open('asy_vs_wbins_lowq2_additionalcuts_bin'+binvar+'_'+field_setting+'.txt','r')
    print(w_binned_asy_lowq2)
    plot_binned_asyV2( can, w_binned_asy_lowq2, output_pdfname, lab, 'BSA vs W via Slices in I.M. of K^{ +}K^{ -}, Q^{2}<2.75' , w_bin_info_lowq2, 
                     xtitle='W (GeV)', ytitle='A^{sin(#phi)}_{LU}', field=field_setting)
    '''

    '''
    #############
    ### print the asymmetry for different scenarios
    # 1 - a_lu^sinphi from sinx / 1 + cosx
    # 2 - a_lu^sinphi from sinx / 1 + cos2x
    # 3 - bp calculated using m1
    # 4 - bp calculated using m2 ( not method 2 or 1 above, but how to determine averaged bp over inb + outb data)
    
    #1
    asy_per_mass = open('asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_bin'+binvar+'_'+field_setting+'_systmodbsa_V2.txt','r')
    plot_asy( can, asy_per_mass, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice, Remove #Lambda(1520), #Lambda(1820), sin(#phi)/(1+cos(#phi))', bin_info, field=field_setting)

    #2
    asy_per_mass = open('asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_bin'+binvar+'_'+field_setting+'_systmodbsa2_V2.txt','r')
    plot_asy( can, asy_per_mass, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice, Remove #Lambda(1520), #Lambda(1820), sin(#phi)/(1+cos(2#phi))', bin_info, field=field_setting)

    #3
    asy_per_mass = open('asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_bin'+binvar+'_'+field_setting+'V2_bpm1.txt','r')
    plot_asy( can, asy_per_mass, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice, Remove #Lambda(1520), #Lambda(1820), BP 87.9%', bin_info, field=field_setting)

    #4
    asy_per_mass = open('asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_bin'+binvar+'_'+field_setting+'V2_bpm2.txt','r')
    plot_asy( can, asy_per_mass, output_pdfname, lab, 'BSA Per I.M. of K^{ +}K^{ -} Slice, Remove #Lambda(1520), #Lambda(1820), BP 88.5%', bin_info, field=field_setting)
    '''

    can.Print('{}]'.format(output_pdfname))
