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
                  gPad, gStyle, TLatex, TGraphErrors, TMath,
                  TLine, kRed, kBlack, kBlue, kWhite,  kTRUE, gROOT, TBox)

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

def get_fit_limits(txt_file_name, sig_range):
    params = {}
    txt_file = open(txt_file_name,'r')
    for ll in txt_file:
        sect = int(ll.split(' ')[0])
        mean = float(ll.split(' ')[1])
        sig = float(ll.split(' ')[2])
        params[sect] = [mean+sig_range*sig, mean-sig_range*sig]
    return params
    
def getFitFractionalHeight( h_temp, h_name, percent_max, save_fits):
    # fitting routine adopted from Andrew Puckett.
    print(' special fitting routine ')
    xlow, xhigh, histmax = 0, 0, 0
    binlow, binhigh, binmax = 0, 0, 0

    
    binmax = h_temp.GetMaximumBin()
    histmax = h_temp.GetMaximum()#GetBinContent(h_temp.GetMaximumBin())
    #print('histmax ' + str(histmax))
    
    binlow = binmax
    binhigh = binmax

    while h_temp.GetBinContent(binhigh) >= percent_max*histmax and binhigh <= h_temp.GetNbinsX() : binhigh+=1
    while h_temp.GetBinContent(binlow) >= percent_max*histmax and binlow > 1 : binlow-=1
    
    xlow = h_temp.GetBinLowEdge(binlow)
    xhigh = h_temp.GetBinLowEdge(binhigh+1)
    
    print(h_name)
    print(' >> bin high values ' + str(binhigh) + ' bin low ' + str(binlow) )
    print(" >> values used " + str(xlow) + " " + str(xhigh) + " " + str(histmax) )
    
    fit_temp = TF1(h_name,'gaus', xlow, xhigh )
    fit_temp.SetParameter(0, histmax)
    fit_temp.SetParameter(1, h_temp.GetBinCenter(h_temp.GetMaximumBin()))
    
                              
    h_temp.Fit(h_name,"R") # add Q for quiet
    save_fits.append(fit_temp)
    
    temp_mean = fit_temp.GetParameter(1)
    temp_rms = fit_temp.GetParameter(2)
    
    return fit_temp#[temp_mean, temp_rms ]


def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None,  log=False,
                     y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None,
                     hline=None, hlines = None ,special_plot = None, 
                     fit_par = None, fit_range=None, sector_fit_info=None, fract_fit=None, save_fit_results=None, vlines=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    fit_result_out = None
    if save_fit_results:
        fit_result_out = open( save_fit_results, 'w')


    for i in range(1,7):
        canvas.cd(i)

        #if isinstance(histos[title_formatter.format(i)], TH1F):
        if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                
            if not special_plot:
                histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
                histos.get(title_formatter.format(i), default_histo).Draw()
                histos.get(title_formatter.format(i), default_histo).SetMinimum(0)

            if y_fit_range:
                my_fit = TF1('myfit',"gaus(0)",y_fit_range[i-1][0], y_fit_range[i-1][1])
                label.DrawLatex(0.55, 0.86, '#mu = {0:6.4f}'.format(fit.GetParameter(1)))
                label.DrawLatex(0.55, 0.8, '#sigma = {0:6.4f}'.format(fit.GetParameter(2)))

            if vline:
                gPad.Update()
                ymax=gPad.GetUymax()                
                line = TLine(vline, 0,
                             vline, ymax)
                print(ymax)

                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line.Draw('same')
                root_garbage_can.append(line)
                global_dump.append(line)

            if vline2:
                gPad.Update()
                ymax=gPad.GetUymax()                
                print(ymax)
                print('drawing line HERE {}'.format(vline2))
                line2 = TLine(vline2, 0,
                              vline2, ymax)
                line2.SetLineColor(kRed)
                line2.SetLineStyle(1)
                line2.Draw('same')

            if vlines:
                gPad.Update()
                ymax=gPad.GetUymax()                
                line = TLine(vlines[i][0], 0,
                             vlines[i][0], ymax)
                line2 = TLine(vlines[i][1], 0,
                              vlines[i][1], ymax)

                print(ymax)

                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line2.SetLineColor(kRed)
                line2.SetLineStyle(1)
                line.Draw('same')
                line2.Draw('same')
                root_garbage_can.append(line)
                global_dump.append(line)            
                root_is_dumb.append(line2)
                global_dump.append(line2)
                
            
        #elif isinstance(histos[title_formatter.format(i)], TH2F):
        elif isinstance(histos.get(title_formatter.format(i), default_histo), TH2F):
            # histos[title_formatter.format(i)].Draw('colz')
            histos.get(title_formatter.format(i), default_histo).Draw('colz')
            if log:
                gPad.SetLogz() 
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass

        if vline:
            line = TLine(vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmin(),
                         vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmax())
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
            root_is_dumb.append(line)

        if vline2:
            print('drawing line')
            line2 = TLine(vline2, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmin(),
                         vline2, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmax())
            line2.SetLineColor(1)
            line2.SetLineStyle(1)
            line2.Draw('same')
            root_is_dumb.append(line2)
            root_garbage_can.append(line2)

        if hline:
            xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 
            line = TLine(xmin, hline, xmax, hline)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
            
        if hlines:
            for hl in hlines:
                xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
                xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 
                line = TLine(xmin, hl, xmax, hl)
                line.SetLineColor(1)
                line.SetLineStyle(1)
                line.Draw('same')
                root_garbage_can.append(line)
            
    
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.025, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)


    if save_fit_results:
        fit_result_out.close()
    canvas.Print(save_name)


def plot_pass_fail_page(canvas, histos, title_formatter, label, save_name,
                        xtitle=None, ytitle=None, title=None, log=False,
                        y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None,
                        hline=None, vlines=None, cut_result_formatter=None, field=None):
    
    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(1,1)
    else:
        canvas.Divide(1,1)
    canvas.cd(1)

    cut_result  = ['raw','fail','pass'] 
    color_result = [1,2,4]#kBlack, kBlue, kRed]
    cut_result_title=None
    if cut_result_formatter:
        cut_result_title = cut_result_formatter
    else:
        cut_result_title=cut_result

    result_index=0
    for rr in cut_result_title:
        result_index+=1                        
        if isinstance(histos.get(title_formatter.format(rr), default_histo), TH1F):
            
            if x_range:
                histos.get(title_formatter.format(rr), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                         
            histos.get(title_formatter.format(rr), default_histo).SetFillColorAlpha(kWhite, 0)
            histos.get(title_formatter.format(rr), default_histo).SetLineColor(color_result[result_index-1])
            histos.get(title_formatter.format(rr), default_histo).Draw('same')
            if log:
                histos.get(title_formatter.format(rr), default_histo).SetMinimum(0.1)
            label.DrawLatex(0.57, 0.8 - 0.05*(result_index-1), '#color[{}]'.format(color_result[result_index-1]) + '{' + cut_result[result_index-1] + '}')
                
            #only plot first time
            if vline and result_index==1:
                gPad.Update()
                ymax=gPad.GetUymax()                
                line = TLine(vline, 0,
                                 vline, ymax)
                print(ymax)
                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line.Draw('same')
                root_garbage_can.append(line)
                global_dump.append(line)

            if vline2 and result_index==1:
                gPad.Update()
                ymax=gPad.GetUymax()                
                print(ymax)
                print('drawing line HERE {}'.format(vline2))
                line2 = TLine(vline2, 0,
                              vline2, ymax)
                line2.SetLineColor(kRed)
                line2.SetLineStyle(1)
                line2.Draw('same')
                root_is_dumb.append(line2)
                global_dump.append(line2)
            

            if vlines and result_index==1:
                gPad.Update()
                ymax=gPad.GetUymax()                
                if log:
                    ymax = histos.get(title_formatter.format(rr), default_histo).GetMaximum()
                line = TLine(vlines[0], 0,
                             vlines[0], ymax)
                line2 = TLine(vlines[1], 0,
                              vlines[1], ymax)
                print(ymax)

                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line2.SetLineColor(kRed)
                line2.SetLineStyle(1)
                
                line.Draw('same')
                line2.Draw('same')                    
                root_garbage_can.append(line)
                global_dump.append(line)            
                root_is_dumb.append(line2)
                global_dump.append(line2)
                

    if title:
        label.DrawLatex(0.1, 0.925, title)
                
    if xtitle:
        label.DrawLatex(0.5, 0.025, xtitle)
                
    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)


    canvas.Print(save_name)

def plot_pass_fail_binned_page(canvas, histos, title_formatter, label, save_name, b_range, plot_details,
                               xtitle=None, ytitle=None, title=None, log=False,
                               y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None,
                               hline=None, vlines=None, cut_result_formatter=None, field=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(5,2)
    else:
        canvas.Divide(2,5)

        
    cc=0
    for bb in range(b_range[0], b_range[1]):
        trento_upper = plot_details['bin_center'][bb-1] + plot_details['bin_width']/2.0
        trento_lower = plot_details['bin_center'][bb-1] - plot_details['bin_width']/2.0

        canvas.cd(cc+1)
        cc+=1
        
        cut_result  = ['raw','fail','pass'] 
        color_result = [1,2,4]#kBlack, kBlue, kRed]
        cut_result_title=None
        if cut_result_formatter:
            cut_result_title = cut_result_formatter
        else:
            cut_result_title=cut_result

        result_index=0
        for rr in cut_result_title:
            result_index+=1                        
            if isinstance(histos.get(title_formatter.format(rr,bb), default_histo), TH1F):
            
                if x_range:
                    histos.get(title_formatter.format(rr,bb), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                         
                histos.get(title_formatter.format(rr,bb), default_histo).SetFillColorAlpha(kWhite, 0)
                histos.get(title_formatter.format(rr,bb), default_histo).SetLineColor(color_result[result_index-1])
                histos.get(title_formatter.format(rr,bb), default_histo).Draw('same')
                if log:
                    histos.get(title_formatter.format(rr,bb), default_histo).SetMinimum(0.1)
                label.DrawLatex(0.57, 0.8 - 0.05*(result_index-1), '#color[{}]'.format(color_result[result_index-1]) + '{' + cut_result[result_index-1] + '}')
                
                #only plot first time
                if vline and result_index==1:
                    gPad.Update()
                    ymax=gPad.GetUymax()                
                    line = TLine(vline, 0,
                                 vline, ymax)
                    print(ymax)
                    line.SetLineColor(kRed)
                    line.SetLineStyle(1)
                    line.Draw('same')
                    root_garbage_can.append(line)
                    global_dump.append(line)

                if vline2 and result_index==1:
                    gPad.Update()
                    ymax=gPad.GetUymax()                
                    print(ymax)
                    print('drawing line HERE {}'.format(vline2))
                    line2 = TLine(vline2, 0,
                                  vline2, ymax)
                    line2.SetLineColor(kRed)
                    line2.SetLineStyle(1)
                    line2.Draw('same')
                    root_is_dumb.append(line2)
                    global_dump.append(line2)
            

                if vlines and result_index==1:
                    gPad.Update()
                    ymax=gPad.GetUymax()                
                    if log:
                        ymax = histos.get(title_formatter.format(rr,bb), default_histo).GetMaximum()
                        line = TLine(vlines[0], 0,
                                     vlines[0], ymax)
                        line2 = TLine(vlines[1], 0,
                                      vlines[1], ymax)
                        print(ymax)
                        
                        line.SetLineColor(kRed)
                        line.SetLineStyle(1)
                        line2.SetLineColor(kRed)
                        line2.SetLineStyle(1)
                        
                        line.Draw('same')
                        line2.Draw('same')                    
                        root_garbage_can.append(line)
                        global_dump.append(line)            
                        root_is_dumb.append(line2)
                        global_dump.append(line2)
                
                
                if title:
                    label.DrawLatex(0.1, 0.925, title)
        
                if xtitle:
                    label.DrawLatex(0.5, 0.025, xtitle)
        
                if ytitle:
                    label.SetTextAngle(90)
                    label.DrawLatex(0.04, 0.5, ytitle)
                    label.SetTextAngle(0)
        

                if result_index == 2:
                    mean = histos.get(title_formatter.format(rr,bb), default_histo).GetMean()
                    rms = histos.get(title_formatter.format(rr,bb), default_histo).GetRMS()
                    #fit_gaus=TF1("pass_gaus","gaus(0)",mean-2*rms, mean+2*rms)
                    #histos.get(title_formatter.format(rr,bb), default_histo).Fit("pass_gaus","R")


        label.DrawLatex(0.6,0.925,'B{0}'.format(bb)+', {0:3.2f}<'.format(trento_lower)+'#phi_{Trento}'+'<{0:3.2f}'.format(trento_upper))

    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)


    canvas.Print(save_name)


def plot_page(canvas, histos, title_formatter, label, save_name,
                        xtitle=None, ytitle=None, title=None, log=False,
                        y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None, hlines = None, shades = None,
                        hline=None, vlines=None, field=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(1,1)
    else:
        canvas.Divide(1,1)
    canvas.cd(1)

    if isinstance(histos.get(title_formatter, default_histo), TH1F):
        
        if x_range:
            histos.get(title_formatter, default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
        if log:
            histos.get(title_formatter, default_histo).SetMinimum(0.1)
        histos.get(title_formatter, default_histo).Draw('same')
            
    if isinstance(histos.get(title_formatter, default_histo), TH2F):
        histos.get(title_formatter, default_histo).Draw('colz')
        
    if title:
        label.DrawLatex(0.1, 0.925, title)
        
    if xtitle:
        label.DrawLatex(0.5, 0.025, xtitle)
        
    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    if hlines:
        for hl in hlines:
            xmin = histos.get(title_formatter).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter).GetXaxis().GetXmax() 
            line = TLine(xmin, hl, xmax, hl)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
        
    if vline:
        line = TLine(vline, histos.get(title_formatter, default_histo).GetYaxis().GetXmin(),
                     vline, histos.get(title_formatter, default_histo).GetYaxis().GetXmax())
        line.SetLineColor(1)
        line.SetLineStyle(1)
        line.Draw('same')
        root_garbage_can.append(line)

    if vlines:
        for vv in vlines:
            gPad.Update()
            ymax=gPad.GetUymax()                
            line = TLine(vv, 0,
                         vv, ymax)
    
            line.SetLineColor(kRed)
            line.SetLineStyle(1)
            line.Draw('same')
            global_dump.append(line)

    if shades:
        ii=0
        while ii < len(shades):
            gPad.Update()
            ymax=gPad.GetUymax()                
            box = TBox(shades[ii], 0.0, shades[ii+1], ymax)
            box.SetFillColorAlpha(4+ii,0.5)
            ii+=2
            box.Draw("same")
            global_dump.append(box)


    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)


    canvas.Print(save_name)

def plot_slices( canvas, histos, title_formatter, save_name, label,
                      title=None, xtitle=None, ytitle=None , plot_details=None, rebinX_option=None, field=None):

    canvas.Clear()
    canvas.Divide(3,4)

    root_is_dumb = []
    page_max = 0
    signal_sum=0
    upper_bound=1.125

    h_phi = TH1F("h_phi","h_phi",10, -180, 180)
    
    for bb in range(1, plot_details['max_bins']):
        canvas.cd(bb)

        trento_upper = plot_details['bin_center'][bb-1] + plot_details['bin_width']/2.0
        trento_lower = plot_details['bin_center'][bb-1] - plot_details['bin_width']/2.0
        if isinstance(histos.get(title_formatter.format(bb), default_histo), TH1F):
            h_temp = histos.get(title_formatter.format(bb), default_histo)            
            
            if rebinX_option:
                h_temp.RebinX(rebinX_option[bb-1])

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
                mass_fit.FixParameter(1, 1.01937)
                mass_fit.FixParameter(2, 0.00458)
            else:
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
            label.DrawLatex(0.6,0.925,'B{0}, {1:3.2f}<#phi<{2:3.2f}'.format(bb, trento_lower, trento_upper ))
            label.DrawLatex(0.6,0.85,'N={0}'.format(int(n_signal)))
            label.DrawLatex(0.6,0.80,'#Delta N={0}'.format(int(d_nsignal)))

            temp_lab = TLatex()
            temp_lab.SetNDC()
            temp_lab.SetTextFont(42)
            temp_lab.SetTextSize(0.082)
            temp_lab.SetTextColor(1)
            #temp_lab.DrawLatex(0.6,0.83,'N={0}'.format(int(n_signal)))
            #temp_lab.DrawLatex(0.6,0.74,'#Delta N={0}'.format(int(d_nsignal)))

            label.DrawLatex(0.6,0.75,'#mu = {0:4.4f}'.format(fit_mean))
            label.DrawLatex(0.6,0.70,'#sigma = {0:5.5f}'.format(fit_sigma))
            label.DrawLatex(0.6,0.65,'#Delta #sigma = {0:1.5f}'.format(fit_sigma_err))
            label.DrawLatex(0.6,0.60,'#chi^{2}'+' / NDF = {0:2.2f}/{1}'.format(fit_chi2,int(fit_ndf)))

            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.65, 0.017, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.055, 0.5, ytitle)
                label.SetTextAngle(0)


            h_phi.SetBinContent(bb,n_signal)
            h_phi.SetBinError(bb,d_nsignal)

            root_is_dumb.append(label)
            root_is_dumb.append(show_back)
            root_is_dumb.append(mass_fit)

        elif isinstance(histos.get(title_formatter.format(bb), default_histo), TH2F):
            histos.get(title_formatter.format(bb), default_histo).Draw("colz")
            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.65, 0.017, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.055, 0.5, ytitle)
                label.SetTextAngle(0)
            label.DrawLatex(0.6,0.925,'B{0}, {1:3.2f}<#phi<{2:3.2f}'.format(bb, trento_lower, trento_upper ))

            root_is_dumb.append(label)

    can.Print(save_name)
    can.Clear()
    
    if isinstance(histos.get(title_formatter.format(1), default_histo), TH1F):
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

    field_setting = 'inbNoutb'
    field_setting=""
    if 'inbNout' in output_pdfname:
        field_setting="inbNoutb"
    elif 'inb' in output_pdfname:
        field_setting="inb"
    elif 'outb' in output_pdfname:
        field_setting="outb"
    else:
        field_setting=""

    # uncomment to view all histograms in file
    #for k,v in histos.items():
    #    print(k, v)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 1000, 900 )#800, 1100)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.04)
    lab.SetTextColor(1)

    chi2_lines=[0]

    can.Print('{}['.format(output_pdfname))
    
    mass_bins = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_mass_rangesVar11_'+field_setting+'.txt')
    mass_low  = mass_bins[1]
    mass_high = mass_bins[2]
    # see the makineFinalPlots.py file for event selection based on calculated missing energy and masses
    kpkm_mass_signal_bin = [mass_low, mass_high] # real values are [1.01, 1.02866]
    
    ### get the bin centers for the phi trento distribution
    #trento_centers = [histos['phitrento_pos_pass_all_lvl1'].GetBinCenter(xx) for xx in range(1, 11)]
    #print(trento_centers)
    '''
    trento_sliced = { 'max_bins':11,
                      'colors': True,                      
                      'bin_center': trento_centers,
                      'bin_width':histos['phitrento_pos_pass_all_lvl1'].GetBinWidth(1),
                      'fix_fit':True
                    }
    '''
    
    # overview of electron kinematics
    # need to name y vs x
    plot_page(can, histos, 'p_ele_theta_pass_all', lab, output_pdfname,
              xtitle='P_{el} (GeV)', ytitle='#Theta_{el} (deg)', title='Electron (FD), #Theta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_ele_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='#Theta_{el} (deg)', title='Electron (FD), #Theta vs #phi, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_ele_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='Vz_{el} (cm)', title='Electron (FD), #phi vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'theta_ele_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{el} (deg)', ytitle='Vz_{el} (cm)', title='Electron (FD), #Theta vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_ele_vz_pass_all', lab, output_pdfname,
              xtitle='P_{el} (GeV)', ytitle='Vz_{el} (cm)', title='Electron (FD), P vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_ele_dvz_elepro_pass_all', lab, output_pdfname,
              xtitle='P (GeV)', ytitle='#Delta Vz_{el} - Vz_{pro} (cm)', title='Electron (FD), P_{el} vs #Delta Vz_{el} - Vz_{pro}, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'xb_q2_pass_all', lab, output_pdfname,
              xtitle='Xb', ytitle='Q^{2} (GeV^{2})', title='(FD), Q^{2} vs Xb, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'w_q2_pass_all', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='Q^{2} (GeV^{2})', title='(FD), Q^{2} vs W, Pass All Cuts', field=field_setting)
    
   


    ## general proton kineatics
    plot_page(can, histos, 'p_pro_theta_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='#Theta_{pr} (deg)', title='Proton (FD), #Theta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_pro_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{pr} (deg)', ytitle='#Theta_{pr} (deg)', title='Proton (FD), #Theta vs #phi, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_pro_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{pr} (deg)', ytitle='Vz_{pr} (cm)', title='Proton (FD), #phi vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'theta_pro_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{pr} (deg)', ytitle='Vz_{pr} (cm)', title='Proton (FD), #Theta vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_pro_vz_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='Vz_{pr} (cm)', title='Proton (FD), P vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_pro_dvz_elepro_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='#Delta Vz (cm)', title='Proton (FD), P vs #Delta Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_pro_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_pro_chi2_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#chi^{2}', title='Proton (FD), #chi^{2} vs P, Pass All Cuts', hlines=chi2_lines, field=field_setting)

    plot_page(can, histos, 'beta_pro_chi2_pass_all', lab, output_pdfname,
                     xtitle='#beta_{pro} (GeV)', ytitle='#chi^{2}', title='Proton (FD), #chi^{2} vs #beta_{pro}, Pass All Cuts',hlines=chi2_lines, field=field_setting)
    


    ## general kaon Plus 
    plot_page(can, histos, 'p_kp_theta_pass_all', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='#Theta_{kp} (deg)', title='Proton (FD), #Theta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_kp_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{kp} (deg)', ytitle='#Theta_{kp} (deg)', title='K^{ +} (FD), #Theta vs #phi, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_kp_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{kp} (deg)', ytitle='Vz_{kp} (cm)', title='K^{ +} (FD), #phi vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'theta_kp_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{kp} (deg)', ytitle='Vz_{kp} (cm)', title='K^{ +} (FD), #Theta vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_kp_vz_pass_all', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='Vz_{kp} (cm)', title='K^{ +} (FD), P vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_kp_dvz_kpkm_pass_all', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='#Delta Vz (cm)', title='K^{ +} (FD), P vs #Delta Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_kp_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_kp_chi2_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#chi^{2}', title='K^{ +} (FD), #chi^{2} vs P, Pass All Cuts', hlines=chi2_lines, field=field_setting)

    plot_page(can, histos, 'beta_kp_chi2_pass_all', lab, output_pdfname,
                     xtitle='#beta_{K^{+}}', ytitle='#chi^{2}', title='K^{ +} (FD), #chi^{2} vs  #beta_{K^{+}}, Pass All Cuts', hlines=chi2_lines, field=field_setting)


    ## general kaon Minus
    plot_page(can, histos, 'p_km_theta_pass_all', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='#Theta_{km} (deg)', title='Proton (FD), #Theta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_km_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{km} (deg)', ytitle='#Theta_{km} (deg)', title='K^{ -} (FD), #Theta vs #phi, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'phi_km_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{km} (deg)', ytitle='Vz_{km} (cm)', title='K^{ -} (FD), #phi vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'theta_km_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{km} (deg)', ytitle='Vz_{km} (cm)', title='K^{ -} (FD), #Theta vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_km_vz_pass_all', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='Vz_{km} (cm)', title='K^{ -} (FD), P vs Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_km_dvz_kpkm_pass_all', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='#Delta Vz (cm)', title='K^{ -} (FD), P vs #Delta Vz, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_km_chi2_pass_all', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#chi^{2}', title='K^{ -} (FD), #chi^{2} vs P, Pass All Cuts', hlines=chi2_lines, field=field_setting)

    plot_page(can, histos, 'beta_km_chi2_pass_all', lab, output_pdfname,
              xtitle='#beta_{K^{-}}', ytitle='#chi^{2}', title='K^{ -} (FD), #chi^{2} vs  #beta_{K^{ -}}, Pass All Cuts', hlines=chi2_lines, field=field_setting)


    ###############
    ## general momentum resolution from reconstructed particles (not simulation resolution) of final final state.
    plot_page(can, histos, 'res_ekpkmX_mntm_vs_p_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#Delta P_{pro} (GeV)', title='(FD),#Delta P_{pro} vs P_{pro}, PassAll', field=field_setting)

    plot_page(can, histos, 'res_epkpX_mntm_vs_p_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{K^{ -}} (GeV)', ytitle='#Delta P_{K^{ -}} (GeV)', title='(FD),#Delta P_{K^{ -}} vs P_{K^{ -}}, PassAll', field=field_setting)

    plot_page(can, histos, 'res_epkmX_mntm_vs_p_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{K^{ +}} (GeV)', ytitle='#Delta P_{K^{ +}} (GeV)', title='(FD),#Delta P_{K^{ +}} vs P_{K^{ +}}, PassAll', field=field_setting)

    add_text_page(can, lab, "Missing Mass Pass All Add.", output_pdfname)
    ###############
    ## general physics kinematics
    plot_page(can,histos, 'im_kpkm_mm_ep_pass_all_pass_additional', lab, output_pdfname,
              xtitle=' MM epX (GeV)', ytitle=' I.M K^{ +}K^{ -} (GeV)', title='(FD), I.M K^{ +}K^{ -} vs MM epX, Pass All Cuts', field=field_setting)

    plot_page(can,histos, 'mm_epkp_mm_epkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle=' MM epK^{ -}X (GeV)', ytitle='MM epK^{ +}X (GeV)', title='(FD),MM epK^{ +}X vs MM epK^{ -}X, Pass All Cuts', field=field_setting)
    
    
    plot_page(can,histos, 'im_kpkm_im_pkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle='I.M. PrK^{ -} (GeV)', ytitle='I.M K^{ +}K^{ -} (GeV)', title='(FD), I.M K^{ +}K^{ -} vs I.M. PrK^{ -}, Pass All Cuts', field=field_setting)

    plot_page(can,histos, 'im_pkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle='I.M. PrK^{ -} (GeV)', ytitle='counts', title='(FD), I.M. PrK^{ -}, Pass All Cuts', shades=[1.5, 1.58, 1.78, 1.9] , field=field_setting)

    plot_page(can, histos, 'phitrento_imkpkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle='#phi_{trento} (deg)', ytitle='I.M. K^{ +}K^{ -} (GeV)', title='Invariant Mass K^{ +}K^{ -} vs #phi_{trento}', field=field_setting)

    
    
    ############################################################################################################################
    # check the calculated missing mass histograms for cuts on sigma from 1 to 4 separatly and together
    # first start with missing energy, simultaneously check the invariant mass distribution of kpkm for each bin in phi
    # 
    '''
    for lvl in range(1,3):
        add_text_page(can, lab, "Pass Exclusivity cuts at #pm {} #sigma".format(lvl), output_pdfname)
        
        plot_slices( can, histos, 'im_kpkm_trento_bin{}_pos_pass_all_lvl'+'{}'.format(lvl) , output_pdfname, lab,
                     title='FD: Inv. Mass, POS', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced)
    
        plot_slices( can, histos, 'im_kpkm_trento_bin{}_neg_pass_all_lvl'+'{}'.format(lvl) , output_pdfname, lab,
                     title='FD: Inv. Mass, NEG', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=trento_sliced)

    
    ##############
    ## check the missing energy and calculated mass2 for each phi trento bin
    ## the only items cut on are the missing energy and missing mass2 
    ## exactly 1 particle of each all in FD

    plot_pass_fail_binned_page(can, histos, 'missing_energy_trentoB{1}_{0}', lab, output_pdfname, [1,11], trento_sliced,
                               xtitle='Missing E (GeV)', ytitle='counts', title='Me, ep#rightarrow epK^{ +}K^{ -}X ', 
                               cut_result_formatter=['combssize1', 'fail_at_least_missing_energy', 'pass_all_but_missing_energy'])

    plot_pass_fail_binned_page(can, histos, 'mm_ekpkmX_trentoB{1}_{0}', lab, output_pdfname, [1,11], trento_sliced,
                               xtitle='Missing Pr (GeV)', ytitle='counts', title='MM^{2} ep#rightarrow eK^{ +}K^{ -}X ', 
                               cut_result_formatter=['combssize1', 'fail_at_least_missing_proton','pass_all_but_missing_proton'])

    plot_pass_fail_binned_page(can, histos, 'mm_epkpX_trentoB{1}_{0}', lab, output_pdfname, [1,11], trento_sliced,
                               xtitle='Missing K^{ +} (GeV)', ytitle='counts', title='MM^{2} ep#rightarrow epK^{ +}X', 
                               cut_result_formatter=['combssize1', 'fail_at_least_missing_kaonM', 'pass_all_but_missing_kaonM'])

    plot_pass_fail_binned_page(can, histos, 'mm_epkmX_trentoB{1}_{0}', lab, output_pdfname, [1,11], trento_sliced,
                               xtitle='Missing K^{ -} (GeV)', ytitle='counts', title='MM^{2} ep#rightarrow epK^{ -}X', 
                               cut_result_formatter=['combssize1', 'fail_at_least_missing_kaonP', 'pass_all_but_missing_kaonP'])

    
    
    ## check the kinematic distributions of electron, proton, kaons in these bins
    plot_slices( can, histos, 'p_ele_theta_trentoB{}_pass_all_pass_additional_pass_phi_mass' , output_pdfname, lab,
                     title='FD: #theta_{e} vs P_{e}', xtitle='P_{e} (GeV)', ytitle='#theta_{e} (deg)' , plot_details=trento_sliced)

    plot_slices( can, histos, 'p_pro_theta_trentoB{}_pass_all_pass_additional_pass_phi_mass' , output_pdfname, lab,
                     title='FD: #theta_{pr} vs P_{pr}', xtitle='P_{pr} (GeV)', ytitle='#theta_{pr} (deg)' , plot_details=trento_sliced)

    plot_slices( can, histos, 'p_kp_theta_trentoB{}_pass_all_pass_additional_pass_phi_mass' , output_pdfname, lab,
                     title='FD: #theta_{K^{ +}} vs P_{K^{ +}}', xtitle='P_{K^{ +}} (GeV)', ytitle='#theta_{K^{ +}} (deg)' , plot_details=trento_sliced)

    plot_slices( can, histos, 'p_km_theta_trentoB{}_pass_all_pass_additional_pass_phi_mass' , output_pdfname, lab,
                     title='FD: #theta_{K^{ -}} vs P_{K^{ -}}', xtitle='P_{K^{ -}} (GeV)', ytitle='#theta_{K^{ -}} (deg)' , plot_details=trento_sliced)
        
    '''
    
    ############################################################################################################################
    add_text_page(can, lab, "Pass All & Phi Cut.", output_pdfname)


    # final kinematic of particles
    plot_page(can, histos, 'p_ele_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{el} (GeV)', ytitle='#Theta_{el} (deg)', title='El,(FD), #Theta vs P,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'phi_ele_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='#Theta_{el} (deg)', title='El,(FD), #Theta vs #phi,PA,1.01<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    plot_page(can, histos, 'p_pro_theta_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='#Theta_{pr} (deg)', title='Pr,(FD), #Theta vs P,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'phi_pro_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{pr} (deg)', ytitle='#Theta_{pr} (deg)', title='Pr,(FD), #Theta vs #phi,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'p_kp_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='#Theta_{kp} (deg)', title='Pr, (FD), #Theta vs P,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'phi_kp_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{kp} (deg)', ytitle='#Theta_{kp} (deg)', title='K^{ +} (FD), #Theta vs #phi,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'p_km_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='#Theta_{km} (deg)', title='K^{ -} (FD), #Theta vs P,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'phi_km_theta_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{km} (deg)', ytitle='#Theta_{km} (deg)', title='K^{ -} (FD), #Theta vs #phi,PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)


    ## look at physics after selecting phi events
    ## this includes cuts on p, theta, and delta vz of kaons now - the additional cuts
    plot_page(can, histos, 'p_pro_beta_pass_all_cut_phi_mass', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'p_kp_beta_pass_all_cut_phi_mass', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all_cut_phi_mass', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    ## binning for 2D to inform for 1D binning of asym
    plot_page(can, histos, 'phitrento_t_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='-t (GeV^{2})', title='(FD), -t vs #phi_{trento},PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    plot_page(can, histos, 'phitrento_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='Q^{2} (GeV^{2})', title='(FD),Q^{2} vs #phi_{trento},PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'phitrento_xb_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='Xb', title='(FD),Xb vs #phi_{trento},PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    plot_page(can, histos, 'phitrento_imkpkm_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='I.M. K^{ +}K^{ -} (GeV)', title='(FD),I.M. K^{ +}K^{ -} vs #phi_{trento}, PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    
    plot_page(can, histos, 'xb_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle="X_{b}", ytitle="Q^{2} (GeV^{2})", title='Q^{2} vs X_{b},'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'w_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle="W (GeV)", ytitle="Q^{2} (GeV^{2})", title='(FD), Q^{2} vs W, PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
        
    plot_page(can, histos, 't_xb_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='X_{b}', title='X_{b} vs -t, '+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    plot_page(can, histos, 't_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='Q^{2} (GeV^{2})', title='(FD),Q^{2} vs -t, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    plot_page(can, histos, 'p_kp_t_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='-t (GeV^{2})', title='(FD), -t vs P_{K^{ +}}, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    plot_page(can, histos, 'p_pro_t_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='-t (GeV^{2})', title='(FD), -t vs P_{pro}, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    plot_page(can, histos, 'cm_cos_theta_kp_vs_t_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='cos(#Theta_{kp})', title='(FD), cos(#Theta_{kp}) vs -t, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)
    
    
    ### 1 D physics plots 
    
    plot_page(can, histos, 'q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='Q^{2} (GeV^{2})', ytitle='counts', title='(FD), Q^{2}, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 't_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='counts', title='(FD), -t, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'xb_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='X_{b}', ytitle='counts', title='(FD), X_{b}, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    plot_page(can, histos, 'phitrento_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento} (deg)', ytitle='counts', title='(FD), #phi_{trento} (deg), Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), field=field_setting)

    add_text_page(can, lab, "Binning Info For Asy.", output_pdfname)

    ### after getting the binning conditions necessary for the asymmetry measurements plot the binning boundaries on the 1D dist and the bin centers
    ### binning variables are aqurired from the monitorAsy.py code in projects/phiasy/monitorAsy.py file
        
    q2_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_q2_'+field_setting+'.txt')
    xb_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_xb_'+field_setting+'.txt')
    t_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_t_'+field_setting+'.txt')
    w_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_w_'+field_setting+'.txt')
    
    q2_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_q2_'+field_setting+'.txt')
    xb_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_xb_'+field_setting+'.txt')
    t_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_t_'+field_setting+'.txt')
    w_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_w_'+field_setting+'.txt')

    plot_page(can, histos, 'q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='Q^{2} (GeV^{2})', ytitle='counts', title='(FD), Q^{2}, Pass All, Final #phi Events',
              vlines=q2_bins, field=field_setting)

    plot_page(can, histos, 't_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='counts', title='(FD), -t, Pass All,Final #phi Events',
              vlines=t_bins, field=field_setting)

    plot_page(can, histos, 'xb_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='X_{b}', ytitle='counts', title='(FD), X_{b}, Pass All,Final #phi Events',
              vlines=xb_bins, field=field_setting)

    plot_page(can, histos, 'w_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='counts', title='(FD), W, Pass All,Final #phi Events',
              vlines=w_bins, field=field_setting)
    
    
    plot_page(can, histos, 'phitrento_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='Q^{2} (GeV^{2})', title='(FD),Q^{2} vs #phi_{trento},PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              hlines=q2_bins, field=field_setting)

    plot_page(can, histos, 'phitrento_t_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='-t (GeV^{2})', title='(FD), -t vs #phi_{trento},PA,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              hlines=t_bins, field=field_setting)

    plot_page(can, histos, 'phitrento_xb_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='Xb', title='(FD),Xb vs #phi_{trento},PA, '+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              hlines=xb_bins, field=field_setting)
    
    plot_page(can, histos, 'phitrento_imkpkm_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='#phi_{trento}', ytitle='I.M. K^{ +}K^{ -} (GeV)', title='(FD),I.M. K^{ +}K^{ -} vs #phi_{trento}, PA, '+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]), 
              hlines=[1.019], field=field_setting)

    ### 2d plots illustrating the equal entries binning procedure (add specialized binning later)
    plot_page(can, histos, 'xb_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle="X_{b}", ytitle="Q^{2} (GeV^{2})", title='Q^{2} vs X_{b},'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              vlines=xb_bins, hlines=q2_bins, field=field_setting)
        
    plot_page(can, histos, 't_xb_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='Xb', title='Xb vs -t, '+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              vlines=t_bins, hlines=xb_bins, field=field_setting)
    
    plot_page(can, histos, 't_q2_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='Q^{2} (GeV^{2})', title='(FD),Q^{2} vs -t, Pass All,'+str(kpkm_mass_signal_bin[0])+'<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              vlines=t_bins, hlines=q2_bins, field=field_setting)

    ## look at the proton km invariant mass spectrum to eventually remove each resonance to study the effect on bsa
    plot_page(can,histos, 'im_pkm_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='I.M. PrK^{ -} (GeV)', ytitle='counts', title='(FD), I.M. PrK^{ -}, Pass All,'+str(kpkm_mass_signal_bin[0]) + '<I.M K^{ +}K^{ -}<'+str(kpkm_mass_signal_bin[1]),
              shades=[1.5, 1.58, 1.78, 1.9], field=field_setting)

    
    can.Print('{}]'.format(output_pdfname))
