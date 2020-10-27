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
                  gPad, gStyle, TLatex, TGraphErrors, TLegend,
                  TLine, kRed, kBlack, kBlue, kWhite, kGreen)

from ROOT import( RooRealVar, RooArgSet, RooArgList, RooDataHist,
                  RooGaussian, RooFit, RooAddPdf, RooPolynomial, RooFormulaVar, RooLinkedList)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
global_dump = []

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
            label.DrawLatex(0.5, 0.015, xtitle)

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
                        hline=None, vlines=None, cut_result_formatter=None):
    
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
                label.DrawLatex(0.5, 0.015, xtitle)
                
            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.04, 0.5, ytitle)
                label.SetTextAngle(0)


    canvas.Print(save_name)

def plot_page(canvas, histos, title_formatter, label, save_name,
                        xtitle=None, ytitle=None, title=None, log=False,
                        y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None, hlines = None,
                        hline=None, vlines=None):

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
        label.DrawLatex(0.5, 0.015, xtitle)
        
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
        '-i1',
        '--input_file1',
        required=True
    )
    ap.add_argument(
        '-i2',
        '--input_file2',
        required=True
    )
    ap.add_argument(
        '-i3',
        '--input_file3',
        required=True
    )
    ap.add_argument(
        '-i4',
        '--input_file4',
        required=True
    )
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    input_rootfile1 = args.input_file1 
    input_rootfile2 = args.input_file2 
    input_rootfile3 = args.input_file3 
    input_rootfile4 = args.input_file4 
    output_pdfname = args.output_prefix + '.pdf'
    rootfile1 = TFile(input_rootfile1)
    rootfile2 = TFile(input_rootfile2)
    rootfile3 = TFile(input_rootfile3)
    rootfile4 = TFile(input_rootfile4)
    histos1 = load_histos(rootfile1)
    histos2 = load_histos(rootfile2)
    histos3 = load_histos(rootfile3)
    histos4 = load_histos(rootfile4)

    # uncomment to view all histograms in file
    #for k,v in histos.items():
    #    print(k, v)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 800, 1100)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)
    
    can.Print('{}['.format(output_pdfname))
    #gPad.SetLogy()

    phi_pass_all = 'im_kpkm_pass_all'

    imkpkm1 = histos1.get(phi_pass_all, default_histo)
    imkpkm2 = histos2.get(phi_pass_all, default_histo)
    imkpkm3 = histos3.get(phi_pass_all, default_histo)
    imkpkm4 = histos4.get(phi_pass_all, default_histo)

    opt='P+E'
    imkpkm1.SetMarkerColor(kRed)
    imkpkm2.SetMarkerColor(kBlack)
    imkpkm3.SetMarkerColor(kGreen)
    imkpkm4.SetMarkerColor(kBlue)

    imkpkm1.SetMarkerStyle(22)
    imkpkm2.SetMarkerStyle(23)
    imkpkm3.SetMarkerStyle(24)
    imkpkm4.SetMarkerStyle(25)

    #imkpkm1.SetMinimum(0.1)
    imkpkm1.Draw(opt)
    imkpkm2.Draw('same+'+opt)
    imkpkm3.Draw('same+'+opt)
    imkpkm4.Draw('same+'+opt)

    leg=TLegend(0.5,0.5, 0.9, 0.9)
    leg.AddEntry(imkpkm1,'El EB, Had. EB','p')
    leg.AddEntry(imkpkm2,'El Mod, Had. EB','p')
    leg.AddEntry(imkpkm3,'El EB, Had. Mod','p')
    leg.AddEntry(imkpkm4,'El Mod, Had. Mod','p')
    leg.Draw('same')

    title='Comparison of PID Methods'
    xtitle='I.M. K^{ +}K^{ -}'
    ytitle='counts'

    if title:
        lab.DrawLatex(0.1, 0.925, title)
        
    if xtitle:
        lab.DrawLatex(0.5, 0.015, xtitle)
        
    if ytitle:
        lab.SetTextAngle(90)
        lab.DrawLatex(0.04, 0.5, ytitle)
        lab.SetTextAngle(0)

    can.Print(output_pdfname)
    can.Print('{}]'.format(output_pdfname))
