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
                  gPad, gStyle, TLatex, TGraphErrors,
                  TLine, kRed, kBlack, kBlue, kWhite, kTRUE, gROOT)

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
                label.DrawLatex(0.5, 0.015, xtitle)
                
            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.04, 0.5, ytitle)
                label.SetTextAngle(0)

            if field and result_index==1:
                flabel=TLatex()
                flabel.SetTextColorAlpha(1, 0.376)
                flabel.SetNDC()
                flabel.SetTextFont(42)
                flabel.SetTextSize(0.08)
                flabel.DrawLatex(0.125, 0.847, field)



    canvas.Print(save_name)

def plot_page(canvas, histos, title_formatter, label, save_name,
                        xtitle=None, ytitle=None, title=None, log=False,
                        y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None, hlines = None,
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
        if vline:
            gPad.Update()
            ymax=gPad.GetUymax()                                        
            line = TLine(vline, 0,
                         vline, ymax)
            line.SetLineColor(kRed)
            line.Draw("same")
            root_is_dumb.append(line)
        if vlines:
            for vl in vlines:
                gPad.Update()
                ymax=gPad.GetUymax()                                        
                line = TLine(vl, 0,
                             vl, ymax)

                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line.SetLineWidth(2)
                line.Draw('same')
                root_garbage_can.append(line)

            
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

    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)
        

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

    can = TCanvas('can', 'can', 1100, 1100)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)
    
    chi2_lines = [-6, -5, -4, -3, -2, -1, 1 , 2, 3, 4, 5, 6]
    field_setting = 'inb'
    

    can.Print('{}['.format(output_pdfname))
    '''

    
    # see the makineFinalPlots.py file for event selection based on calculated missing energy and masses

    # overview of final state events kinematics 
    # pass_all suffix shows the result after massing the missing masses and missing energy cuts
    # pass_all_pass_additional pass the missing masses and the coplanarity and delta vz kpkm cuts. 

    
    plot_page(can, histos, 'p_pro_theta_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#Theta_{pr} (deg)', title='Proton (FD), #Theta vs P, Pass All Cuts', field=field_setting)

    plot_page(can, histos, 'p_kp_theta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#Theta_{K^{+}} (deg)', title='K^{ +}(FD), #Theta vs P, Pass All Cuts' ,field=field_setting)

    plot_page(can, histos, 'p_km_theta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{-}} (GeV)', ytitle='#Theta_{K^{-}} (deg)', title='K^{ -}(FD), #Theta vs P, Pass All Cuts',field=field_setting)

    # can we further improve with cuts on beta vs p?
    # check beta vs p at beginning passing all exclusivity cuts
    # beta vs p
    plot_page(can, histos, 'p_pro_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, Pass All Cuts',field=field_setting)

    plot_page(can, histos, 'p_kp_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, Pass All Cuts',field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, Pass All Cuts',field=field_setting)

    # check relationship between chi2 and p
    
    plot_page(can, histos, 'p_pro_chi2_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#chi^{2}', title='Proton (FD), #chi^{2} vs P, Pass All Cuts', hlines=chi2_lines,field=field_setting)

    plot_page(can, histos, 'p_kp_chi2_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#chi^{2}', title='K^{ +} (FD), #chi^{2} vs P, Pass All Cuts', hlines=chi2_lines,field=field_setting)

    plot_page(can, histos, 'p_km_chi2_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{-}} (GeV)', ytitle='#chi^{2}', title='K^{ -}(FD), #chi^{2} vs P, Pass All Cuts', hlines=chi2_lines,field=field_setting)

    # check relationship between chi2 and beta
    plot_page(can, histos, 'beta_pro_chi2_pass_all', lab, output_pdfname,
                     xtitle='#beta_{pro} (GeV)', ytitle='#chi^{2}', title='Proton (FD), #chi^{2} vs #beta_{pro}, Pass All Cuts',hlines=chi2_lines,field=field_setting)

    plot_page(can, histos, 'beta_kp_chi2_pass_all', lab, output_pdfname,
                     xtitle='#beta_{K^{+}}', ytitle='#chi^{2}', title='K^{ +} (FD), #chi^{2} vs  #beta_{K^{+}}, Pass All Cuts', hlines=chi2_lines,field=field_setting)

    plot_page(can, histos, 'beta_km_chi2_pass_all', lab, output_pdfname,
                     xtitle='#beta_{K^{-}}', ytitle='#chi^{2}', title='K^{ -}(FD), #chi^{2} vs  #beta_{K^{-}}, Pass All Cuts', hlines=chi2_lines,field=field_setting)
    
    
    #now show invariant mass before and after cutting on chi2 from -6 to 6 in intervals of 1
    plot_pass_fail_page(can, histos, 'im_kpkm_pass_all_{}_6sigcut',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, w/ +-6#sigma #chi^{2}',field=field_setting)

    #show effect of sigma cut on beta vs p. clearly it is a direct cut on beta as function of p
    plot_page(can, histos, 'p_pro_beta_pass_all_pass_6sigcut', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, w/ +-6#sigma #chi^{2}',field=field_setting)
    
    plot_page(can, histos, 'p_kp_beta_pass_all_pass_6sigcut', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, w/ +-6#sigma #chi^{2}',field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all_pass_6sigcut', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, w/+-6#sigma #chi^{2}',field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_pass_all_{}_5sigcut',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, w/ +-5#sigma #chi^{2}',field=field_setting)

    #show effect of sigma cut on beta vs p. clearly it is a direct cut on beta as function of p
    plot_page(can, histos, 'p_pro_beta_pass_all_pass_5sigcut', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, w/ +-5#sigma #chi^{2}',field=field_setting)
    
    plot_page(can, histos, 'p_kp_beta_pass_all_pass_5sigcut', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, w/ +-5#sigma #chi^{2}',field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all_pass_5sigcut', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, w/+-5#sigma #chi^{2}',field=field_setting)



    plot_pass_fail_page(can, histos, 'im_kpkm_pass_all_{}_4sigcut',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail +-4#sigma #chi^{2}',field=field_setting)

    #show effect of sigma cut on beta vs p. clearly it is a direct cut on beta as function of p
    plot_page(can, histos, 'p_pro_beta_pass_all_pass_4sigcut', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, w/ +-4#sigma #chi^{2}',field=field_setting)
    
    plot_page(can, histos, 'p_kp_beta_pass_all_pass_4sigcut', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, w/ +-4#sigma #chi^{2}',field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all_pass_4sigcut', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, w/+-4#sigma #chi^{2}',field=field_setting)


    plot_pass_fail_page(can, histos, 'im_kpkm_pass_all_{}_3sigcut',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail +-3#sigma #chi^{2}',field=field_setting)

    #show effect of sigma cut on beta vs p. clearly it is a direct cut on beta as function of p
    plot_page(can, histos, 'p_pro_beta_pass_all_pass_3sigcut', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, w/ +-3#sigma #chi^{2}',field=field_setting)
    
    plot_page(can, histos, 'p_kp_beta_pass_all_pass_3sigcut', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, w/ +-3#sigma #chi^{2}',field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all_pass_3sigcut', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, w/+-3#sigma #chi^{2}',field=field_setting)
    '''
    '''
    plot_pass_fail_page(can, histos, 'im_kpkm_pass_all_{}_2sigcut',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail +-2#sigma #chi^{2}')

    #show effect of sigma cut on beta vs p. clearly it is a direct cut on beta as function of p
    plot_page(can, histos, 'p_pro_beta_pass_all_pass_2sigcut', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, w/ +-2#sigma #chi^{2}')
    
    plot_page(can, histos, 'p_kp_beta_pass_all_pass_2sigcut', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, w/ +-2#sigma #chi^{2}')

    plot_page(can, histos, 'p_km_beta_pass_all_pass_2sigcut', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, w/+-2#sigma #chi^{2}')


    plot_pass_fail_page(can, histos, 'im_kpkm_pass_all_{}_1sigcut',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail +-1#sigma #chi^{2}')

    #show effect of sigma cut on beta vs p. clearly it is a direct cut on beta as function of p
    plot_page(can, histos, 'p_pro_beta_pass_all_pass_1sigcut', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, w/ +-1#sigma #chi^{2}')
    
    plot_page(can, histos, 'p_kp_beta_pass_all_pass_1sigcut', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, w/ +-1#sigma #chi^{2}')

    plot_page(can, histos, 'p_km_beta_pass_all_pass_1sigcut', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, w/+-1#sigma #chi^{2}')
    '''
    
    '''


    ## cut on the vertex difference between kp and km
    plot_page(can, histos, 'p_kp_dvz_kpkm_pass_all', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#Delta Vz_{K^{+}} - Vz_{K^{-}} (cm)', title='K^{ +}(FD), #Delta Vz_{K^{+}} - Vz_{K^{-}} vs P_{K^{+}}, Pass', hlines=[-7,5],field=field_setting)

    plot_page(can, histos, 'dvz_kpkm_pass_all', lab, output_pdfname,
              xtitle='#Delta Vz_{K^{+}} - Vz_{K^{-}} (cm)', ytitle='counts', title='(FD), #Delta Vz_{K^{+}} - Vz_{K^{-}}, Pass', vlines=[-7,5],field=field_setting)


    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail #Delta Vz_{K^{ +}} - Vz_{K^{ -}}',
                        cut_result_formatter=['pass_all','pass_all_fail_dvz_kpkm','pass_all_pass_dvz_kpkm'],field=field_setting)

    ## cut on the momentum of the proton 
    ## motivated by beta vs p and high P separation of proton , kaon, pion
    plot_page(can, histos, 'p_pro_beta_pass_all', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#beta_{pro}', title='Proton (FD), #beta_{pro} vs P_{pro}, Pass All Cuts', vline=3.5,field=field_setting)
    
    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail P_{pro} Cut',
                        cut_result_formatter=['pass_all','pass_all_high_pro_p','pass_all_low_pro_p'],field=field_setting)

    ## cut on the momentum of the both charged kaons
    ## motivated by beta vs p and high P separation of proton , kaon, pion
    plot_page(can, histos, 'p_kp_beta_pass_all', lab, output_pdfname,
              xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{ +}}', title='K^{ +} (FD), #beta_{K^{ +}} vs P_{K^{ +}}, Pass All Cuts', vline=3.5,field=field_setting)

    plot_page(can, histos, 'p_km_beta_pass_all', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{ -}}', title='K^{ -} (FD), #beta_{K^{ -}} vs P_{K^{ -}}, Pass All Cuts', vline=3.5,field=field_setting)
    
    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail P_{K^{ +}},P_{K^{ -}} Cut',
                        cut_result_formatter=['pass_all','pass_all_high_kaonPM_p','pass_all_low_kaonPM_p'],field=field_setting)
    
    
    
    ## cut on the angle of the proton, and kaonP and kaonM to limit tracks strictly to FD coverage, 35 degrees    
    plot_page(can, histos, 'p_pro_theta_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#Theta_{pr} (deg)', title='Proton (FD), #Theta vs P, Pass All Cuts', hlines=[35],field=field_setting)

    plot_page(can, histos, 'p_kp_theta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#Theta_{K^{+}} (deg)', title='K^{ +}(FD), #Theta vs P, Pass All Cuts', hlines=[35],field=field_setting)

    plot_page(can, histos, 'p_km_theta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{-}} (GeV)', ytitle='#Theta_{K^{-}} (deg)', title='K^{ -}(FD), #Theta vs P, Pass All Cuts', hlines=[35],field=field_setting)
    
    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail #theta<35 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_thetaFD','pass_all_pass_thetaFD'],field=field_setting)

    plot_page(can, histos, 'cpl_pro_pass_all', lab, output_pdfname,
                     xtitle='Coplanarity #Theta_{pr}', ytitle='counts', title='Proton (FD), Coplanarity #Theta_{pr}, Pass All Cuts', vline=9,field=field_setting)

    plot_page(can, histos, 'cpl_kp_pass_all', lab, output_pdfname,
                     xtitle='Coplanarity #Theta_{K^{ +}}', ytitle='counts', title='K^{ +} (FD), Coplanarity #Theta_{K^{ +}}, Pass All Cuts', vline=9,field=field_setting)

    plot_page(can, histos, 'cpl_km_pass_all', lab, output_pdfname,
                     xtitle='Coplanarity #Theta_{K^{ -}}', ytitle='counts', title='K^{ -} (FD), Coplanarity #Theta_{K^{ -}}, Pass All Cuts', vline=9,field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail Coplan. #theta_{pr}<9 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_cpl_pro','pass_all_pass_cpl_pro'],field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail Coplan. #theta_{K^{ +}}<9 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_cpl_kp','pass_all_pass_cpl_kp'],field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail Coplan. #theta_{K^{ -}}<9 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_cpl_km','pass_all_pass_cpl_km'],field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail All Coplan. Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_cpl_all','pass_all_pass_cpl_all'],field=field_setting)
    '''
    ### check delta theta of pro, kp, km (calculated - measured)
    plot_page(can, histos, 'delta_theta_pro_pass_all', lab, output_pdfname,
                     xtitle='#Delta #Theta_{pr}', ytitle='counts', title='Proton (FD), #Delta #Theta_{pr}, Pass Me,MM2 Cuts', vlines=[-6,6],field=field_setting)

    plot_page(can, histos, 'delta_theta_kp_pass_all', lab, output_pdfname,
                     xtitle='#Delta #Theta_{K^{ +}}', ytitle='counts', title='K^{ +} (FD), #Delta #Theta_{K^{ +}}, Pass Me,MM2 Cuts', vlines=[-6,6],field=field_setting)

    plot_page(can, histos, 'delta_theta_km_pass_all', lab, output_pdfname,
                     xtitle='#Delta #Theta_{K^{ -}}', ytitle='counts', title='K^{ -} (FD), #Delta #Theta_{K^{ -}}, Pass Me,MM2 Cuts', vlines=[-6,6],field=field_setting)

    plot_page(can, histos, 'delta_theta_pro_pass_all_pass_additional', lab, output_pdfname,
                     xtitle='#Delta #Theta_{pr}', ytitle='counts', title='Proton (FD), #Delta #Theta_{pr}, Pass All Cuts', vlines=[-6,6],field=field_setting)

    plot_page(can, histos, 'delta_theta_kp_pass_all_pass_additional', lab, output_pdfname,
                     xtitle='#Delta #Theta_{K^{ +}}', ytitle='counts', title='K^{ +} (FD), #Delta #Theta_{K^{ +}}, Pass All Cuts', vlines=[-6,6],field=field_setting)

    plot_page(can, histos, 'delta_theta_km_pass_all_pass_additional', lab, output_pdfname,
                     xtitle='#Delta #Theta_{K^{ -}}', ytitle='counts', title='K^{ -} (FD), #Delta #Theta_{K^{ -}}, Pass All Cuts', vlines=[-6,6],field=field_setting)


    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail |#Delta #theta_{pr}|<6 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_delta_theta_pro','pass_all_pass_delta_theta_pro'],field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail |#Delta #theta_{K^{ +}}|<6 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_delta_theta_kp','pass_all_pass_delta_theta_kp'],field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail |#Delta #theta_{K^{ -}}|<6 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_delta_theta_km','pass_all_pass_delta_theta_km'],field=field_setting)

    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail All |#Delta #theta|<6 deg Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_delta_theta_all','pass_all_pass_delta_theta_all'],field=field_setting)
    

    '''
    ## with additional cuts
    ## chi2 max of 6
    ## delta vz of kpkm
    plot_pass_fail_page(can, histos, 'im_kpkm_{}',  lab, output_pdfname,
                        xtitle='I.M. K^{ +}K^{ -} (GeV)', ytitle='counts', title='I.M. K^{ +}K^{ -}, Pass-Fail Additional Cuts',
                        cut_result_formatter=['pass_all','pass_all_fail_additional','pass_all_pass_additional'],field=field_setting)
    
    
    '''
    can.Print('{}]'.format(output_pdfname))
