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
                  TLine, kRed, kBlack, kBlue, kWhite,  kTRUE, gROOT)

from ROOT import( RooRealVar, RooArgSet, RooArgList, RooDataHist,
                  RooGaussian, RooFit, RooAddPdf, RooPolynomial, RooFormulaVar, RooLinkedList)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
global_dump = []

gROOT.SetBatch(kTRUE)


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
                     hline=None, special_plot = None, 
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
 
            if y_fit_range:
                fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                if special_plot : 
                    print('special limits')
                    #histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', 0.79, 0.97)
                    y_fit_range[0] = 0.79
                    y_fit_range[1] = 0.97
                    #histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])
                    x = RooRealVar("x","x",0.6,1.7)
                    l = RooArgList(x)
                    data = RooDataHist("data", "data set with x1 s{}".format(i), l, histos.get(title_formatter.format(i), default_histo) )
            
                    init_mean = sector_fit_info['mean']
                    init_mean_min = sector_fit_info['S{}_mean_min'.format(i)]
                    init_mean_max = sector_fit_info['S{}_mean_max'.format(i)]
                    init_sig = sector_fit_info['sig']
                    init_sig_min = sector_fit_info['sig_min']
                    init_sig_max = sector_fit_info['sig_max']

                    mean = RooRealVar("mean","Mean of Gaussian", init_mean,init_mean_min, init_mean_max)#0.938, 0.85, 1.0 ) # fit_par[0], 0.9*fit_par[0], 1.13*fit_par[0] ) #0.938, 0.85, 1.0 )
                    sigma = RooRealVar("sigma","Width of Gaussian", init_sig, init_sig_min, init_sig_max)# 0.1, 0.009, 0.9)    #fit_par[1], 0.95*fit_par[1], 1.05*fit_par[1] )#0.1, 0.009, 0.9)
                    gauss = RooGaussian("gauss","gauss(x,mean,sigma)",x,mean,sigma)
                
                    p1 = RooRealVar("p1","coeff #1", -10, 10)
                    p2 = RooRealVar("p2","coeff #2", -10, 10) # -10, 10) -100., 100.)
                    p3 = RooRealVar("p3","coeff #3", -10, 10) #1, -1000., 1000.)
                
                    poly_bck = RooPolynomial("px", "px", x, RooArgList(p1,p2,p3))
                    #if fit_range:
                    #x.setRange("signal", 0.83, 1.05);

                    if sector_fit_info:
                        #print(' fit range specified as : {} and {} '.format( sector_fit_range[i][0], sector_fit_range[i][1] ))
                        x.setRange("signal", sector_fit_info['S{}_range'.format(i)][0], sector_fit_info['S{}_range'.format(i)][1] )#0.83, 1.05);
                    
                    Nsig = RooRealVar("Nsig","Nsig",0.0,data.numEntries())
                    Nbkg = RooRealVar("Nbkg","Nbkg",0.0,data.numEntries())
                                
                    model = RooAddPdf("model","sig+bkgd",RooArgList(gauss,poly_bck), RooArgList(Nsig, Nbkg))
                    result = model.fitTo(data, RooFit.PrintLevel(1),  RooFit.Range("signal"), RooFit.Extended(0), RooFit.Minimizer("Minuit2","migrad"))
                    #store the fit values in a TF1 but the PDF is plotted instead.               
                    fit=TF1("gaus{}".format(i),"gaus")
                    fit.SetParameter(1, mean.getVal())
                    fit.SetParameter(2, sigma.getVal())

                    mean.Print()
                    sigma.Print()   
                    
                    xframe=x.frame()
                    xframe.SetTitle(" ; ; ")

                    data.plotOn(xframe, RooFit.MarkerSize(0.35) )
                    model.plotOn(xframe, RooFit.LineColor(kRed))
                    xframe.Draw()
                    can.Draw()
                    
                    if save_fit_results:
                        fit_result_out.write(str(i) + ' ' + str(mean.getVal()) + ' ' + str(sigma.getVal()) + '\n')

                    root_is_dumb.append(result)
                    root_is_dumb.append(data)
                    root_is_dumb.append(model)
                    root_is_dumb.append(xframe)

                elif fract_fit:
                    fit = getFitFractionalHeight( histos.get(title_formatter.format(i), default_histo), histos.get(title_formatter.format(i), default_histo).GetTitle(), fract_fit, root_is_dumb)
                    if save_fit_results:
                        fit_result_out.write(str(i) + ' ' + str(fit.GetParameter(1)) + ' ' + str(fit.GetParameter(2))  + '\n')
                else:
                    histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1]) 
                    if save_fit_results:
                        fit_result_out.write(str(i) + ' ' + str(fit.GetParameter(1)) + ' ' + str(fit.GetParameter(2))  + '\n')

            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                
            if not special_plot:
                histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
                histos.get(title_formatter.format(i), default_histo).Draw()
                histos.get(title_formatter.format(i), default_histo).SetMinimum(0)

            if y_fit_range:
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


def plot_sector_pass_fail_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False,
                     y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None,
                               hline=None, vlines=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    cut_result  = ['raw','fail','pass']
    color_result = [1,2,4]#kBlack, kBlue, kRed]
    
    result_index=0
    for rr in cut_result:
        result_index+=1        
        for i in range(1,7):
            canvas.cd(i)
            if log and result_index==1:
                print(' setting log Y ')
                gPad.SetLogy() 

            if isinstance(histos.get(title_formatter.format(rr,i), default_histo), TH1F):
    
                if x_range:
                    histos.get(title_formatter.format(rr,i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                         
                histos.get(title_formatter.format(rr,i), default_histo).SetFillColorAlpha(kWhite, 0)
                histos.get(title_formatter.format(rr,i), default_histo).SetLineColor(color_result[result_index-1])
                histos.get(title_formatter.format(rr,i), default_histo).Draw('same')
                if log:
                    histos.get(title_formatter.format(rr,i), default_histo).SetMinimum(0.1)
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
                        ymax = histos.get(title_formatter.format(rr,i), default_histo).GetMaximum()
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


                if title:
                    label.DrawLatex(0.1, 0.925, title)
                        
                if xtitle:
                    label.DrawLatex(0.5, 0.015, xtitle)

                if ytitle:
                    label.SetTextAngle(90)
                    label.DrawLatex(0.04, 0.5, ytitle)
                    label.SetTextAngle(0)


    canvas.Print(save_name)


def plot_page(canvas, histos, histo_title, label, save_name,
              xtitle=None, ytitle=None, title=None, log=False,
              fits=None, x_bin_step=None, x_range=None, y_range = None, y_fit_range=None, 
              fit_method=None, merge_bins_status=None, slice_fit_info = None, manual_gfit_par=None):
    
    canvas.Clear() 
    root_is_dumb=[]
    if isinstance(histos[histo_title], TH1F):
        histos[histo_title].Draw()
    elif isinstance(histos[histo_title], TH2F):
        histos[histo_title].Draw('colz')
        if log:
            gPad.SetLogz() 

        if y_range:
               hist.GetYaxis().SetRangeUser(y_range[0], y_range[1])                        
               #hist.RebinX(x_bin_step)#int(math.floor(x_bin_step/2)))
               hist.Draw("colz")
            

        if fits:
           htitle = histo_title
           hist = histos.get(htitle, default_histo2d).Clone(htitle+'_clone')
           hist.RebinX(x_bin_step)
           graph, slices, fits, bin_info = fit_slices(hist, x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, 0)
           graph.SetMarkerStyle(8)
           graph.SetMarkerColor(2)
           graph.SetMarkerSize(1)
               
               
           graph.GetHistogram().SetMinimum(y_range[0])
           graph.GetHistogram().SetMaximum(y_range[1])
           h_x_min = hist.GetBinCenter(1)
           h_x_max = hist.GetBinCenter(hist.GetNbinsX())
           graph.GetXaxis().SetRangeUser(x_range[0], x_range[1])
           graph.SetTitle('')
           graph.Draw('P')
           root_is_dumb.append(graph)            
           root_is_dumb.append(hist)
           else:
               graph.Draw('AP')
               root_is_dumb.append(graph)
            
        if manual_gfit_par:
            print('>>> Fitting Data Points of Graph')
            
            x_min, x_max = manual_gfit_par['S{}_range'.format(0)][0], manual_gfit_par['S{}_range'.format(0)][1]
            print('>>> Using manual fit range {} {}'.format(x_min, x_max))
            poly_fit = TF1(histo_title, "[0] + [1]*x + [2]*x*x + [3]*x*x*x",x_min, x_max)            
            graph.Fit(poly_fit,'R')
            poly_fit.SetLineColor(kBlack)
            poly_fit.Draw('same')
            title_temp = manual_gfit_par['file_out_name']
            f_out = open('/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/python/fit_parameters/'+title_temp+'.txt','w')
            f_out.write('a {}'.format(poly_fit.GetParameter(0)) + ' \n')
            f_out.write('b {}'.format(poly_fit.GetParameter(1)) + ' \n')
            f_out.write('c {}'.format(poly_fit.GetParameter(2)) + ' \n')
            f_out.write('d {}'.format(poly_fit.GetParameter(3)))
            f_out.close()
            root_is_dumb.append(poly_fit)



    else:
        #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
        pass
    
    if title:
        label.DrawLatex(0.1, 0.925, title)

    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

        
    canvas.Print(save_name)
    
    if fits:
        # For the slices 
        slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
        slice_pdfname = histo_title + '_slices.pdf'
        slice_can.Print(slice_pdfname + '[')
        htitle = histo_title
        hist = histos.get(htitle, default_histo2d).Clone(htitle+'_clone2') 
        hist.RebinX(x_bin_step)
        graph, slices, fits, bin_info = fit_slices(hist, x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, 0)

        ncols = 3
        #ncols = int(np.ceil(len(slices) / nrows) + 1)
        nrows = int((slice_fit_info['S{}_max_bin_ignore'.format(0)] - slice_fit_info['S{}_min_bin_ignore'.format(0)])/ncols) + 1
        if nrows <= 2:
            nrows+=1

        print(' --------> number of cols {} '.format(ncols))
        print(' --------> number of rows {} '.format(nrows))
        

        slice_can.Clear() 
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices, fits)):
            bin_center = bin_info['bin_center_{}'.format(j)]
            if j<slice_fit_info['S0_min_bin_ignore'] or j>slice_fit_info['S0_max_bin_ignore']: continue
            print(j+1-slice_fit_info['S0_min_bin_ignore'])
            slice_can.cd(j+1-slice_fit_info['S0_min_bin_ignore'])
            #slice_can.cd(j+1)
            s.Draw()
            lab.DrawLatex(0.15, 0.85, '#mu = {0:6.4f}'.format(f.GetParameter(1)))
            lab.DrawLatex(0.15, 0.80, '#sigma = {0:6.4f}'.format(f.GetParameter(2)))
            lab.DrawLatex(0.15, 0.75, 'Bin:{}'.format(j))
            lab.DrawLatex(0.15, 0.70, 'Bin Cntr:{}'.format(bin_center))
            
        slice_can.Print(slice_pdfname)
        slice_can.Print(slice_pdfname + ']')


def plot_page_2d_overlap(canvas, histos, title_formatter, loop_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False):
    
    canvas.Clear() 
    jj=0
    for rr in loop_formatter:
        if isinstance(histos[title_formatter.format(rr)], TH2F):            
            if jj == 0:
                histos.get(title_formatter.format(rr), default_histo).Draw()
            if jj == 1:
                histos.get(title_formatter.format(rr), default_histo).Draw('colz+same')
                
            if log:
                gPad.SetLogz() 
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass
    
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)
        jj+=1
    canvas.Print(save_name)

def plot_sector_2d_page(canvas, histos, title_formatter, loop_formatter, label, save_name,
                        xtitle=None, ytitle=None, title=None,  log=False,
                        y_fit_range=None, landscape=False, x_range=None ):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    jj=0
    for rr in loop_formatter:
        for i in range(1,7):
            canvas.cd(i)
            if isinstance(histos.get(title_formatter.format(i,rr), default_histo), TH2F):
                
                if jj == 0:
                    histos.get(title_formatter.format(i,rr), default_histo).Draw()
                if jj == 1:
                    histos.get(title_formatter.format(i,rr), default_histo).Draw('colz+same')
                if log:
                    gPad.SetLogz() 
            else:
                #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
                pass


            if title:
                label.DrawLatex(0.1, 0.925, title)
                
            if xtitle:
                label.DrawLatex(0.5, 0.015, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.04, 0.5, ytitle)
                label.SetTextAngle(0)
        jj+=1

    canvas.Print(save_name)


def fit_slices(histo, x_range, x_bin_step, y_fit_range, label, xlabel, fit_method, merge_bins_status, slice_fit_info, sect=None):

    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    junk=[]
    
    x_bin_start=1
    x_bin_end=1
    ev_per_hist=0
    fraction_events_per_hist=0.2
    n_ev_per_hist=histo.GetEntries()*fraction_events_per_hist
    done = False
    i=0
    #for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
    dict_projec = {}
    bin_info = {}
    if merge_bins_status:
        print('--> performing bin merging procedure to fill slices with equal events')
        while not done:
            projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin_start, x_bin_end)
            #x_bin, x_bin + x_bin_step)
            n_ev_projec = projec.GetEntries()
            if n_ev_projec <= n_ev_per_hist:
                x_bin_end+=1
                ev_per_hist+=n_ev_projec
            else:
                print(histo.GetXaxis().GetBinCenter(x_bin_end))
                print(histo.GetXaxis().GetBinCenter(x_bin_start))
                bin_info[i] = [histo.GetXaxis().GetBinCenter(x_bin_start),  histo.GetXaxis().GetBinCenter(x_bin_end)]
                avg_bin_center = histo.GetXaxis().GetBinCenter(x_bin_start) + (histo.GetXaxis().GetBinCenter(x_bin_end) - histo.GetXaxis().GetBinCenter(x_bin_start))/2            
                projec.SetTitle('{} Bin {}; {}; counts'.format(xlabel,i,xlabel))
                temp_store = [projec,avg_bin_center]
                dict_projec[i] = temp_store
                print dict_projec
                print('projec {} bin rane {} to {} with {} events, avg bin center {}'.format(i,x_bin_start, x_bin_end, n_ev_projec, avg_bin_center))
                x_bin_start=x_bin_end
                i+=1            
                
            if x_bin_end == histo.GetNbinsX():
                done=True

    if merge_bins_status:
        for  key, value in dict_projec.items():
            fit = TF1(histo.GetTitle() + '_fit{}'.format(key), 'gaus')
            fit.SetTitle(histo.GetTitle() + '_fit{}'.format(key))
            fit.SetName(histo.GetTitle() + '_fit{}'.format(key))
        
            projec = value[0]
            avg_bin_center = value[1]
            if y_fit_range:
                fit.SetParameter(1, 0.5 * projec.GetMean() )#(y_fit_range[0] + y_fit_range[1]))
            
                #projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])
                #projec.Rebin(2)
                hist_mean = projec.GetMean()
                hist_sig = projec.GetRMS()
                y_fit_min = hist_mean - 2*hist_sig
                y_fit_max = hist_mean + 2*hist_sig
                if fit_method == 2 :
                    fit = getFitFractionalHeight( projec, projec.GetTitle(), 0.6, junk)
                elif fit_method == 1:
                    projec.Fit(fit, 'R', '',  y_fit_min,y_fit_max)

            slices.append(projec)
            fits.append(fit)
            
        x_values.append(avg_bin_center)#0.5 * (x_high + x_low))

    if not merge_bins_status:
        for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
            projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin)# + x_bin_step)
            fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus')
            fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
            fit.SetName(histo.GetTitle() + '_fit{}'.format(i))
            
            if y_fit_range:
                fit.SetParameter(1, 0.5 * (y_fit_range[0] + y_fit_range[1]))
            
                
            if 'rebinX' in slice_fit_info.keys():
                projec.Rebin(slice_fit_info['rebinX'])
            else:
                projec.Rebin(2)
            
            hist_mean = projec.GetMean()
            hist_sig = projec.GetRMS()
            y_fit_min = hist_mean - 2*hist_sig
            y_fit_max = hist_mean + 2*hist_sig
            if fit_method == 2 :
                if 'fit_fract' in slice_fit_info.keys():
                    fit = getFitFractionalHeight( projec, projec.GetTitle(), slice_fit_info['fit_fract'], junk)
                else:
                    fit = getFitFractionalHeight( projec, projec.GetTitle(), 0.6, junk)
                    #fit = getFitFractionalHeight( projec, projec.GetTitle(), 0.6, junk)
            elif fit_method == 1:
                projec.Fit(fit, 'R', '',  y_fit_min,y_fit_max)
            elif fit_method == 3 and slice_fit_info and sect:
                print('---> Using ROOFit <---')
                h_min = projec.GetXaxis().GetXmin()
                h_max = projec.GetXaxis().GetXmax()
                print(' hist x range is {} to {}'.format(h_min, h_max))

                rf_x = RooRealVar("xs","xs",h_min, h_max)
                rf_l = RooArgList(rf_x)
                rf_data = RooDataHist("slice_data","data_set_slice{}".format(i), rf_l, projec)
                init_mean = slice_fit_info['mean']
                init_mean_min = slice_fit_info['S{}_mean_min'.format(sect)]
                init_mean_max = slice_fit_info['S{}_mean_max'.format(sect)]
                init_sig = slice_fit_info['sig']
                init_sig_min = slice_fit_info['sig_min']
                init_sig_max = slice_fit_info['sig_max']
                
                slice_mean = RooRealVar("slice_mean","Mean of Gaussian", init_mean,init_mean_min, init_mean_max)#0.938, 0.85, 1.0 ) # fit_par[0], 0.9*fit_par[0], 1.13*fit_par[0] ) #0.938, 0.85, 1.0 )
                slice_sigma = RooRealVar("slice_sigma","Width of Gaussian", init_sig, init_sig_min, init_sig_max)# 0.1, 0.009, 0.9)    #fit_par[1], 0.95*fit_par[1], 1.05*fit_par[1] )#0.1, 0.009, 0.9)
                slice_gauss = RooGaussian("slice_gauss","gauss(x,slice_mean,slice_sigma)",rf_x,slice_mean,slice_sigma)
                slice_result = slice_gauss.fitTo(rf_data,RooFit.PrintLevel(1),RooFit.Extended(0), RooFit.Minimizer("Minuit2","migrad"))
                fit.SetParameter(1, slice_mean.getVal())
                fit.SetParameter(2, slice_sigma.getVal())
                slice_mean.Print()
                slice_sigma.Print()
                projec.Fit(fit, 'R', '',  slice_mean.getVal() - 3*slice_sigma.getVal(), slice_mean.getVal() + 3*slice_sigma.getVal())
                
            
                
            print('-----------> fit result chi2 ' + str(fit.GetChisquare()))
            #projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])
            
            x_low = histo.GetXaxis().GetBinCenter(x_bin)
            x_high = histo.GetXaxis().GetBinCenter(x_bin + x_bin_step)

            projec.SetTitle(xlabel +' bin ' +str(i) + ' ' + str(0.5 * (x_high + x_low)) + ' ; ' + xlabel + ' ;  counts')

            x_center = histo.GetXaxis().GetBinCenter(x_bin) #+ 0.5*(histo.GetXaxis().GetBinCenter(1) + histo.GetXaxis().GetBinCenter(2))
            print(' x bin low {} and x bin high {} -> bin center {}'.format(x_low, x_high, 0.5*(x_high + x_low)))
            print(' x center is {}'.format(x_center))
            x_values.append(x_center) #0.5 * (x_high + x_low))
            slices.append(projec)
            fits.append(fit)
            bin_info['bin_center_{}'.format(i)] = x_center #0.5*(x_high + x_low)

    means = array('d')
    means_err = array('d')
    stds = array('d')
    stds_err = array('d')
    zeros = array('d')
    x_values_temp = array('d')
    
    xx_counter  = 0
    for f in fits:
        xx_counter+=1
        if (xx_counter-1)<slice_fit_info['S{}_min_bin_ignore'.format(sect)] or (xx_counter-1)>slice_fit_info['S{}_max_bin_ignore'.format(sect)]: continue

        x_values_temp.append(x_values[xx_counter-1])
        means.append(f.GetParameter(1))
        means_err.append(f.GetParError(1))
        stds.append(f.GetParameter(2))
        stds_err.append(f.GetParError(2))
        zeros.append(0.0)
      
    #graph = TGraphErrors(len(x_values_temp), x_values_temp, means, zeros, stds)
    graph = TGraphErrors(len(x_values_temp), x_values_temp, means, zeros, means_err)
    graph.SetName('g_' + histo.GetName())

    return graph, slices, fits, bin_info
    # return np.array(x_values), np.array(means), np.array(stds), slices, fits 

 
def plot_fits(canvas, histos, x_range, x_bin_step, title_formatter,
              save_name, label, y_fit_range, fit_method, y_range=None,
              title=None, xtitle=None, ytitle=None, hline=None, vlines=None,
              merge_bins_status = False, slice_fit_info = None, manual_gfit_par=None):

    canvas.Clear()
    canvas.Divide(2,3)

    root_is_dumb = []
    for i in range(1,7):
        canvas.cd(i)

        htitle = title_formatter.format(i)
        hist = histos.get(htitle, default_histo2d).Clone(htitle+'_clone')
        hist.RebinX(x_bin_step)
        #graph, slices, fits, bin_info = fit_slices(histos.get(htitle, default_histo2d), x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, i)
        graph, slices, fits, bin_info = fit_slices(hist, x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, i)
        graph.SetMarkerStyle(8)
        graph.SetMarkerColor(2)
        graph.SetMarkerSize(1)
        #hist = histos.get(htitle, default_histo2d).Clone(htitle+'_clone')
        if y_range:
            hist.GetYaxis().SetRangeUser(y_range[0], y_range[1])                        
            #hist.RebinX(x_bin_step)#int(math.floor(x_bin_step/2)))
            hist.Draw("colz")

            graph.GetHistogram().SetMinimum(y_range[0])
            graph.GetHistogram().SetMaximum(y_range[1])
            h_x_min = hist.GetBinCenter(1)
            h_x_max = hist.GetBinCenter(hist.GetNbinsX())
            graph.GetXaxis().SetRangeUser(x_range[0], x_range[1])
            graph.SetTitle('')
            graph.Draw('P')
            root_is_dumb.append(graph)            
            root_is_dumb.append(hist)
        else:
            graph.Draw('AP')
            root_is_dumb.append(graph)
            
        if hline:
            line = TLine(x_range[0], hline, x_range[1], hline)
            line.SetLineStyle(8)
            line.SetLineWidth(1)
            line.Draw()
            root_is_dumb.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title + ' S ' + str(i))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

        if vlines:
            for sect, ll in bin_info.items():
                lmin = TLine(ll[0], y_range[0], ll[0], y_range[1])
                lmax = TLine(ll[1], y_range[0], ll[1], y_range[1])
                lmin.SetLineColor(kRed)
                lmax.SetLineColor(kRed)
                lmin.Draw('same')
                if ll[1] < graph.GetXaxis().GetXmax():
                    lmax.Draw('same')
                root_is_dumb.append(lmin)
                root_is_dumb.append(lmax)

        if manual_gfit_par:
            print('>>> Fitting Data Points of Graph')
            
            x_min, x_max = manual_gfit_par['S{}_range'.format(i)][0], manual_gfit_par['S{}_range'.format(i)][1]
            print('>>> Using manual fit range {} {}'.format(x_min, x_max))
            poly_fit = TF1(title + str(i), "[0] + [1]*x + [2]*x*x + [3]*x*x*x",x_min, x_max)            
            graph.Fit(poly_fit,'R')
            poly_fit.SetLineColor(kBlack)
            poly_fit.Draw('same')
            title_temp = manual_gfit_par['file_out_name'] + '_S' + str(i)
            f_out = open('/w/hallb-scifs17exp/clas12/bclary/CLAS12/david_elastic/elastic-clas12/python/fit_parameters/'+title_temp+'.txt','w')
            f_out.write('a {}'.format(poly_fit.GetParameter(0)) + ' \n')
            f_out.write('b {}'.format(poly_fit.GetParameter(1)) + ' \n')
            f_out.write('c {}'.format(poly_fit.GetParameter(2)) + ' \n')
            f_out.write('d {}'.format(poly_fit.GetParameter(3)))
            f_out.close()
            root_is_dumb.append(poly_fit)




    canvas.Print(save_name)

    # For the slices 
    slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
    slice_pdfname = title_formatter.split('_{}')[0] + '_slices.pdf'
    slice_can.Print(slice_pdfname + '[')
    for i in range(1,7):
        htitle = title_formatter.format(i)
        hist = histos.get(htitle, default_histo2d).Clone(htitle+'_clone2') 
        hist.RebinX(x_bin_step)
        #graph, slices, fits, bin_info = fit_slices(histos.get(htitle, default_histo2d), x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, i)
        graph, slices, fits, bin_info = fit_slices(hist, x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, i)

        # Size of slices page
        ncols = 3
        #ncols = int(np.ceil(len(slices) / nrows) + 1)
        nrows = int((slice_fit_info['S{}_max_bin_ignore'.format(i)] - slice_fit_info['S{}_min_bin_ignore'.format(i)])/ncols) + 1
        if nrows <= 2:
            nrows+=1

        print(' --------> number of cols {} '.format(ncols))
        print(' --------> number of rows {} '.format(nrows))
        

        slice_can.Clear() 
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices, fits)):
            bin_center = bin_info['bin_center_{}'.format(j)]
            if j<slice_fit_info['S{}_min_bin_ignore'.format(i)] or j>slice_fit_info['S{}_max_bin_ignore'.format(i)]: continue
            print(j+1-slice_fit_info['S{}_min_bin_ignore'.format(i)])
            slice_can.cd(j+1-slice_fit_info['S{}_min_bin_ignore'.format(i)])
            #slice_can.cd(j+1)
            s.Draw()
            lab.DrawLatex(0.15, 0.85, '#mu = {0:6.4f}'.format(f.GetParameter(1)))
            lab.DrawLatex(0.15, 0.80, '#sigma = {0:6.4f}'.format(f.GetParameter(2)))
            lab.DrawLatex(0.15, 0.75, 'Bin:{}'.format(j))
            lab.DrawLatex(0.15, 0.70, 'Bin Cntr:{}'.format(bin_center))
            
        slice_can.Print(slice_pdfname)
    slice_can.Print(slice_pdfname + ']')

                
        

def plot_overlap_1d( canvas, histos, title_formatter, hnames, save_name, label, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, logy=None, set_marker_style=None):

    canvas.Clear()
    canvas.Divide(2,3)

    root_is_dumb = []
    
    if isinstance(histos.get(title_formatter.format(hnames[0]), default_histo), TH1F) and isinstance(histos.get(title_formatter.format(hnames[1]), default_histo), TH1F):
        opt=''
        if set_marker_style:            
            opt='P+E'
            histos[title_formatter.format(hnames[0])].SetMarkerStyle(22)
            histos[title_formatter.format(hnames[1])].SetMarkerStyle(23)
            histos[title_formatter.format(hnames[0])].SetMarkerSize(1)
            histos[title_formatter.format(hnames[1])].SetMarkerSize(1)
            histos[title_formatter.format(hnames[0])].SetMarkerColor(kRed)
            histos[title_formatter.format(hnames[1])].SetMarkerColor(kBlack)
            

        if logy:
            print(' setting log Y ')
            gPad.SetLogy() 
            histos.get(title_formatter.format(hnames[0]), default_histo).SetMinimum(0.1)



        histos[title_formatter.format(hnames[0])].Draw(opt)
        histos[title_formatter.format(hnames[1])].Draw(opt+' SAME')
        
        root_is_dumb.append(histos[title_formatter.format(hnames[0])])
        root_is_dumb.append(histos[title_formatter.format(hnames[1])])

        if hline:
            line = TLine(x_range[0], hline, x_range[1], hline)
            line.SetLineStyle(8)
            line.SetLineWidth(1)
            line.Draw()
            root_is_dumb.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

        if vlines:
            for sect, ll in bin_info.items():
                lmin = TLine(ll[0], y_range[0], ll[0], y_range[1])
                lmax = TLine(ll[1], y_range[0], ll[1], y_range[1])
                lmin.SetLineColor(kRed)
                lmax.SetLineColor(kRed)
                lmin.Draw('same')
                if ll[1] < graph.GetXaxis().GetXmax():
                    lmax.Draw('same')
                root_is_dumb.append(lmin)
                root_is_dumb.append(lmax)


    else:
        print('not th1f')
    

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

    can = TCanvas('can', 'can', 800, 1100)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    ## first plot all negative and positives beta vs p
    '''
    plot_sector_page(can, histos, 'p_beta_pos_dcr2_s{}_FTOF', lab, output_pdfname,
                     xtitle='P (GeV)', ytitle='#beta', title='El_{eb}, #beta vs p, All Positives',  log=True)

    plot_sector_page(can, histos, 'p_beta_neg_dcr2_s{}_FTOF', lab, output_pdfname,
                     xtitle='P (GeV)', ytitle='#beta', title='El_{eb}, #beta vs p, All Negatives',  log=True)
    
    plot_sector_page(can, histos, 'vz_theta_pos_dcr2_s{}_FTOF', lab, output_pdfname,
                     xtitle='#theta (deg)', ytitle='vz (cm)', title='El_{eb}, Vertex_{z} vs #theta, All Positives',  log=True)

    plot_sector_page(can, histos, 'vz_theta_neg_dcr2_s{}_FTOF', lab, output_pdfname,
                     xtitle='#theta (deg)', ytitle='vz (cm)', title='El_{eb}, Vertex_{z} vs #theta, All Negatives',  log=True)

    plot_page(can, histos, 'vz_phi_pos_dcr2_FTOF', lab, output_pdfname,
                     xtitle='#phi (deg)', ytitle='vz (cm)', title='El_{eb}, Vertex_{z} vs #phi, All Positives',  log=True)

    plot_page(can, histos, 'vz_phi_neg_dcr2_FTOF', lab, output_pdfname,
                     xtitle='#phi (deg)', ytitle='vz (cm)', title='El_{eb}, Vertex_{z} vs #phi, All Negatives',  log=True)

    plot_page(can, histos, 'delta_vz_phi_pos_dcr2_FTOF', lab, output_pdfname,
                     xtitle='#phi (deg)', ytitle='#Delta vz (cm)', title='El_{eb}, #Delta Vertex_{z} vs #phi, All Positives',  log=True)

    plot_page(can, histos, 'delta_vz_phi_neg_dcr2_FTOF', lab, output_pdfname,
                     xtitle='#phi (deg)', ytitle='#Delta vz (cm)', title='El_{eb}, #Delta Vertex_{z} vs #phi, All Negatives',  log=True)

    '''

    plot_page_special(can, histos, 'p_beta_pos_FTOF', lab, output_pdfname, 
                      xtitle='Momentum (GeV)', ytitle='#beta', title='#beta vs Momentum, Positive Tracks')
    
    

    '''
    ## look at the hit position for all charged partices in dc r1, r2, r3

    plot_page(can, histos, 'dc1hit_FTOF_pass_pos', lab, output_pdfname,
              xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R1 Hit Y vs X, All Positives',  log=True)
    
    plot_page(can, histos, 'dc1hit_FTOF_pass_neg', lab, output_pdfname,
              xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R1 Hit Y vs X, All Negatives',  log=True)
    
    plot_page(can, histos, 'dc2hit_FTOF_pass_pos', lab, output_pdfname,
              xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R2 Hit Y vs X, All Positives',  log=True)
    
    plot_page(can, histos, 'dc2hit_FTOF_pass_neg', lab, output_pdfname,
              xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R2 Hit Y vs X, All Negatives',  log=True)
    
    plot_page(can, histos, 'dc3hit_FTOF_pass_pos', lab, output_pdfname,
              xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R3 Hit Y vs X, All Positives',  log=True)
    
    plot_page(can, histos, 'dc3hit_FTOF_pass_neg', lab, output_pdfname,
              xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R3 Hit Y vs X, All Negatives',  log=True)
    
    ## now look at the local hit position which we plan to cut on
    plot_sector_page(can, histos, 'local_hit_dcr1_s{}_FTOF_pass_pos', lab, output_pdfname,
                     xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R1, Local #theta vs #phi, All Positives ',  log=True)

    plot_sector_page(can, histos, 'local_hit_dcr1_s{}_FTOF_pass_neg', lab, output_pdfname,
                     xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R1, Local #theta vs #phi, All Negatives ',  log=True)
    
    #### show cut results on pip as illustrative purpose    
    plot_page_2d_overlap(can, histos, 'local_hit_dcr1_s1_FTOF_pass_pip_{}', ['fail_dcr1','pass_dcr1'], lab, output_pdfname,
                         xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R1, Local #theta vs #phi, EB #pi^{+}',  log=True)

    plot_page_2d_overlap(can, histos, 'local_hit_dcr1_s1_FTOF_pass_pim_{}', ['fail_dcr1','pass_dcr1'], lab, output_pdfname,
                         xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R1, Local #theta vs #phi, EB #pi^{-}',  log=True)
    
    #### show cut results on the DC lab frame coordinates
    plot_page_2d_overlap(can, histos, 'dc1hit_FTOF_pass_pim_{}', ['fail_dcr1','pass_dcr1'], lab, output_pdfname,
                         xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R1 Hit Y vs X, EB #pi^{-}',  log=True)
    
    
    #### now look at the remaining regions
    ############ REGION 2 ---------------------
    plot_sector_page(can, histos, 'local_hit_dcr2_s{}_FTOF_pass_pos', lab, output_pdfname,
                     xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R2, Local #theta vs #phi, All Positives ',  log=True)

    plot_sector_page(can, histos, 'local_hit_dcr2_s{}_FTOF_pass_neg', lab, output_pdfname,
                     xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R2, Local #theta vs #phi, All Negatives ',  log=True)

    #### show cut results on pip as illustrative purpose    
    plot_page_2d_overlap(can, histos, 'local_hit_dcr2_s1_FTOF_pass_pip_{}', ['fail_dcr2','pass_dcr2'], lab, output_pdfname,
                         xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R2, Local #theta vs #phi, EB #pi^{+}',  log=True)

    plot_page_2d_overlap(can, histos, 'local_hit_dcr2_s1_FTOF_pass_pim_{}', ['fail_dcr2','pass_dcr2'], lab, output_pdfname,
                         xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R2, Local #theta vs #phi, EB #pi^{-}',  log=True)
    
    #### show cut results on the DC lab frame coordinates
    plot_page_2d_overlap(can, histos, 'dc2hit_FTOF_pass_pip_{}', ['fail_dcr2','pass_dcr2'], lab, output_pdfname,
                         xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R2 Hit Y vs X, EB #pi^{-}',  log=True)
    
    ############ REGION 3 -----------------------------
    plot_sector_page(can, histos, 'local_hit_dcr3_s{}_FTOF_pass_pos', lab, output_pdfname,
                     xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R3, Local #theta vs #phi, All Positives ',  log=True)

    plot_sector_page(can, histos, 'local_hit_dcr3_s{}_FTOF_pass_neg', lab, output_pdfname,
                     xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R3, Local #theta vs #phi, All Negatives ',  log=True)

    #### show cut results on pip as illustrative purpose    
    plot_page_2d_overlap(can, histos, 'local_hit_dcr3_s1_FTOF_pass_pip_{}', ['fail_dcr3','pass_dcr3'], lab, output_pdfname,
                         xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R3, Local #theta vs #phi, EB #pi^{+}',  log=True)

    plot_page_2d_overlap(can, histos, 'local_hit_dcr3_s1_FTOF_pass_pim_{}', ['fail_dcr3','pass_dcr3'], lab, output_pdfname,
                         xtitle='Local #phi (deg)', ytitle='Local #phi (deg)', title='El_{eb}, R3, Local #theta vs #phi, EB #pi^{-}',  log=True)
    
    #### show cut results on the DC lab frame coordinates
    plot_page_2d_overlap(can, histos, 'dc3hit_FTOF_pass_pim_{}', ['fail_dcr3','pass_dcr3'], lab, output_pdfname,
                         xtitle='DC X (cm)', ytitle='DC Y (cm)', title='El_{eb}, DC R3 Hit Y vs X, EB #pi^{-}',  log=True)
    '''
    
    
    

    can.Print('{}]'.format(output_pdfname))
        
