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
## note - this is used for analyzing the result of different magnetic fields for 
## the phi analysis

from array import array 
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors, TLegend,
                  TLine, kRed, kBlack, kBlue, kWhite,  kTRUE, gROOT)

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
    
                              
    h_temp.Fit(h_name,"RN") # add Q for quiet
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
                    # plot this shit
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
            label.DrawLatex(0.5, 0.025, xtitle)

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
                    label.DrawLatex(0.5, 0.025, xtitle)

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

        if fits:
           htitle = histo_title
           hist = histos.get(htitle, default_histo2d).Clone(htitle+'_clone')
           hist.RebinX(x_bin_step)
           graph, slices, fits, bin_info = fit_slices(hist, x_range, 1, y_fit_range, label, ytitle, fit_method, merge_bins_status, slice_fit_info, 0)
           graph.SetMarkerStyle(8)
           graph.SetMarkerColor(2)
           graph.SetMarkerSize(1)
               
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
        label.DrawLatex(0.5, 0.025, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

        
    canvas.Print(save_name)
    
    if fits:
        # For the slices 
        slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
        slice_can.SetBatch(kTRUE)

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
            label.DrawLatex(0.5, 0.025, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)
        jj+=1
    canvas.Print(save_name)

def plot_page_overlap_grid(canvas, histos, histos2, histo_titles, label, save_name, nrow, ncol,
                           xtitle=None, ytitle=None, titles=None, log=False, rebin2D=None,
                           fits=None, x_bin_step=None, x_range=None, y_range = None, y_fit_range=None, 
                           fit_method=None, merge_bins_status=None, slice_fit_info = None, manual_gfit_par=None):
    
    canvas.Clear() 
    canvas.Divide(nrow,ncol)
    root_is_dumb=[]
    
    cc=0
    for histo_title in histo_titles:
        canvas.cd(cc+1)
        if isinstance(histos[histo_title], TH1F):
            opt='P+E'
            histos[histo_title].SetMarkerStyle(22)
            histos2[histo_title].SetMarkerStyle(23)
            histos[histo_title].SetMarkerSize(1)
            histos2[histo_title].SetMarkerSize(1)
            histos[histo_title].SetMarkerColor(kRed)
            histos2[histo_title].SetMarkerColor(kBlack)
            histos[histo_title].SetLineColor(kRed)
            histos2[histo_title].SetLineColor(kBlack)
            # scale to outbending data
            histos[histo_title].Scale( histos2[histo_title].GetMaximum()/histos[histo_title].GetMaximum())
            if histos2.get(histo_title, default_histo).GetMaximum() > histos.get(histo_title, default_histo).GetMaximum():
                histos2[histo_title].Draw(opt)
                histos[histo_title].Draw(opt+' SAME')
            else:
                histos[histo_title].Draw(opt)
                histos2[histo_title].Draw(opt+' SAME')

            leg = TLegend(0.7, 0.7, 0.9, 0.9)
            leg.AddEntry( histos[histo_title],'inb','pl')
            leg.AddEntry( histos2[histo_title],'out.','pl')
            leg.Draw('same')
            root_is_dumb.append(leg)


        elif isinstance(histos[histo_title], TH2F):            
            if log:
                print(' setting log Z ')
                gPad.SetLogz() 

            opt='CONT2'
            histos[histo_title].SetLineColor(kRed)
            histos2[histo_title].SetLineColor(kBlack)
            if rebin2D:
                histos[histo_title].Rebin2D(rebin2D[0], rebin2D[1])
                histos2[histo_title].Rebin2D(rebin2D[0], rebin2D[1])
            histos[histo_title].Draw(opt)
            histos2[histo_title].Draw(opt+' SAME')
        
            leg = TLegend(0.7, 0.7, 0.9, 0.9)
            leg.AddEntry( histos[histo_title],'inb','l')
            leg.AddEntry( histos2[histo_title],'out.','l')
            leg.Draw('same')
            root_is_dumb.append(leg)
            
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass
    
        if titles:
            label.DrawLatex(0.1, 0.925, titles[cc][0])

            label.DrawLatex(0.5, 0.025, titles[cc][2])

            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, titles[cc][1])
            label.SetTextAngle(0)
        cc+=1

        
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
                label.DrawLatex(0.5, 0.025, xtitle)

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
            label.DrawLatex(0.5, 0.025, xtitle)
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


def plot_2d_sidebyside(canvas,  histos_inb, histos_outb, title_formatter, save_name, label, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, logy=None, set_marker_style=None, scale=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    
    canvas.Divide(1,2)
    canvas.cd(1)
    histos_inb[title_formatter].Draw("colz")

    if hline:
        line = TLine(x_range[0], hline, x_range[1], hline)
        line.SetLineStyle(8)
        line.SetLineWidth(1)
        line.Draw()
        root_is_dumb.append(line)
        
    if title:
        label.DrawLatex(0.1, 0.925, 'INB:'+title)

    if xtitle:
        label.DrawLatex(0.5, 0.025, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.035, 0.5, ytitle)
        label.SetTextAngle(0)
    
    canvas.cd(2)
    histos_outb[title_formatter].Draw("colz")

    if hline:
        line = TLine(x_range[0], hline, x_range[1], hline)
        line.SetLineStyle(8)
        line.SetLineWidth(1)
        line.Draw()
        root_is_dumb.append(line)
        
    if title:
        label.DrawLatex(0.1, 0.925, 'OUTB:'+title)

    if xtitle:
        label.DrawLatex(0.5, 0.025, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.035, 0.5, ytitle)
        label.SetTextAngle(0)

    canvas.Print(save_name)
    

def plot_overlap_1d( canvas, histos_inb, histos_outb, title_formatter, save_name, label, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, logy=None, set_marker_style=None, scale=None, fit=None):

    canvas.Clear()
    canvas.Divide(2,3)

    root_is_dumb = []
    
    if True:# isinstance((histos_data.get(title_formatter), default_histo), TH1F) and isinstance((histos_sim.get(title_formatter), default_histo), TH1F):
        opt=''
        if set_marker_style:            
            opt='P+E'
            histos_inb[title_formatter].SetMarkerStyle(22)
            histos_outb[title_formatter].SetMarkerStyle(23)
            histos_inb[title_formatter].SetMarkerSize(1)
            histos_outb[title_formatter].SetMarkerSize(1)
            histos_inb[title_formatter].SetMarkerColor(kRed)
            histos_outb[title_formatter].SetMarkerColor(kBlack)
            histos_inb[title_formatter].SetLineColor(kRed)
            histos_outb[title_formatter].SetLineColor(kBlack)
        else:
            histos_inb[title_formatter].SetLineColor(kRed)
            histos_outb[title_formatter].SetLineColor(kBlack)

                        
            
        if logy:
            print(' setting log Y ')
            gPad.SetLogy() 
            if histos_outb.get(title_formatter, default_histo).GetMaximum() > histos_outb.get(title_formatter, default_histo).GetMaximum():
                histos_outb.get(title_formatter, default_histo).SetMinimum(0.1)
            else:
                histos_inb.get(title_formatter, default_histo).SetMinimum(0.1)


        if scale:
            if scale=="to_inb":
                histos_outb[title_formatter].Scale( histos_inb[title_formatter].GetEntries()/histos_outb[title_formatter].GetEntries())
            elif scale=="integral":                
                dxmin = histos_inb[title_formatter].GetXaxis().GetXmin()
                dxmax = histos_inb[title_formatter].GetXaxis().GetXmax()
                sxmin = histos_outb[title_formatter].GetXaxis().GetXmin()
                sxmax = histos_outb[title_formatter].GetXaxis().GetXmax()
                dbmin = histos_inb[title_formatter].GetXaxis().FindBin(dxmin)
                dbmax = histos_inb[title_formatter].GetXaxis().FindBin(dxmax)
                sbmin = histos_outb[title_formatter].GetXaxis().FindBin(sxmin)
                sbmax = histos_outb[title_formatter].GetXaxis().FindBin(sxmax)                
                histos_inb[title_formatter].Scale( histos_outb[title_formatter].GetIntegral(sbmin, sbmax)/histos_inb[title_formatter].GetIntegral(dbmin,dbmax))
            else:
                histos_inb[title_formatter].Scale( histos_outb[title_formatter].GetMaximum()/histos_inb[title_formatter].GetMaximum())


        print(' maximum of outb hist is {}'.format(histos_outb.get(title_formatter, default_histo).GetMaximum()))
        
        if histos_outb.get(title_formatter, default_histo).GetMaximum() > histos_inb.get(title_formatter, default_histo).GetMaximum():
            histos_outb[title_formatter].Draw(opt)
            histos_inb[title_formatter].Draw(opt+' SAME')
        else:
            histos_inb[title_formatter].Draw(opt)
            histos_outb[title_formatter].Draw(opt+' SAME')
            
        root_is_dumb.append(histos_inb[title_formatter])
        root_is_dumb.append(histos_outb[title_formatter])

        if hline:
            line = TLine(x_range[0], hline, x_range[1], hline)
            line.SetLineStyle(8)
            line.SetLineWidth(1)
            line.Draw()
            root_is_dumb.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.45, 0.025, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

        if vlines:
            vl=0
            for ll in vlines:
                gPad.Update()
                ymax=gPad.GetUymax()                                        
                print('{} {}'.format(vl,ll))
                line = TLine(ll, 0,
                             ll, ymax)
                line.SetLineColor(kRed)
                if vl > 1:
                    line.SetLineColor(kBlack)
                line.SetLineStyle(1)
                line.SetLineWidth(2)
                line.Draw('same')
                vl+=1
                root_is_dumb.append(line)

        leg = TLegend(0.7, 0.7, 0.9, 0.9)
        leg.AddEntry(histos_outb[title_formatter],'outb','p')
        leg.AddEntry(histos_inb[title_formatter],'inb','p')
        leg.Draw('same')
        root_is_dumb.append(leg)



        if fit:
            fit_gaus_inb = getFitFractionalHeight( histos_inb[title_formatter], histos_inb[title_formatter].GetTitle(), fit, root_is_dumb)
            fit_gaus_outb = getFitFractionalHeight( histos_outb[title_formatter], histos_outb[title_formatter].GetTitle(), fit, root_is_dumb)            
            mean_inb = histos_inb[title_formatter].GetMean()
            rms_inb = histos_inb[title_formatter].GetRMS()
            fit_gaus_poly_inb = TF1(histos_inb[title_formatter].GetTitle()+"gaus_poly","gaus(0)+pol2(3)", mean_inb-1.5*rms_inb, mean_inb+1.5*rms_inb)
            fit_gaus_poly_inb.SetParameter(0,fit_gaus_inb.GetParameter(0))
            fit_gaus_poly_inb.SetParameter(1,fit_gaus_inb.GetParameter(1))
            fit_gaus_poly_inb.SetParameter(2,fit_gaus_inb.GetParameter(2))

            mean_outb = histos_outb[title_formatter].GetMean()
            rms_outb = histos_outb[title_formatter].GetRMS()
            fit_gaus_poly_outb = TF1(histos_outb[title_formatter].GetTitle()+"outb_gaus_poly","gaus(0)+pol2(3)", mean_outb-1.5*rms_outb, mean_outb+1.5*rms_outb)
            fit_gaus_poly_outb.SetParameter(0,fit_gaus_outb.GetParameter(0))
            fit_gaus_poly_outb.SetParameter(1,fit_gaus_outb.GetParameter(1))
            fit_gaus_poly_outb.SetParameter(2,fit_gaus_outb.GetParameter(2))

            # fit with gaus + pol2 instead - result is more accurate
            fit_gaus_poly_inb.SetLineColor(kRed)
            fit_gaus_poly_outb.SetLineColor(kBlack)

            histos_inb[title_formatter].Fit(histos_inb[title_formatter].GetTitle()+"gaus_poly","R")
            histos_outb[title_formatter].Fit(histos_outb[title_formatter].GetTitle()+"outb_gaus_poly","R")

            label.DrawLatex(0.53, 0.60, 'inb: #mu = {0:.4f}'.format(fit_gaus_inb.GetParameter(1)))
            label.DrawLatex(0.53, 0.55, 'inb: #sigma = {0:.4f}'.format(fit_gaus_inb.GetParameter(2)))

            label.DrawLatex(0.53, 0.50, 'outb: #mu = {0:.4f}'.format(fit_gaus_outb.GetParameter(1)))
            label.DrawLatex(0.53, 0.45, 'out: #sigma = {0:.4f}'.format(fit_gaus_outb.GetParameter(2)))

            root_is_dumb.append(fit_gaus_outb)
            root_is_dumb.append(fit_gaus_inb)
            root_is_dumb.append(label)

    else:
        print('not th1f')
    

    canvas.Print(save_name)



def plot_overlap_2d( canvas, histos_data, histos_sim, title_formatter, save_name, label, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, logy=None, set_marker_style=None, rebin2D = None):

    canvas.Clear()
    canvas.Divide(2,3)

    root_is_dumb = []
    
    if True:# isinstance((histos_data.get(title_formatter), default_histo), TH1F) and isinstance((histos_sim.get(title_formatter), default_histo), TH1F):
        opt=''
        if set_marker_style and rebin2D:            
            opt='CONT2'
            histos_data[title_formatter].SetLineColor(kRed)
            histos_sim[title_formatter].SetLineColor(kBlack)
            histos_data[title_formatter].Rebin2D(rebin2D[0], rebin2D[1])
            histos_sim[title_formatter].Rebin2D(rebin2D[0], rebin2D[1])

        if logy:
            print(' setting log Z ')
            gPad.SetLogz() 

        histos_data[title_formatter].Draw(opt)
        histos_sim[title_formatter].Draw(opt+' SAME')
        
        leg = TLegend(0.7, 0.7, 0.9, 0.9)
        leg.AddEntry( histos_data[title_formatter],'inb','l')
        leg.AddEntry( histos_sim[title_formatter],'out.','l')
        leg.Draw('same')
        root_is_dumb.append(leg)

        root_is_dumb.append(histos_data[title_formatter])
        root_is_dumb.append(histos_sim[title_formatter])

        if hline:
            line = TLine(x_range[0], hline, x_range[1], hline)
            line.SetLineStyle(8)
            line.SetLineWidth(1)
            line.Draw()
            root_is_dumb.append(line)
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.025, xtitle)

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

def plot_overlap_1d_multiplots( canvas, histos_data, histos_sim, title_format, save_name, label, 
                                title=None, xtitle=None, ytitle=None, 
                                hline=None, vlines=None, logy=None, set_marker_style=None, scale=None, ncols=None, nbins=None):
    
    canvas.Clear()
    # Size of slices page
    #ncols = 3
    
    nrows = int(np.ceil((nbins[1]-nbins[0]) / ncols) + 1)
    if nrows <= 2:
        nrows+=1

    print(' --------> number of cols {} '.format(ncols))
    print(' --------> number of rows {} '.format(nrows))
    
    canvas.Divide(ncols, nrows)
    
    root_is_dumb = []
    for ii in range(nbins[0], nbins[1]):
        canvas.cd(ii)
        if True:# isinstance((histos_data.get(title_formatter), default_histo), TH1F) and isinstance((histos_sim.get(title_formatter), default_histo), TH1F):
            title_formatter = title_format.format(ii)
            opt=''
            if set_marker_style:            
                opt='P+E'
                histos_data[title_formatter].SetMarkerStyle(22)
                histos_sim[title_formatter].SetMarkerStyle(23)
                histos_data[title_formatter].SetMarkerSize(1)
                histos_sim[title_formatter].SetMarkerSize(1)
                histos_data[title_formatter].SetMarkerColor(kRed)
                histos_sim[title_formatter].SetMarkerColor(kBlack)
                histos_data[title_formatter].SetLineColor(kRed)
                histos_sim[title_formatter].SetLineColor(kBlack)
            else:
                histos_data[title_formatter].SetLineColor(kRed)
                histos_sim[title_formatter].SetLineColor(kBlack)
            
            
            if logy:
                print(' setting log Y ')
                gPad.SetLogy() 
                if histos_sim.get(title_formatter, default_histo).GetMaximum() > histos_sim.get(title_formatter, default_histo).GetMaximum():
                    histos_sim.get(title_formatter, default_histo).SetMinimum(0.1)
                else:
                    histos_data.get(title_formatter, default_histo).SetMinimum(0.1)


            if scale:
                if scale=="to_data":
                    histos_sim[title_formatter].Scale( histos_data[title_formatter].GetEntries()/histos_sim[title_formatter].GetEntries())
                elif scale=="integral":                
                    dxmin = histos_data[title_formatter].GetXaxis().GetXmin()
                    dxmax = histos_data[title_formatter].GetXaxis().GetXmax()
                    sxmin = histos_sim[title_formatter].GetXaxis().GetXmin()
                    sxmax = histos_sim[title_formatter].GetXaxis().GetXmax()
                    dbmin = histos_data[title_formatter].GetXaxis().FindBin(dxmin)
                    dbmax = histos_data[title_formatter].GetXaxis().FindBin(dxmax)
                    sbmin = histos_sim[title_formatter].GetXaxis().FindBin(sxmin)
                    sbmax = histos_sim[title_formatter].GetXaxis().FindBin(sxmax)                
                    histos_data[title_formatter].Scale( histos_sim[title_formatter].GetIntegral(sbmin, sbmax)/histos_data[title_formatter].GetIntegral(dbmin,dbmax))
                else:
                    histos_data[title_formatter].Scale( histos_sim[title_formatter].GetMaximum()/histos_data[title_formatter].GetMaximum())

                    
            print(' maximum of simulation hist is {}'.format(histos_sim.get(title_formatter, default_histo).GetMaximum()))
        
            if histos_sim.get(title_formatter, default_histo).GetMaximum() > histos_data.get(title_formatter, default_histo).GetMaximum():
                histos_sim[title_formatter].Draw(opt)
                histos_data[title_formatter].Draw(opt+' SAME')
            else:
                histos_data[title_formatter].Draw(opt)
                histos_sim[title_formatter].Draw(opt+' SAME')
        
            root_is_dumb.append(histos_data[title_formatter])
            root_is_dumb.append(histos_sim[title_formatter])

            if hline:
                line = TLine(x_range[0], hline, x_range[1], hline)
                line.SetLineStyle(8)
                line.SetLineWidth(1)
                line.Draw()
                root_is_dumb.append(line)
            
            if title:
                label.DrawLatex(0.1, 0.925, title.format(ii))
                
            if xtitle:
                label.DrawLatex(0.5, 0.025, xtitle)

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

            leg = TLegend(0.7, 0.7, 0.9, 0.9)
            leg.AddEntry(histos_sim[title_formatter],'gemc','p')
            leg.AddEntry(histos_data[title_formatter],'data','p')
            leg.Draw('same')
            root_is_dumb.append(leg)

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
        '-iinb',
        '--input_inb_file',
        required=True
    )
    ap.add_argument(
        '-ioutb',
        '--input_outb_file',
        required=True
    )
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    input_inb_rootfile = args.input_inb_file
    input_outb_rootfile = args.input_outb_file

    output_pdfname = args.output_prefix + '.pdf'
    rootfile_inb = TFile(input_inb_rootfile)
    rootfile_outb = TFile(input_outb_rootfile)
    histos_inb = load_histos(rootfile_inb)
    histos_outb = load_histos(rootfile_outb)

    # uncomment to view all histograms in file
    #for k,v in histos.items():
    #    print(k, v)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 1100, 1100) # was 800 by 1100
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    kpkm_mass_limit=[1.0107, 1.0287]

    
    
    #electron
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Inb vs Outb, P_{e} Final #phi (FD) Events' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Inb vs Outb, #Theta_{e} Final #phi (FD) Events' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Inb vs Outb, #Phi_{e} Final #phi (FD) Events' , xtitle='#Phi_{el} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_ele_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{el} vs P_{el}, PassAll, '+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{el} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_ele_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{el} vs #phi_{el}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{el} (deg)', ytitle='#theta (deg)')


    #proton
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, P_{p} Final #phi (FD) Events' , xtitle='P_{p} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #Theta_{p} Final #phi (FD) Events' , xtitle='#Theta_{p} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #Phi_{p} Final #phi (FD) Events' , xtitle='#Phi_{p} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{p} vs P_{p}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{p} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_pro_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{p} vs #phi_{p}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{p} (deg)', ytitle='#theta (deg)')



    #kaons
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, P_{K^{ +}} Final #phi (FD) Events' , xtitle='P_{K^{ +}} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #Theta_{K^{ +}} Final #phi (FD) Events' , xtitle='#Theta_{K^{ +}} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #Phi_{K^{ +}} Final #phi (FD) Events' , xtitle='#Phi_{K^{ +}} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_kp_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ +}} vs P_{K^{ +}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{K^{ +}} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_kp_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ +}} vs #phi_{K^{ +}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{K^{ +}} (deg)', ytitle='#theta (deg)')



    #kaonM
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, P_{K^{ -}} Final #phi (FD) Events' , xtitle='P_{K^{ -}} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #Theta_{K^{ -}} Final #phi (FD) Events' , xtitle='#Theta_{K^{ -}} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #Phi_{K^{ -}} Final #phi (FD) Events' , xtitle='#Phi_{K^{ -}} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_km_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ -}} vs P_{K^{ -}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{K^{ -}} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_km_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ -}} vs #phi_{K^{ -}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{K^{ -}} (deg)', ytitle='#theta (deg)')


    # variables that are used for events selection
    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_pro_pass_all', output_pdfname, lab, 
                     title='Inb vs Outb, #theta_{Copl} Pr. Pass Me, MM2' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True, vlines =[9,9])

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_kp_pass_all', output_pdfname, lab, 
                     title='Inb vs Outb, #theta_{Copl} K^{ +}. Pass Me, MM2' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True, vlines =[9,9])

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_km_pass_all', output_pdfname, lab, 
                     title='Inb vs Outb, #theta_{Copl} K^{ -} Pass Me, MM2' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True, vlines =[9,9])


    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #theta_{Copl} Pr. Final #phi (FD) Events' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #theta_{Copl} K^{ +}. Final #phi (FD) Events' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #theta_{Copl} K^{ -} Final #phi (FD) Events' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True)


    ##### now look at q2, xb, -t    
    plot_overlap_1d( can, histos_inb, histos_outb, 'w_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, W Final #phi (FD) Events' , xtitle='W (GeV)', ytitle='counts', set_marker_style=True, scale=True, vlines=[2])

    plot_overlap_1d( can, histos_inb, histos_outb, 'q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True, vlines=[1])

    plot_overlap_1d( can, histos_inb, histos_outb, 'xb_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, X_{b} Final #phi (FD) Events' , xtitle='X_{b}', ytitle='counts', set_marker_style=True, scale=True)
        
    plot_overlap_1d( can, histos_inb, histos_outb, 't_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, -t Final #phi (FD) Events' , xtitle='-t (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phitrento_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, #phi_{trento} Final #phi (FD) Events' , xtitle='#phi_{trento} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    ###########
    ## check with cut removing low P protons 
    #plot_overlap_1d( can, histos_data, histos_sim, 'p_ele_pass_all_pass_additional_pass_phi_mass_lowproPcut', output_pdfname, lab, 
    #                 title='FD:REC Data vs Sim, P_{e} Final #phi, P_{pr}>0.8' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, scale=True)


    #plot_overlap_1d( can, histos_data, histos_sim, 't_pass_all_pass_additional_pass_phi_mass_lowproPcut', output_pdfname, lab, 
    #                 title='FD:REC Data vs Sim, -t Final #phi, P_{pr}>0.8' , xtitle='-t (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True)


    ###########
    ## check with cuts on proton angle
    #plot_overlap_1d_multiplots(can,histos_data, histos_sim, 't_pass_all_pass_additional_pass_phi_mass_thetacut{}', output_pdfname, lab, 
    #title='FD,REC, -t,Final #phi, Pro. #Theta<35-{0}' , xtitle='-t (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True, ncols=2, nbins=[1,8])


    #plot_overlap_1d_multiplots(can,histos_data, histos_sim, 'p_ele_pass_all_pass_additional_pass_phi_mass_thetacut{}', output_pdfname, lab, 
    #                           title='FD,REC, El. P,Final #phi, Pro. #Theta<35-{0}' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, scale=True, ncols=2, nbins=[1,8])


    ###########
    # look at the relationship between -t and proton angle and electron momentum
    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_ele_t_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='-t vs P_{el}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{el} (GeV)', ytitle='-t (GeV^{2})')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'theta_pro_t_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='-t vs #Theta_{pro}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#Theta_{pro} (deg)', ytitle='-t (GeV^{2})')


    plot_2d_sidebyside(can, histos_inb, histos_outb, 'w_q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='W vs Q^{2}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='W (GeV)', ytitle='Q^{2} (GeV^{2})')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'xb_q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='Q^{2} vs X_{b} PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='X_{b}', ytitle='Q^{2} (GeV^{2})')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 't_q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='-t vs Q^{2}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]), xtitle='-t (GeV^{2})', ytitle='Q^{2} (GeV^{2})')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 't_xb_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='-t vs X_{b} PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]), xtitle='-t (GeV^{2})', ytitle='X_{b}')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 't_w_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='-t vs W PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]), xtitle='-t (GeV^{2})', ytitle='W (GeV)')
    

    ############
    ## compare missing energy and missing mass plots for cuts
    ## all in FD
    ## load cut boundaries from fits from exclusive_phi/epkpkm_top/makeFinalPlots.py output
    f_me_inb=open('../exclusive_phi/epkpkm_top/excl_phi_me_limits_pidtype1_inb.txt','r')
    f_me_outb=open('../exclusive_phi/epkpkm_top/excl_phi_me_limits_pidtype1_outb.txt','r')
    me_cuts = []
    for ff in [f_me_inb, f_me_outb]:
        for ll in ff:
            l = ll.split(' ')
            mean = float(l[0])
            sig = float(l[1])
            me_cuts.append(mean-4*sig)
            me_cuts.append(mean+4*sig)

    f_pro_inb=open('../exclusive_phi/epkpkm_top/excl_phi_mm2_pro_limits_pidtype1_inb.txt','r')
    f_pro_outb=open('../exclusive_phi/epkpkm_top/excl_phi_mm2_pro_limits_pidtype1_outb.txt','r')
    pr_cuts = []
    for ff in [f_pro_inb, f_pro_outb]:
        for ll in ff:
            l = ll.split(' ')
            print(l)
            mean = float(l[0])
            sig = float(l[1])
            pr_cuts.append(mean-4*sig)
            pr_cuts.append(mean+4*sig)

    f_kp_inb=open('../exclusive_phi/epkpkm_top/excl_phi_mm2_kp_limits_pidtype1_inb.txt','r')
    f_kp_outb=open('../exclusive_phi/epkpkm_top/excl_phi_mm2_kp_limits_pidtype1_outb.txt','r')
    kp_cuts = []
    for ff in [f_kp_inb, f_kp_outb]:
        for ll in ff:
            l = ll.split(' ')
            mean = float(l[0])
            sig = float(l[1])
            kp_cuts.append(mean-4*sig)
            kp_cuts.append(mean+4*sig)

    f_km_inb=open('../exclusive_phi/epkpkm_top/excl_phi_mm2_km_limits_pidtype1_inb.txt','r')
    f_km_outb=open('../exclusive_phi/epkpkm_top/excl_phi_mm2_km_limits_pidtype1_outb.txt','r')
    km_cuts = []
    for ff in [f_km_inb, f_km_outb]:
        for ll in ff:
            l = ll.split(' ')
            mean = float(l[0])
            sig = float(l[1])
            km_cuts.append(mean-4*sig)
            km_cuts.append(mean+4*sig)

    print(pr_cuts)

    plot_overlap_1d( can, histos_inb, histos_outb, 'missing_energy_pass_all_but_missing_energy', output_pdfname, lab, 
                     title='Inb vs Outb, Missing E,But MissE' , xtitle='Missing Energy (GeV)', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.65, vlines=me_cuts)

    plot_overlap_1d( can, histos_inb, histos_outb, 'mm_ekpkmX_pass_all_but_missing_proton', output_pdfname, lab, 
                     title='Inb vs Outb, MM^{2} Pr.,But MM^{2}_{p}' , xtitle='Missing Mass^{2} Pro. (GeV^{2})', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.65, vlines=pr_cuts)

    plot_overlap_1d( can, histos_inb, histos_outb, 'mm_epkpX_pass_all_but_missing_kaonM', output_pdfname, lab, 
                     title='Inb vs Outb, MM^{2} K^{ -}.,But MM^{2}_{K^{ -}}' , xtitle='Missing Mass^{2} K^{ -} (GeV^{2})', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.65, vlines=km_cuts)

    plot_overlap_1d( can, histos_inb, histos_outb, 'mm_epkmX_pass_all_but_missing_kaonP', output_pdfname, lab, 
                     title='Inb vs Outb, MM^{2} K^{ +}.,But MM^{2}_{K^{ +}}' , xtitle='Missing Mass^{2} K^{ +} (GeV^{2})', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.65, vlines=kp_cuts)


    ## pass all final events 
    plot_overlap_1d( can, histos_inb, histos_outb, 'missing_energy_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, Missing E,Final #phi' , xtitle='Missing Energy (GeV)', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.5)

    plot_overlap_1d( can, histos_inb, histos_outb, 'mm_ekpkmX_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, MM^{2} Pr, Final #phi' , xtitle='Missing Mass^{2} Pro. (GeV^{2})', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.5)

    plot_overlap_1d( can, histos_inb, histos_outb, 'mm_epkpX_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, MM^{2} K^{ -}.,Final #phi' , xtitle='Missing Mass^{2} K^{ -} (GeV^{2})', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.5)

    plot_overlap_1d( can, histos_inb, histos_outb, 'mm_epkmX_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='Inb vs Outb, MM^{2} K^{ +}.,Final #phi' , xtitle='Missing Mass^{2} K^{ +} (GeV^{2})', ytitle='counts', 
                     set_marker_style=True, scale=True, fit = 0.5)
    
    ## compare 2D distributions
    plot_overlap_2d(  can, histos_inb, histos_outb, 'p_ele_theta_pass_all_cut_phi_mass', output_pdfname, lab,
                      title='Ele: Inb vs Outb, #Theta vs P Final #phi (FD) Events' , xtitle='P (GeV)', ytitle='#Theta (deg)', set_marker_style=True, rebin2D=[5,5])

    plot_overlap_2d(  can, histos_inb, histos_outb, 'p_pro_theta_pass_all_cut_phi_mass', output_pdfname, lab,
                      title='Pro: Inb vs Outb, #Theta vs P Final #phi (FD) Events' , xtitle='P (GeV)', ytitle='#Theta (deg)', set_marker_style=True, rebin2D=[5,5])

    plot_overlap_2d(  can, histos_inb, histos_outb, 'p_kp_theta_pass_all_cut_phi_mass', output_pdfname, lab,
                      title='K^{ +}: Inb vs Outb, #Theta vs P Final #phi (FD) Events' , xtitle='P (GeV)', ytitle='#Theta (deg)', set_marker_style=True, rebin2D=[5,5])

    plot_overlap_2d(  can, histos_inb, histos_outb, 'p_km_theta_pass_all_cut_phi_mass', output_pdfname, lab,
                      title='K^{ -}: Inb vs Outb, #Theta vs P Final #phi (FD) Events' , xtitle='P (GeV)', ytitle='#Theta (deg)', set_marker_style=True, rebin2D=[5,5])
    
    
    ##############################################################################3
    ###########################################
    # grid the results for easy display 
    histo_titles = [ 'p_ele_pass_all_pass_additional_pass_phi_mass',
                     'theta_ele_pass_all_pass_additional_pass_phi_mass',
                     'phi_ele_pass_all_pass_additional_pass_phi_mass',
                     'p_pro_pass_all_pass_additional_pass_phi_mass',
                     'theta_pro_pass_all_pass_additional_pass_phi_mass',
                     'phi_pro_pass_all_pass_additional_pass_phi_mass',
                     'p_kp_pass_all_pass_additional_pass_phi_mass',
                     'theta_kp_pass_all_pass_additional_pass_phi_mass',
                     'phi_kp_pass_all_pass_additional_pass_phi_mass',
                     'p_km_pass_all_pass_additional_pass_phi_mass',
                     'theta_km_pass_all_pass_additional_pass_phi_mass',
                     'phi_km_pass_all_pass_additional_pass_phi_mass']

    titles_gen = [ ['Inb vs Outb: Electron P ', 'counts','P (GeV)'],
                   ['Inb vs Outb: Electron #theta', 'counts','#theta (deg)'],
                   ['Inb vs Outb: Electron #phi ', 'counts','#phi (deg)'],
                   ['Inb vs Outb: Proton P', 'counts','P (GeV)'],
                   ['Inb vs Outb: Proton #theta', 'counts','#theta (deg)'],
                   ['Inb vs Outb: Proton #phi', 'counts','#phi (deg)'],
                   ['Inb vs Outb: K^{ +} P', 'counts','P (GeV)'],
                   ['Inb vs Outb: K^{ +} #theta', 'counts','#theta (deg)'],
                   ['Inb vs Outb: K^{ +} #phi', 'counts','#phi (deg)'],
                   ['Inb vs Outb: K^{ -} P ', 'counts','P (GeV)'],
                   ['Inb vs Outb: K^{ -} #theta', 'counts','#theta (deg)'],
                   ['Inb vs Outb: K^{ -} #phi ', 'counts','#phi (deg)']]
               


    plot_page_overlap_grid(can, histos_inb, histos_outb, histo_titles, lab, output_pdfname, 3, 4,
                   titles=titles_gen, log=False)

    ###### 
    ## 2d grid overlap
    
    histo_titles = [ 'p_ele_theta_pass_all_pass_additional_pass_phi_mass',
                     'phi_ele_theta_pass_all_pass_additional_pass_phi_mass',
                     'p_pro_theta_pass_all_pass_additional_pass_phi_mass',
                     'phi_pro_theta_pass_all_pass_additional_pass_phi_mass',
                     'p_kp_theta_pass_all_pass_additional_pass_phi_mass',
                     'phi_kp_theta_pass_all_pass_additional_pass_phi_mass',
                     'p_km_theta_pass_all_pass_additional_pass_phi_mass',
                     'phi_km_theta_pass_all_pass_additional_pass_phi_mass']
    titles_gen = [ ['Inb vs Outb: Electron P vs #theta', '#theta (deg)','P (GeV)'],
                   ['Inb vs Outb: Electron #phi vs #theta', '#theta (deg)','#phi (deg)'],
                   ['Inb vs Outb: Proton P vs #theta', '#theta (deg)','P (GeV)'],
                   ['Inb vs Outb: Proton #phi vs #theta', '#theta (deg)','#phi (deg)'],
                   ['Inb vs Outb: K^{ +} vs #theta', '#theta (deg)','P (GeV)'],
                   ['Inb vs Outb: K^{ +} #phi vs #theta', '#theta (deg)','#phi (deg)'],
                   ['Inb vs Outb: K^{ -} vs #theta', '#theta (deg)','P (GeV)'],
                   ['Inb vs Outb: K^{ -} #phi vs #theta', '#theta (deg)','#phi (deg)']]
               

    plot_page_overlap_grid(can, histos_inb, histos_outb, histo_titles, lab, output_pdfname, 2, 4,
                   titles=titles_gen, log=False)
    

    histo_titles_kins = ['xb_q2_pass_all_pass_additional_pass_phi_mass',
                         't_q2_pass_all_pass_additional_pass_phi_mass',
                         'w_q2_pass_all_pass_additional_pass_phi_mass',
                         't_xb_pass_all_pass_additional_pass_phi_mass',
                         't_w_pass_all_pass_additional_pass_phi_mass']
    # title , y title , x title
    titles_kin_gen = [ ['Inb vs Outb: Q^{2} vs x_{B}','Q^{2} (GeV^{2})','x_{B}'],
                       ['Inb vs Outb: Q^{2} vs -t',' Q^{2} (GeV^{2})','-t (GeV^{2})'],
                       ['Inb vs Outb: Q^{2} vs W', ' Q^{2} (GeV^{2})', 'W (GeV)'],
                       ['Inb vs Outb: -t vs x_{B}','x_{B}','-t (GeV^{2})'],
                       ['Inb vs Outb: W vs -t', 'W (GeV)','-t (GeV^{2})'] ]

    plot_page_overlap_grid(can, histos_inb, histos_outb, histo_titles_kins, lab, output_pdfname, 2, 3,
                   titles=titles_kin_gen, log=False)

    histo_titles_1dkin = [ 'w_pass_all_pass_additional_pass_phi_mass',
                           'q2_pass_all_pass_additional_pass_phi_mass',
                           'xb_pass_all_pass_additional_pass_phi_mass',
                           't_pass_all_pass_additional_pass_phi_mass',
                           'phitrento_pass_all_pass_additional_pass_phi_mass']
    titles_kin1d = [ ['Inb vs Outb: W','counts','W (GeV)'],
                     ['Inb vs Outb: Q^{2}','counts','Q^{2} (GeV)'],
                     ['Inb vs Outb: x_{B}','counts','x_{B}'],
                     ['Inb vs Outb: -t','counts','-t (GeV^{2})'],
                     ['Inb vs Outb: #phi_{trento}','counts','#phi_{trento} (deg)']
                   ]
    plot_page_overlap_grid(can, histos_inb, histos_outb, histo_titles_1dkin, lab, output_pdfname, 2, 3,
                           titles=titles_kin1d, log=False)


    can.Print('{}]'.format(output_pdfname))
        
