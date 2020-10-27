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
        label.DrawLatex(0.5, 0.015, xtitle)

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

def plot_overlap_1d( canvas, histos_data, histos_sim, title_formatter, save_name, label, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, logy=None, set_marker_style=None, scale=None, fit=None):

    canvas.Clear()
    canvas.Divide(1,1)

    root_is_dumb = []
    
    if True:# isinstance((histos_data.get(title_formatter), default_histo), TH1F) and isinstance((histos_sim.get(title_formatter), default_histo), TH1F):
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

        leg = TLegend(0.7, 0.7, 0.9, 0.9)
        leg.AddEntry(histos_sim[title_formatter],'data2','p')
        leg.AddEntry(histos_data[title_formatter],'data1','p')
        leg.Draw('same')
        root_is_dumb.append(leg)



        if fit:
            #mean_data = histos_data[title_formatter].GetMean()
            #sig_data = histos_data[title_formatter].GetRMS()
            fit_gaus_data = getFitFractionalHeight( histos_data[title_formatter], histos_data[title_formatter].GetTitle(), fit, root_is_dumb)
            fit_gaus_sim = getFitFractionalHeight( histos_sim[title_formatter], histos_sim[title_formatter].GetTitle(), fit, root_is_dumb)            
            mean_data = histos_data[title_formatter].GetMean()
            rms_data = histos_data[title_formatter].GetRMS()
            fit_gaus_poly_data = TF1(histos_data[title_formatter].GetTitle()+"gaus_poly","gaus(0)+pol2(3)", mean_data-1.5*rms_data, mean_data+1.5*rms_data)
            fit_gaus_poly_data.SetParameter(0,fit_gaus_data.GetParameter(0))
            fit_gaus_poly_data.SetParameter(1,fit_gaus_data.GetParameter(1))
            fit_gaus_poly_data.SetParameter(2,fit_gaus_data.GetParameter(2))

            mean_sim = histos_sim[title_formatter].GetMean()
            rms_sim = histos_sim[title_formatter].GetRMS()
            fit_gaus_poly_sim = TF1(histos_sim[title_formatter].GetTitle()+"sim_gaus_poly","gaus(0)+pol2(3)", mean_sim-1.5*rms_sim, mean_sim+1.5*rms_sim)
            fit_gaus_poly_sim.SetParameter(0,fit_gaus_sim.GetParameter(0))
            fit_gaus_poly_sim.SetParameter(1,fit_gaus_sim.GetParameter(1))
            fit_gaus_poly_sim.SetParameter(2,fit_gaus_sim.GetParameter(2))

            # fit with gaus + pol2 instead - result is more accurate
            fit_gaus_poly_data.SetLineColor(kRed)
            fit_gaus_poly_sim.SetLineColor(kBlack)

            histos_data[title_formatter].Fit(histos_data[title_formatter].GetTitle()+"gaus_poly","R")
            histos_sim[title_formatter].Fit(histos_sim[title_formatter].GetTitle()+"sim_gaus_poly","R")
            

            label.DrawLatex(0.56, 0.60, '#color[2]{#mu_{data}}'+' = {0:.4f}'.format(fit_gaus_poly_data.GetParameter(1)))
            label.DrawLatex(0.56, 0.55, '#color[2]{#sigma_{data}}'+' = {0:.4f}'.format(fit_gaus_poly_data.GetParameter(2)))

            label.DrawLatex(0.56, 0.50, '#color[1]{#mu_{gemc}}'+' = {0:.4f}'.format(fit_gaus_poly_sim.GetParameter(1)))
            label.DrawLatex(0.56, 0.45, '#color[1]{#sigma_{gemc}}'+' = {0:.4f}'.format(fit_gaus_poly_sim.GetParameter(2)))
            #fit_gaus_data.Draw("same")
            #fit_gaus_sim.Draw("same")
            fit_gaus_poly_data.Draw("same")
            fit_gaus_poly_sim.Draw("same")
            root_is_dumb.append(fit_gaus_poly_sim)
            root_is_dumb.append(fit_gaus_poly_data)
            root_is_dumb.append(label)

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
        '-id1',
        '--input_data1_file',
        required=True
    )
    ap.add_argument(
        '-id2',
        '--input_data2_file',
        required=True
    )
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    input_data1_rootfile = args.input_data1_file
    input_data2_rootfile = args.input_data2_file

    output_pdfname = args.output_prefix + '.pdf'
    rootfile_data1 = TFile(input_data1_rootfile)
    rootfile_data2 = TFile(input_data2_rootfile)
    histos_data1 = load_histos(rootfile_data1)
    histos_data2 = load_histos(rootfile_data2)

    # uncomment to view all histograms in file
    #for k,v in histos.items():
    #    print(k, v)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 1100, 1100)#800, 1100)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    kpkm_mass_limit=[1.0107, 1.0287]

    
    ## first compare the final state selected phi reconstructed phi events in GEMC to generated distribution
    # electrons
        
    #electron
    # Validate output of codes between two different run through of the code on data
    plot_overlap_1d( can, histos_data1, histos_data2, 'p_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Data1 vs Data2, P_{e} Final #phi (FD) Events' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can,  histos_data1, histos_data2, 'theta_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Data1 vs Data2, #Theta_{e} Final #phi (FD) Events' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can,  histos_data1, histos_data2, 'phi_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Data1 vs Data2, #Phi_{e} Final #phi (FD) Events' , xtitle='#Phi_{el} (deg)', ytitle='counts', set_marker_style=True, scale=True)


    ##
    plot_overlap_1d( can,  histos_data1, histos_data2, 'q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Data1 vs Data2, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can,  histos_data1, histos_data2, 'xb_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Data1 vs Data2, xb Final #phi (FD) Events' , xtitle='xb', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can,  histos_data1, histos_data2, 't_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC Data1 vs Data2, -t Final #phi (FD) Events' , xtitle='-t (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True)

    #plot_overlap_1d( can,  histos_data1, histos_data2, 'w_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
     #                title='REC Data1 vs Data2, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV)', ytitle='counts', set_marker_style=True, scale=True)


    




    can.Print('{}]'.format(output_pdfname))
