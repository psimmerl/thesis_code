#!/usr/bin/env python 
import argparse 
import numpy as np
import math
import os
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
                  gPad, gStyle, TLatex, TGraphErrors, TMultiGraph,
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
        if y_fit_range:
            projec.Fit(fit, 'R', '',  y_fit_range[0],y_fit_range[1])
        else:
            projec.Fit(fit, 'R', '',  y_fit_min,y_fit_max)

        slices.append(projec)
        fits.append(fit)
            
        


def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None,  log=False, particle_name=None, sig_frac=None, 
                     y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None,
                     hline=None, special_plot = None, y_range= None, calc_sigma = None,
                     fit_par = None, fit_range=None, sector_fit_info=None, fract_fit=None, save_fit_results=None, vlines=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 

        
    x_start=x_range[0]
    x_stop=x_range[1]

    slices_sect = []
    l_mg = []
    for ss in range(1,7):
    
        root_is_dumb.append(histos.get(title_formatter.format(ss), default_histo))        
        htitle = title_formatter.format(ss)

        histo = histos.get(htitle, default_histo2d).Clone(htitle+'_clone')
        histo.RebinX(22)
        
        means=array('d')
        sigma=array('d')
        psigma=array('d')
        x_centers=array('d')
        x_err =array('d')
        root_is_dumb.append(histo)

        #####################
        ## parameter for pol3 for calculating sigma to initilaize the fit functions
        f_pos_in=None
        f_neg_in=None
	f_mean = None
        if calc_sigma and os.path.exists('{0}_pos_out_poly3_par_s{1}.txt'.format(particle_name, ss)):
            f_pos_in = open('{0}_pos_out_poly3_par_s{1}.txt'.format(particle_name, ss), 'r')
            
        if calc_sigma and os.path.exists('{0}_neg_out_poly3_par_s{1}.txt'.format(particle_name, ss)):
            f_neg_in = open('{0}_neg_out_poly3_par_s{1}.txt'.format(particle_name, ss), 'r')
       
	pos_sigma_func = TF1("pos_sigma_fit_pol3_s".format(ss),"pol3") 
       
        if calc_sigma and os.path.exists('{0}_mean_out_poly3_par_s{1}.txt'.format(particle_name, ss)):
            f_mean_in = open('{0}_mean_out_poly3_par_s{1}.txt'.format(particle_name, ss), 'r')
        mean_func = TF1("mean_fit_pol3_s".format(ss),"pol1") 
 	
	# set parameter values
	if calc_sigma and os.path.exists('{0}_pos_out_poly3_par_s{1}.txt'.format(particle_name, ss)):
	    cc=0
            for ll in f_pos_in:
            	print("ll in {}".format(ll))
            	pos_sigma_func.SetParameter(cc, float(ll))
            	cc+=1

	if calc_sigma and os.path.exists('{0}_mean_out_poly3_par_s{1}.txt'.format(particle_name, ss)):
            cc=0
            for ll in f_mean_in:
        	mean_func.SetParameter(cc, float(ll))
                cc+=1
        ##### finish reading in parameters for setting sigma
        
        slices = []
        x_bin_step=1
        for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
            projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin)# + x_bin_step)
            projec.Rebin(2)
            # center of momentum bin
            x_center = histo.GetXaxis().GetBinCenter(x_bin)
	    
		
            y_min = -0.55
            y_max = 0.55
            if y_range:
                y_min = y_range[0]
                y_max = y_range[1]
            if i < 2:
                y_min = -0.5
                y_max = 0.5
                
	    pos_cal_sigma = 0
	    calc_mean = 0
            if calc_sigma:
                #f_cal_pos_sig = TF1("cal_pos_sig_s{}".format(ss))
                #f_cal_neg_sig = TF1("cal_neg_sig_s{}".format(ss))
                pos_cal_sigma = pos_sigma_func.Eval(x_center)
		calc_mean = mean_func.Eval(x_center)	
	    else:
		pos_cal_sigma=0.15

 	    # add code to change n sigma range to be smaller at larger P	   
	    frac_sigma = 1.0
	     
            fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus(0)+pol3(3)',y_min, y_max )# -1.0*pos_cal_sigma, 1.0*pos_cal_sigma)#y_min, y_max )
            #fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus(0)+pol3(3)', -sig_frac*pos_cal_sigma, sig_frac*pos_cal_sigma)#y_min, y_max )
            fit.SetParameter(1, -0.05)
            fit.SetParameter(2, 0.15)            

              

            fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
            fit.SetName(histo.GetTitle() + '_fit{}'.format(i))
            
	    if  calc_sigma:
		print(' using fitted value to init value of {}'.format(pos_cal_sigma))
		print(' using fitted value of MEAN value of {}'.format(calc_mean))
		fit.SetParameter(1, calc_mean)
		fit.SetParameter(2, pos_cal_sigma / 3 )
	        

            projec.Fit(fit,'R')#, -0.5, 0.5)

            if( fit.GetParameter(2) < 0 ):
                fit.SetParameter(1, fit.GetParameter(1))
                fit.SetParameter(2, np.abs(fit.GetParameter(2)))
                projec.Fit(fit,'R')#, -0.5, 0.5)

            
                
                

            #x_center = histo.GetXaxis().GetBinCenter(x_bin) 
            mean = fit.GetParameter(1)
            sig = fit.GetParameter(2)
            mean_3sig = mean - 3*sig
            mean_p3sig = mean + 3*sig
            means.append(mean)
            sigma.append(mean_3sig)
            psigma.append(mean_p3sig)
            x_centers.append(x_center)
            x_err.append(0)
            slices.append(projec)
        slices_sect.append(slices)
            
        mg = TMultiGraph()
        graph0 = TGraphErrors(len(x_centers), x_centers, means, x_err, x_err)
        graph0.SetName('g0_' + histo.GetName())
        fit_mean = TF1("fit_mean_s{}".format(ss),"pol3",x_centers[0], x_centers[-1])
	graph0.Fit("fit_mean_s{}".format(ss),"R")       

	graph1 = TGraphErrors(len(x_centers), x_centers, sigma, x_err, x_err)        
        graph1.SetName('g1_' + histo.GetName())
        fit_neg_sigma = TF1("fit_neg_sigma_s{}".format(ss),"pol3",x_centers[0], x_centers[-1])
        graph1.Fit("fit_neg_sigma_s{}".format(ss),"R")

        graph2 = TGraphErrors(len(x_centers), x_centers, psigma, x_err, x_err)
        graph2.SetName('g2_' + histo.GetName())
        fit_pos_sigma = TF1("fit_pos_sigma_s{}".format(ss),"pol3",x_centers[0], x_centers[-1])
        graph2.Fit("fit_pos_sigma_s{}".format(ss),"R")


	# write out calcualted fit parameters
	f_pos_out = open('{0}_pos_out_poly3_par_s{1}.txt'.format(particle_name, ss), 'w')	
	f_neg_out = open('{0}_neg_out_poly3_par_s{1}.txt'.format(particle_name, ss), 'w')
	f_mean_out = open('{0}_mean_out_poly3_par_s{1}.txt'.format(particle_name, ss), 'w')
	
	for par_i in range(0, fit_pos_sigma.GetNpar()):
		print(" writing out value of par {}".format(fit_pos_sigma.GetParameter(par_i)))
		f_pos_out.write("{} \n".format(fit_pos_sigma.GetParameter(par_i)))

	for par_i in range(0, fit_neg_sigma.GetNpar()):
		print(" writing out value of par {}".format(fit_neg_sigma.GetParameter(par_i)))
		f_neg_out.write("{} \n".format(fit_neg_sigma.GetParameter(par_i)))

	for par_i in range(0, fit_mean.GetNpar()):
		print(" writing out value of par mean {}".format(fit_mean.GetParameter(par_i)))
		f_mean_out.write("{} \n".format(fit_mean.GetParameter(par_i)))

	f_pos_out.close()	
	f_neg_out.close()
	f_mean_out.close()

        graph0.SetMarkerStyle(8)
        graph0.SetMarkerColor(kBlack)
        graph0.SetMarkerSize(1)
        
        graph1.SetMarkerStyle(8)
        graph1.SetMarkerColor(kBlack)
        graph1.SetMarkerSize(1)

        graph2.SetMarkerStyle(8)
        graph2.SetMarkerColor(kBlack)
        graph2.SetMarkerSize(1)
        
        mg.Add(graph0)#,"p")
        mg.Add(graph1)#,"p")
        mg.Add(graph2)#,"p")
        root_is_dumb.append(graph0)
        root_is_dumb.append(graph1)
        root_is_dumb.append(graph2)
        root_is_dumb.append(mg)
        l_mg.append(mg)
        #mg.Draw("AP")

        
    canvas.Clear()
    canvas.Divide(2,3)
    for ss in range(1,7):
        canvas.cd(ss)
        canvas.SetLogz()
    
        root_is_dumb.append(histos.get(title_formatter.format(ss), default_histo))
        
        htitle = title_formatter.format(ss)
        histos[htitle].Draw("colz")    
        l_mg[ss-1].Draw("P")

        if title:
            label.DrawLatex(0.1, 0.925, title)
            
        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)
            
        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)    
    
        root_is_dumb.append(label)


    canvas.Print(save_name)


    
    slice_can = TCanvas("s{}".format(0),"s{}".format(0),900,1800)
    slice_can.SetBatch(kTRUE)
    slice_can.Divide(3,9)
    
    slice_can.Print('{}['.format("slice_"+title_formatter.format(0)+save_name))
    
    for i in range(1,7):  

        ff = 1
        for pf in slices_sect[i-1]:
            slice_can.cd(ff)
            ff+=1
            pf.GetXaxis().SetRangeUser(-1.0, 1.0)
            pf.Draw()
            label.DrawLatex(0.5, 0.015,ytitle)
            label.DrawLatex(0.1, 0.925, 'Momentum Bin {}'.format(ff))
	                

        slice_can.Print("slice_"+title_formatter.format(0)+save_name)
    slice_can.Print('{}]'.format("slice_"+title_formatter.format(0)+save_name))
    
    root_is_dumb.append(label)
    #canvas.Print(save_name)
            
        

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

    plot_sector_page(can, histos, "ftof_deltat_evst_pr_s{}", lab, output_pdfname,
                particle_name = "pr", xtitle="Momentum (GeV)", ytitle="#Delta t (ns)", title="#Delta t vs Momentum FTOF Panel 1B, Proton",  log=True,
                     x_range=[6,23], y_range = [-0.3, 0.3], sig_frac = 1.1, calc_sigma=False)


    plot_sector_page(can, histos, "ftof_deltat_evst_kp_s{}", lab, output_pdfname,
                 particle_name = "kp", xtitle="Momentum (GeV)", ytitle="#Delta t (ns)", title="#Delta t vs Momentum FTOF Panel 1B, Kaon^{+}",  log=True,
                     x_range=[5,14], y_range=[-0.2,0.2], sig_frac=2, calc_sigma=False)


    can.Print('{}]'.format(output_pdfname))
