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

from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TMath, TFormula
from ROOT import TLegend, TLatex, TLine
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad, kBlack, kRed, kBlue, kMagenta, kPink
from array import array

import math
import sys

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


gStyle.SetOptStat(0000)


## method to get fits to phi and the S/B ratio
def plot_phi_mass_page(canvas, histos, title_formatter, label, save_name, plot_details,
                       xtitle=None, ytitle=None, title=None, fout=None):

    
    
    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    canvas.Divide(1,2)
    root_is_dumb.append(canvas)

    #fit form
    gaus = "[0]*TMath::Gaus(x,[1],[2])"
    bck = "[3]*(x>2*0.4937)*TMath::Power(x-2*0.4937, [4])*TMath::Exp(-[5]*(x - 2*0.4937) )"

    form = gaus + " + " + bck
    bckform = bck
    fitparm3 = plot_details['fitparms']  
    field = plot_details['field_setting']

    f_out = None
    if fout:
        f_out = open(fout,'w')

    for bb in range(1, len(plot_details['dimension'])+1):
        canvas.cd(bb)
        
        if 'fitparms_binned' in plot_details.keys():
            fitparm3 = plot_details['fitparms_binned'][bb-1]
        print(' fit parameters {}'.format(fitparm3))

        phim = histos[title_formatter.format(bb)]
        phim.SetLineColor(kBlack)
        phim.GetXaxis().SetRangeUser(0.95, 1.140)
        phim.Draw()
        print(" my function to fit phi is " + form)
        max_fit_range = 1.065
        phif = TF1("fitPhi{}".format(bb),form, 2*0.4937, max_fit_range)

        for pp in range(0,len(fitparm3)):
            print("par number %d, value of %f" %(pp, fitparm3[pp]) )
            phif.SetParameter(pp,fitparm3[pp])

        phif.SetNpx(10000)    
        phim.Fit("fitPhi{}".format(bb),"BR")


        fbck = TF1("bck",form, 2*0.4937, max_fit_range)
        for pp in range(3,len(fitparm3)):
            print("par number %d, value of %f" %(pp, phif.GetParameter(pp) ))
            fbck.SetParameter(pp, phif.GetParameter(pp))
            print(' fbck par {} is {}'.format(pp, fbck.GetParameter(pp)))
        #fbck.SetNpx(000)
        fbck.SetLineColor(kBlue)
        fbck.SetLineStyle(2)
        fbck.SetLineWidth(3)
        fbck.FixParameter(0,0)
        
        bin_width = phim.GetXaxis().GetBinWidth(1)
        print(' bin width to normalize {}'.format(bin_width))
        fsignal=TF1('fsignal',form, 0.98, max_fit_range)
        #fsignal.SetNpx(10000)
        for ii in range(0, len(fitparm3)):
            fsignal.SetParameter(ii,phif.GetParameter(ii))

        root_is_dumb.append(fsignal)
        root_is_dumb.append(fbck)
        root_is_dumb.append(phif)
        root_is_dumb.append(phim)

        ## get the phi signal with N - dN AND S - dS and max with (+) signs
        fit_sigma_err = phif.GetParError(2)
        fit_norm_err = phif.GetParError(0)
        fit_sigma = phif.GetParameter(2)
        fit_norm  = phif.GetParameter(0)

        sig_gaus = TF1('sig_gaus{}'.format(bb),'gaus(0)',0.98,max_fit_range)
        sig_gaus.SetParameter(0,phif.GetParameter(0))
        sig_gaus.SetParameter(1,phif.GetParameter(1))
        sig_gaus.SetParameter(2,phif.GetParameter(2))
        sig_gaus.SetLineColor(kMagenta)
        root_is_dumb.append(sig_gaus)
        
        fbck.Draw('same')
        sig_gaus.Draw('same')

        normC = np.sqrt( 2.0 * TMath.Pi())/bin_width
        d_nsignal = np.sqrt( pow( fit_norm_err*fit_sigma ,2) + pow( fit_norm*fit_sigma_err ,2)  ) * normC # use specific range, FX integrals -inf to inf.
        Nsignal = phif.GetParameter(0)*phif.GetParameter(2)*(TMath.Sqrt(2*TMath.Pi())/bin_width)
        label.DrawLatex(0.52, 0.82, 'Number of Phi Events: {0} #pm {1}'.format(int(Nsignal), int(d_nsignal) ) )
        label.DrawLatex(0.52, 0.76, '#mu = {0:6.5} GeV #pm {1:2.2f} MeV'.format(phif.GetParameter(1), phif.GetParError(1)*1000.0))
        label.DrawLatex(0.52, 0.7,  '#sigma = {0:2.2f} #pm {1:2.2f} MeV'.format(phif.GetParameter(2)*1000.0, phif.GetParError(2)*1000.0))
        label.DrawLatex(0.52, 0.64, 'FWHM = {0:2.2f} MeV'.format(2*TMath.Sqrt(2.0*TMath.Log(2.0))*phif.GetParameter(2)*1000.0))

        ## special calculation of S/B done here
        ## this is for method 1 in the thesis
        min_range_to_use = 1.0107
        max_range_to_use = 1.0287
        int_cal_bck = fbck.Integral(min_range_to_use, max_range_to_use)
        int_cal_sig = sig_gaus.Integral(min_range_to_use, max_range_to_use)
        sb_ratio_to_use = int_cal_sig/int_cal_bck
        print('defined signal range from {} - {} has S/B of {}'.format(min_range_to_use, max_range_to_use,sb_ratio_to_use))
        label.DrawLatex(0.52, 0.58, 'S/B in {0:.4f}'.format(min_range_to_use) +'< K^{ +}K^{ -}'+'<{0:0.4f}'.format(max_range_to_use)+' = {0:.4f}'.format(sb_ratio_to_use))

        label.DrawLatex(0.1, 0.925,'Invariant Mass of Charged Kaons from epK^{+}K^{-}')
        label.DrawLatex(0.65, 0.925,'Bin Range:{0:.4f}< {1} < {2:.4f}'.format(plot_details['binning'][bb-1], plot_details['varname'], plot_details['binning'][bb]))
        label.DrawLatex(0.5, 0.015,'Invariant Mass K^{+}K^{-} (GeV)')
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5,'Counts')
        label.SetTextAngle(0)
        
        if fout:
            # just write s/b ratio to file 
            f_out.write(str(sb_ratio_to_use) + ' \n')
        
        flab = TLatex()
        flab.SetTextColorAlpha(1, 0.376)
        flab.SetNDC()
        flab.DrawLatex(0.125, 0.847, field_setting)
        flab.SetTextFont(42)
        flab.SetTextSize(0.08)

        flab.SetTextSize(0.05)
        flab.SetTextColorAlpha(1, 1)
        root_is_dumb.append(label)
        root_is_dumb.append(flab)
                        
    if fout:
        f_out.close()
    canvas.Print(save_name)

#########################################################

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

    can = TCanvas('can', 'can', 1000, 900 )#800, 1100)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.04)
    lab.SetTextColor(1)

    chi2_lines=[0]

    can.Print('{}['.format(output_pdfname))

            
    field_setting=""
    if 'inbNout' in output_pdfname:
        field_setting="inbNoutb"
    elif 'inb' in output_pdfname:
        field_setting="inb"
    elif 'outb' in output_pdfname:
        field_setting="outb"
    else:
        field_setting=""

    q2_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_q2_'+field_setting+'.txt')
    xb_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_xb_'+field_setting+'.txt')
    t_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_t_'+field_setting+'.txt')
    w_bins=readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_limits_w_'+field_setting+'.txt')

    q2_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_q2_'+field_setting+'.txt')
    xb_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_xb_'+field_setting+'.txt')
    t_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_t_'+field_setting+'.txt')
    w_centers = readInBinningInfo('/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/phiasy/bin_centers_w_'+field_setting+'.txt')
    
    t_plot_detail={ 'fitparms': [ 1000 , 1.0197 , 0.00438 , 2452.97 , 0.339580 , 22],#0.1 ],
                    'field_setting':field_setting,
                    'dimension': t_centers,
                    'binning':t_bins,
                    'varname':'-t'
    }


    plot_phi_mass_page(can, histos, 't_bin_{}_mass_small_bin', lab, output_pdfname , t_plot_detail,
                       fout='t_bin_mass_sb_ratio_'+field_setting+'_final.txt')

    ################
    ## q2 binning 
    q2_plot_detail={ 'fitparms':  [ 350 , 1.0200 , 0.004207, 2452.97 , 0.6753 , 22 ], #[ 1000 , 1.01937 , 0.00457582 , 2452.97 , 0.339580 , 0.1 ],
                    'field_setting':field_setting,
                    'dimension': q2_centers,
                    'binning':q2_bins,
                    'varname':'Q^{2}'
    }


    plot_phi_mass_page(can, histos, 'q2_bin_{}_mass_small_bin', lab, output_pdfname , q2_plot_detail,
                       fout='q2_bin_mass_sb_ratio_'+field_setting+'_final.txt')

    ################
    ## xb binning 
    xb_plot_detail={ 'fitparms': [ 1000, 1.0198 , 0.0043 , 2452.97 , 0.339580 , 0.1 ],
                     'fitparms_binned': [[ 1700, 1.0198 , 0.0043 , 2452.97 , 0.339580 , 0.1 ],[ 1000, 1.0196 , 0.00407 , 2452.97 , 0.339580 , 0.1 ]],
                     'field_setting':field_setting,
                     'dimension': xb_centers,
                     'binning':xb_bins,
                     'varname':'X_{b}'
    }


    plot_phi_mass_page(can, histos, 'xb_bin_{}_mass_small_bin', lab, output_pdfname , xb_plot_detail,
                       fout='xb_bin_mass_sb_ratio_'+field_setting+'_final.txt')

    
    ################
    ## W binning 
    w_plot_detail={ 'fitparms': [ 3000 , 1.0195 , 0.0042 , 2452.97 , 0.339580 , 0.1 ], 
                    'field_setting':field_setting,
                    'dimension': w_centers,
                    'binning':w_bins,
                    'varname':'W'
    }


    plot_phi_mass_page(can, histos, 'w_bin_{}_mass_small_bin', lab, output_pdfname , w_plot_detail,
                       fout='w_bin_mass_sb_ratio_'+field_setting+'_final.txt')

    '''
    #################
    ## W binning but for low q2 < 2.75
    w_plot_detail_lowq2={ 'fitparms': [ 2700 , 1.02 , 0.0048, 2452.97 , 0.339580 , 0.1 ], 
                          'field_setting':field_setting,
                          'dimension': w_centers,
                          'binning':w_bins,
                          'varname':'W'
    }

    plot_phi_mass_page(can, histos, 'w_bin_{}_mass_small_bin_lowq2', lab, output_pdfname , w_plot_detail_lowq2,
                       fout='w_bin_mass_lowq2_sb_ratio_'+field_setting+'_final.txt')
    '''
    

    can.Print('{}]'.format(output_pdfname))
