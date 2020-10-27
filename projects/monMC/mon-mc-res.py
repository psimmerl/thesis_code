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
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas, TLegend,
                  gPad, gStyle, TLatex, TGraphErrors,
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
    
                              
    h_temp.Fit(h_name,"R") # add Q for quiet
    save_fits.append(fit_temp)
    
    temp_mean = fit_temp.GetParameter(1)
    temp_rms = fit_temp.GetParameter(2)
    
    return fit_temp#[temp_mean, temp_rms ]


def plot_page(canvas, histos, histo_title, label, save_name,
              xtitle=None, ytitle=None, title=None, log=False,
              fit=None, x_bin_step=None, x_range=None, y_range = None, y_fit_range=None, 
              fit_method=None, merge_bins_status=None, slice_fit_info = None, manual_gfit_par=None):
    
    canvas.Clear() 
    root_is_dumb=[]
    if isinstance(histos[histo_title], TH1F):
        histos[histo_title].Draw()

        if fit:
            fit_meth = TF1('fit'+histos[histo_title].GetTitle(),'gaus(0)',x_range[0], x_range[1])
            histos[histo_title].Fit('fit'+histos[histo_title].GetTitle(),'R')
            fit_meth.SetLineColor(kRed)
            fit_meth.Draw('same')
            label.DrawLatex(0.6, 0.85, 'Counts: {0}'.format(histos[histo_title].GetEntries()))
            label.DrawLatex(0.6, 0.8, '#mu: {0:6.4f}'.format(fit_meth.GetParameter(1)))
            label.DrawLatex(0.6, 0.75, '#sigma: {0:6.4f}'.format(fit_meth.GetParameter(2)))
            
            
            root_is_dumb.append(fit_meth)

    elif isinstance(histos[histo_title], TH2F):
        histos[histo_title].Draw('colz')
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

        
    canvas.Print(save_name)


def plot_page_overlap1D(canvas, histos, histo_titles, form, label, save_name,
                        xtitle=None, ytitle=None, title=None, log=False,
                        fit=None, x_bin_step=None, x_range=None, y_range = None, y_fit_range=None, 
                        fit_method=None, merge_bins_status=None, slice_fit_info = None, manual_gfit_par=None):
    
    canvas.Clear() 
    root_is_dumb=[]
    if isinstance(histos[histo_titles.format(form[0])], TH1F):
        histo_title = histo_titles.format(form[0])
        histo_title_smear = histo_titles.format(form[1])
        histos[histo_title].SetLineColor(kBlack)
        histos[histo_title_smear].SetLineColor(kRed)
        
        histos[histo_title].Draw()
        histos[histo_title_smear].Draw("same")

        leg = TLegend(0.7, 0.7, 0.9, 0.9)
        leg.AddEntry(histos[histo_title],'No Smear','l')
        leg.AddEntry(histos[histo_title_smear],'With Smear','l')
        leg.Draw('same')
        root_is_dumb.append(leg)

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

    ## 2d resolution for electron
    plot_page(can, histos, 'p_ele_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{el} (GeV)', ytitle='counts', title='FD:El, Final Phi, #Delta P_{el}')

    plot_page(can, histos, 'theta_ele_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' #theta_{el} (deg)', ytitle='counts', title='FD:El, Final Phi, #Delta #theta_{el}')

    plot_page(can, histos, 'phi_ele_res_sim_resolutions' , lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='counts', title='FD:El, Final Phi, #Delta #phi_{el}')

    plot_page(can, histos, 'p_ele_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{el} (GeV)', ytitle='#Delta P_{el} (GeV)', title='FD:El, Final Phi, #Delta P_{el} vs P_{el}')

    plot_page(can, histos, 'theta_ele_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{el} (GeV)', ytitle='#Delta P_{el} (GeV)', title='FD:El, Final Phi, #Delta P_{el} vs #Theta_{el}')

    plot_page(can, histos, 'theta_ele_delta_theta_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{el} (GeV)', ytitle='#Delta #Theta_{el} (deg)', title='FD:El, Final Phi, #Delta #Theta_{el} vs #Theta_{el}')

    # 2d resolution for proton
    plot_page(can, histos, 'p_pro_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{p} (GeV)', ytitle='counts', title='FD:Pr, Final Phi, #Delta P_{p}')
    plot_page(can, histos, 'theta_pro_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' #theta_{p} (deg)', ytitle='counts', title='FD:Pr, Final Phi, #Delta #theta_{p}')
    plot_page(can, histos, 'phi_pro_res_sim_resolutions' , lab, output_pdfname,
              xtitle='#phi_{p} (deg)', ytitle='counts', title='FD:Pr, Final Phi, #Delta #phi_{p}')

    plot_page(can, histos, 'p_pro_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{pr} (GeV)', ytitle='#Delta P_{pr} (GeV)', title='FD:Pr, Final Phi, #Delta P_{pr} vs P_{pr}')

    plot_page(can, histos, 'theta_pro_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{pr} (GeV)', ytitle='#Delta P_{pr} (GeV)', title='FD:Pr, Final Phi, #Delta P_{pr} vs #Theta_{pr}')

    plot_page(can, histos, 'theta_pro_delta_theta_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{pr} (GeV)', ytitle='#Delta #Theta_{pr} (deg)', title='FD:Pr, Final Phi, #Delta #Theta_{pr} vs #Theta_{pr}')

    # 2d resolution for K+
    plot_page(can, histos, 'p_kp_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{K^{ +}} (GeV)', ytitle='counts', title='FD:K^{ +}, Final Phi, #Delta P_{K^{ +}}')
    plot_page(can, histos, 'theta_kp_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' #theta_{K^{ +}} (deg)', ytitle='counts', title='FD:K^{ +}, Final Phi, #Delta #theta_{K^{ +}}')
    plot_page(can, histos, 'phi_kp_res_sim_resolutions' , lab, output_pdfname,
              xtitle='#phi_{K^{ +}} (deg)', ytitle='counts', title='FD:K^{ +}, Final Phi, #Delta #phi_{K^{ +}}')

    plot_page(can, histos, 'p_kp_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{K^{ +}} (GeV)', ytitle='#Delta P_{K^{ +}} (GeV)', title='FD: K^{ +}, Final Phi, #Delta P_{K^{ +}} vs P_{K^{ +}}')

    plot_page(can, histos, 'theta_kp_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{K^{ +}} (GeV)', ytitle='#Delta P_{K^{ +}} (GeV)', title='FD: K^{ +}, Final Phi, #Delta P_{K^{ +}} vs #Theta_{K^{ +}}')

    plot_page(can, histos, 'theta_kp_delta_theta_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{K^{ +}} (GeV)', ytitle='#Delta #Theta_{K^{ +}} (deg)', title='FD: K^{ +}, Final Phi, #Delta #Theta_{K^{ +}} vs #Theta_{K^{ +}}')

    # 2d resolution for K-
    plot_page(can, histos, 'p_km_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{K^{ -}} (GeV)', ytitle='counts', title='FD:K^{ -}, Final Phi, #Delta P_{K^{ -}}')
    plot_page(can, histos, 'theta_km_res_sim_resolutions' , lab, output_pdfname,
              xtitle=' #theta_{K^{ -}} (deg)', ytitle='counts', title='FD:K^{ -}, Final Phi, #Delta #theta_{K^{ -}}')
    plot_page(can, histos, 'phi_km_res_sim_resolutions' , lab, output_pdfname,
              xtitle='#phi_{K^{ -}} (deg)', ytitle='counts', title='FD:K^{ -}, Final Phi, #Delta #phi_{K^{ -}}')

    plot_page(can, histos, 'p_km_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle=' P_{K^{ -}} (GeV)', ytitle='#Delta P_{K^{ -}} (GeV)', title='FD: K^{ -}, Final Phi, #Delta P_{K^{ -}} vs P_{K^{ -}}')

    plot_page(can, histos, 'theta_km_delta_p_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{K^{ -}} (GeV)', ytitle='#Delta P_{K^{ -}} (GeV)', title='FD: K^{ -}, Final Phi, #Delta P_{K^{ -}} vs #Theta_{K^{ -}}')

    plot_page(can, histos, 'theta_km_delta_theta_sim_resolutions' , lab, output_pdfname,
              xtitle='#Theta_{K^{ -}} (GeV)', ytitle='#Delta #Theta_{K^{ -}} (deg)', title='FD: K^{ -}, Final Phi, #Delta #Theta_{K^{ -}} vs #Theta_{K^{ -}}')


    # 2d resolutionfor Q2, -t , xb 
    plot_page(can, histos, 'q2_delta_q2_sim_resolutions' , lab, output_pdfname,
              xtitle='Q^{2} (GeV^{2})', ytitle='#Delta Q^{2} (GeV^{2})', title='FD: Q^{2}, Final Phi, #Delta Q^{2} vs Q^{2}')

    plot_page(can, histos, 't_delta_t_sim_resolutions' , lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='#Delta -t (GeV^{2})', title='FD: -t, Final Phi, #Delta -t vs -t')

    plot_page(can, histos, 'xb_delta_xb_sim_resolutions' , lab, output_pdfname,
              xtitle='X_{b}', ytitle='#Delta X_{b}', title='FD: X_{b}, Final Phi, #Delta X_{b} vs X_{b}')


    ## single 1d for q2, -t, xb, phi trento

    plot_page(can, histos, 'delta_q2_sim_resolutions' , lab, output_pdfname,
              xtitle='#Delta Q^{2} (GeV^{2})', ytitle='counts', title='FD: Q^{2}, Final Phi, #Delta Q^{2}', fit=True,x_range=[-0.03, 0.03])

    plot_page(can, histos, 'delta_t_sim_resolutions' , lab, output_pdfname,
              xtitle='#Delta -t (GeV^{2})', ytitle='counts', title='FD: -t, Final Phi, #Delta -t', fit=True,x_range=[-0.02,0])

    plot_page(can, histos, 'delta_xb_sim_resolutions' , lab, output_pdfname,
              xtitle='#Delta X_{b}', ytitle='counts', title='FD: X_{b}, Final Phi, #Delta X_{b}', fit=True,x_range=[-0.01, 0.01])

    plot_page(can, histos, 'delta_phi_trento_sim_resolutions' , lab, output_pdfname,
              xtitle='#Delta #phi_{Trento}', ytitle='counts', title='FD: #phi_{Trento}, Final Phi, #Delta #phi_{Trento}', fit=True,x_range=[-3, 3])




    #####################################################################################################################################################################
    add_text_page(can, lab, 'smearing final state particles', output_pdfname)

    
    plot_page(can, histos, 'p_ele_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{el} (GeV)', ytitle='counts', title='FD:El, Final Phi,Smear, #Delta P_{el}')

    plot_page(can, histos, 'theta_ele_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' #theta_{el} (deg)', ytitle='counts', title='FD:El, Final Phi,Smear, #Delta #theta_{el}')

    plot_page(can, histos, 'phi_ele_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='counts', title='FD:El, Final Phi,Smear, #Delta #phi_{el}')

    plot_page(can, histos, 'p_ele_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{el} (GeV)', ytitle='#Delta P_{el} (GeV)', title='FD:El, Final Phi,Smear, #Delta P_{el} vs P_{el}')

    plot_page(can, histos, 'theta_ele_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{el} (GeV)', ytitle='#Delta P_{el} (GeV)', title='FD:El, Final Phi,Smear #Delta P_{el} vs #Theta_{el}')

    plot_page(can, histos, 'theta_ele_delta_theta_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{el} (GeV)', ytitle='#Delta #Theta_{el} (deg)', title='FD:El, Final Phi,Smear #Delta #Theta_{el} vs #Theta_{el}')

    ## proton smeared
    plot_page(can, histos, 'p_pro_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{p} (GeV)', ytitle='counts', title='FD:Pr, Final Phi,Smear, #Delta P_{p}')
    plot_page(can, histos, 'theta_pro_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' #theta_{p} (deg)', ytitle='counts', title='FD:Pr, Final Phi,Smear, #Delta #theta_{p}')
    plot_page(can, histos, 'phi_pro_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#phi_{p} (deg)', ytitle='counts', title='FD:Pr, Final Phi,Smear, #Delta #phi_{p}')


    plot_page(can, histos, 'p_pro_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{pr} (GeV)', ytitle='#Delta P_{pr} (GeV)', title='FD:Pr, Final Phi,Smear #Delta P_{pr} vs P_{pr}')

    plot_page(can, histos, 'theta_pro_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{pr} (GeV)', ytitle='#Delta P_{pr} (GeV)', title='FD:Pr, Final Phi,Smear #Delta P_{pr} vs #Theta_{pr}')

    plot_page(can, histos, 'theta_pro_delta_theta_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{pr} (GeV)', ytitle='#Delta #Theta_{pr} (deg)', title='FD:Pr, Final Phi,Smear #Delta #Theta_{pr} vs #Theta_{pr}')


    # 2d resolution for K+ SMEAR
    plot_page(can, histos, 'p_kp_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{K^{ +}} (GeV)', ytitle='counts', title='FD:K^{ +}, Final Phi,Smear, #Delta P_{K^{ +}}')
    plot_page(can, histos, 'theta_kp_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' #theta_{K^{ +}} (deg)', ytitle='counts', title='FD:K^{ +}, Final Phi,Smear, #Delta #theta_{K^{ +}}')
    plot_page(can, histos, 'phi_kp_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#phi_{K^{ +}} (deg)', ytitle='counts', title='FD:K^{ +}, Final Phi,Smear, #Delta #phi_{K^{ +}}')

    plot_page(can, histos, 'p_kp_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{K^{ +}} (GeV)', ytitle='#Delta P_{K^{ +}} (GeV)', title='FD: K^{ +},Smear Final Phi, #Delta P_{K^{ +}} vs P_{K^{ +}}')

    plot_page(can, histos, 'theta_kp_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{K^{ +}} (GeV)', ytitle='#Delta P_{K^{ +}} (GeV)', title='FD: K^{ +},Smear Final Phi, #Delta P_{K^{ +}} vs #Theta_{K^{ +}}')

    plot_page(can, histos, 'theta_kp_delta_theta_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{K^{ +}} (GeV)', ytitle='#Delta #Theta_{K^{ +}} (deg)', title='FD: K^{ +},Smear Final Phi, #Delta #Theta_{K^{ +}} vs #Theta_{K^{ +}}')

    # 2d resolution for K- SMEAR
    plot_page(can, histos, 'p_km_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{K^{ -}} (GeV)', ytitle='counts', title='FD:K^{ -}, Final Phi,Smear, #Delta P_{K^{ -}}')
    plot_page(can, histos, 'theta_km_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' #theta_{K^{ -}} (deg)', ytitle='counts', title='FD:K^{ -}, Final Phi,Smear, #Delta #theta_{K^{ -}}')
    plot_page(can, histos, 'phi_km_res_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#phi_{K^{ -}} (deg)', ytitle='counts', title='FD:K^{ -}, Final Phi,Smear, #Delta #phi_{K^{ -}}')

    plot_page(can, histos, 'p_km_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle=' P_{K^{ -}} (GeV)', ytitle='#Delta P_{K^{ -}} (GeV)', title='FD: K^{ -}, Final Phi,Smear #Delta P_{K^{ -}} vs P_{K^{ -}}')

    plot_page(can, histos, 'theta_km_delta_p_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{K^{ -}} (GeV)', ytitle='#Delta P_{K^{ -}} (GeV)', title='FD: K^{ -}, Final Phi,Smear #Delta P_{K^{ -}} vs #Theta_{K^{ -}}')

    plot_page(can, histos, 'theta_km_delta_theta_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Theta_{K^{ -}} (GeV)', ytitle='#Delta #Theta_{K^{ -}} (deg)', title='FD: K^{ -}, Final Phi,Smear #Delta #Theta_{K^{ -}} vs #Theta_{K^{ -}}')

    # 2d resolutionfor Q2, -t , xb 
    plot_page(can, histos, 'q2_delta_q2_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='Q^{2} (GeV^{2})', ytitle='#Delta Q^{2} (GeV^{2})', title='FD: Q^{2}, Final Phi, Smear, #Delta Q^{2} vs Q^{2}')

    plot_page(can, histos, 't_delta_t_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='#Delta -t (GeV^{2})', title='FD: -t, Final Phi,Smear #Delta -t vs -t')

    plot_page(can, histos, 'xb_delta_xb_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='X_{b}', ytitle='#Delta X_{b}', title='FD: X_{b}, Final Phi,Smear #Delta X_{b} vs X_{b}')


    # 1d SMEARed Q2, -t, xb
    plot_page(can, histos, 'delta_q2_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Delta Q^{2} (GeV^{2})', ytitle='counts', title='FD: Q^{2}, Final Phi,Smear, #Delta Q^{2}', fit=True,x_range=[-0.03, 0.03])

    plot_page(can, histos, 'delta_t_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Delta -t (GeV^{2})', ytitle='counts', title='FD: -t, Final Phi,Smear, #Delta -t', fit=True,x_range=[-0.02,0])

    plot_page(can, histos, 'delta_xb_sim_resolutions_with_smear' , lab, output_pdfname,
              xtitle='#Delta X_{b}', ytitle='counts', title='FD: X_{b}, Final Phi, Smear, #Delta X_{b}', fit=True,x_range=[-0.01, 0.01])
    

    #####################################################################################################################################################################
    # compare with and without smear 
    add_text_page(can, lab, 'compare w/ w/out smear', output_pdfname)
    plot_page_overlap1D(can, histos, 'p_ele_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' P_{el} (GeV)', ytitle='counts', title='FD:El, Final Phi,Comp., #Delta P_{el}')

    plot_page_overlap1D(can, histos, 'theta_ele_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' #theta_{el} (deg)', ytitle='counts', title='FD:El, Final Phi,Comp, #Delta #theta_{el}')

    plot_page_overlap1D(can, histos, 'phi_ele_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle='#phi_{el} (deg)', ytitle='counts', title='FD:El, Final Phi,Comp, #Delta #phi_{el}')


    # compare proton 
    plot_page_overlap1D(can, histos, 'p_pro_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' P_{p} (GeV)', ytitle='counts', title='FD:Pr, Final Phi,Comp., #Delta P_{p}')

    plot_page_overlap1D(can, histos, 'theta_pro_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' #theta_{p} (deg)', ytitle='counts', title='FD:Pr, Final Phi,Comp, #Delta #theta_{p}')

    plot_page_overlap1D(can, histos, 'phi_pro_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle='#phi_{p} (deg)', ytitle='counts', title='FD:Pr, Final Phi,Comp, #Delta #phi_{p}')


    # compare kaons plus
    plot_page_overlap1D(can, histos, 'p_kp_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' P_{K^{ +}} (GeV)', ytitle='counts', title='FD:K^{ +}, Final Phi,Comp, #Delta P_{K^{ +}}')

    plot_page_overlap1D(can, histos, 'theta_kp_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' #theta_{K^{ +}} (deg)', ytitle='counts', title='FD:K^{ +}, Final Phi,Comp, #Delta #theta_{K^{ +}}')

    plot_page_overlap1D(can, histos, 'phi_kp_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle='#phi_{K^{ +}} (deg)', ytitle='counts', title='FD:K^{ +}, Final Phi,Smear, #Delta #phi_{K^{ +}}')    

    # compare kaons Minus
    plot_page_overlap1D(can, histos, 'p_km_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' P_{K^{ -}} (GeV)', ytitle='counts', title='FD:K^{ -}, Final Phi,Comp, #Delta P_{K^{ -}}')

    plot_page_overlap1D(can, histos, 'theta_km_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle=' #theta_{K^{ -}} (deg)', ytitle='counts', title='FD:K^{ -}, Final Phi,Comp, #Delta #theta_{K^{ -}}')

    plot_page_overlap1D(can, histos, 'phi_km_res_sim_resolutions{}', ['','_with_smear'], lab, output_pdfname,
                        xtitle='#phi_{K^{ -}} (deg)', ytitle='counts', title='FD:K^{ -}, Final Phi,Smear, #Delta #phi_{K^{ -}}')    

    can.Print('{}]'.format(output_pdfname))
    
