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


def plot_2d_sidebyside(canvas,  histos_inb, histos_outb, title_formatter, label, save_name, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, hlines=None, 
                       logy=None, set_marker_style=None, scale=None, shades=None):

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
        label.DrawLatex(0.1, 0.925, 'PIDTYPE1:'+title)

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
        label.DrawLatex(0.1, 0.925, 'PIDTYPE2:'+title)

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
            label.DrawLatex(0.5, 0.025, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)


        leg = TLegend(0.7, 0.7, 0.9, 0.9)
        leg.AddEntry(histos_outb[title_formatter],'PIDTYPE2','p')
        leg.AddEntry(histos_inb[title_formatter],'PIDTYPE1','p')
        leg.Draw('same')
        root_is_dumb.append(leg)

        
    canvas.Print(save_name)

def add_text_page(can, label, text, save_name):
    """ Write some text onto a page. """
    can.Clear()
    can.cd(1)
    label.DrawLatex(0.1, 0.5, text)
    can.Print(save_name)
    
    
if __name__ == '__main__':


    # inb referes to pidtype1 
    # out referes to pitdtyp2
    # pidtype 1 is just the first file 
    # can compare any variation of pids two at a time
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-ipid1',
        '--input_inb_file',
        required=True
    )
    ap.add_argument(
        '-ipid2',
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

    can = TCanvas('can', 'can', 800, 1100)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    kpkm_mass_limit=[1.01, 1.02]

    add_text_page(can, lab, "COMPARE PIDS", output_pdfname)
    
    ## look at the events from inbending and outbending before making too many cuts
    add_text_page(can, lab, "COMP-El Before Final Sel.", output_pdfname)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_ele_theta_pass_all', lab,  output_pdfname,
                       xtitle='P_{el} (GeV)', ytitle='#Theta_{el} (deg)', title='Electron (FD), #Theta vs P, Pass All')
    

    plot_2d_sidebyside(can, histos_inb, histos_outb,  'phi_ele_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='#Theta_{el} (deg)', title='Electron (FD), #Theta vs #phi, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_ele_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{el} (deg)', ytitle='Vz_{el} (cm)', title='Electron (FD), #phi vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'theta_ele_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{el} (deg)', ytitle='Vz_{el} (cm)', title='Electron (FD), #Theta vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb,  'p_ele_vz_pass_all', lab, output_pdfname,
              xtitle='P_{el} (GeV)', ytitle='Vz_{el} (cm)', title='Electron (FD), P vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_ele_dvz_elepro_pass_all', lab, output_pdfname,
              xtitle='P (GeV)', ytitle='#Delta Vz_{el} - Vz_{pro} (cm)', title='Electron (FD), P_{el} vs #Delta Vz_{el} - Vz_{pro}, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'xb_q2_pass_all', lab, output_pdfname,
              xtitle='Xb', ytitle='Q^{2} (GeV^{2})', title='(FD), Q^{2} vs Xb, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb,  'w_q2_pass_all', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='Q^{2} (GeV^{2})', title='(FD), Q^{2} vs W, Pass All')
    
   


    ## general proton kineatics
    add_text_page(can, lab, "COMP- Pr Before Final Sel.", output_pdfname)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_theta_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='#Theta_{pr} (deg)', title='Proton (FD), #Theta vs P, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb,  'phi_pro_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{pr} (deg)', ytitle='#Theta_{pr} (deg)', title='Proton (FD), #Theta vs #phi, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_pro_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{pr} (deg)', ytitle='Vz_{pr} (cm)', title='Proton (FD), #phi vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'theta_pro_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{pr} (deg)', ytitle='Vz_{pr} (cm)', title='Proton (FD), #Theta vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_vz_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='Vz_{pr} (cm)', title='Proton (FD), P vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_dvz_elepro_pass_all', lab, output_pdfname,
              xtitle='P_{pr} (GeV)', ytitle='#Delta Vz (cm)', title='Proton (FD), P vs #Delta Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{pro} (GeV)', ytitle='#beta_{pr}', title='Proton (FD), #beta vs P, Pass All')

    chi2_lines=[-6,-6]
    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_chi2_pass_all', lab, output_pdfname,
                       xtitle='P_{pro} (GeV)', ytitle='#chi^{2}', title='Proton (FD), #chi^{2} vs P, Pass All', hlines=[-6,6])

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'beta_pro_chi2_pass_all', lab, output_pdfname,
                       xtitle='#beta_{pro} (GeV)', ytitle='#chi^{2}', title='Proton (FD), #chi^{2} vs #beta_{pro}, Pass All',hlines=[-6,6])
    


    ## general kaon Plus 
    add_text_page(can, lab, "COMP -KaonP Before Final Sel.", output_pdfname)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_kp_theta_pass_all', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='#Theta_{kp} (deg)', title='Proton (FD), #Theta vs P, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_kp_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{kp} (deg)', ytitle='#Theta_{kp} (deg)', title='K^{ +} (FD), #Theta vs #phi, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_kp_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{kp} (deg)', ytitle='Vz_{kp} (cm)', title='K^{ +} (FD), #phi vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'theta_kp_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{kp} (deg)', ytitle='Vz_{kp} (cm)', title='K^{ +} (FD), #Theta vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_kp_vz_pass_all', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='Vz_{kp} (cm)', title='K^{ +} (FD), P vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_kp_dvz_kpkm_pass_all', lab, output_pdfname,
              xtitle='P_{kp} (GeV)', ytitle='#Delta Vz (cm)', title='K^{ +} (FD), P vs #Delta Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_kp_beta_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#beta_{K^{+}}', title='K^{ +}(FD), #beta vs P, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb,  'p_kp_chi2_pass_all', lab, output_pdfname,
                     xtitle='P_{K^{+}} (GeV)', ytitle='#chi^{2}', title='K^{ +} (FD), #chi^{2} vs P, Pass All', hlines=chi2_lines)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'beta_kp_chi2_pass_all', lab, output_pdfname,
                     xtitle='#beta_{K^{+}}', ytitle='#chi^{2}', title='K^{ +} (FD), #chi^{2} vs  #beta_{K^{+}}, Pass All', hlines=chi2_lines)


    ## general kaon Minus
    add_text_page(can, lab, "COMP -KaonM Before Final Sel.", output_pdfname)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_km_theta_pass_all', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='#Theta_{km} (deg)', title='Proton (FD), #Theta vs P, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_km_theta_pass_all', lab, output_pdfname,
              xtitle='#phi_{km} (deg)', ytitle='#Theta_{km} (deg)', title='K^{ -} (FD), #Theta vs #phi, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_km_vz_pass_all', lab, output_pdfname,
              xtitle='#phi_{km} (deg)', ytitle='Vz_{km} (cm)', title='K^{ -} (FD), #phi vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'theta_km_vz_pass_all', lab, output_pdfname,
              xtitle='#theta_{km} (deg)', ytitle='Vz_{km} (cm)', title='K^{ -} (FD), #Theta vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_km_vz_pass_all', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='Vz_{km} (cm)', title='K^{ -} (FD), P vs Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_km_dvz_kpkm_pass_all', lab, output_pdfname,
              xtitle='P_{km} (GeV)', ytitle='#Delta Vz (cm)', title='K^{ -} (FD), P vs #Delta Vz, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb,  'p_km_beta_pass_all', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#beta_{K^{-}}', title='K^{ -}(FD), #beta vs P, Pass All')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_km_chi2_pass_all', lab, output_pdfname,
              xtitle='P_{K^{-}} (GeV)', ytitle='#chi^{2}', title='K^{ -} (FD), #chi^{2} vs P, Pass All', hlines=chi2_lines)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'beta_km_chi2_pass_all', lab, output_pdfname,
              xtitle='#beta_{K^{-}}', ytitle='#chi^{2}', title='K^{ -} (FD), #chi^{2} vs  #beta_{K^{ -}}, Pass All', hlines=chi2_lines)

    
    add_text_page(can, lab, "COMP - Masses Pass All+Add.", output_pdfname)

    ###############
    ## general physics kinematics
    plot_2d_sidebyside(can, histos_inb, histos_outb,'im_kpkm_mm_ep_pass_all_pass_additional', lab, output_pdfname,
              xtitle=' MM epX (GeV)', ytitle=' I.M K^{ +}K^{ -} (GeV)', title='(FD), I.M K^{ +}K^{ -} vs MM epX, Pass All Cuts')

    plot_2d_sidebyside(can, histos_inb, histos_outb,'mm_epkp_mm_epkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle=' MM epK^{ -}X (GeV)', ytitle='MM epK^{ +}X (GeV)', title='(FD),MM epK^{ +}X vs MM epK^{ -}X, Pass All Cuts')
    
    
    plot_2d_sidebyside(can, histos_inb, histos_outb,'im_kpkm_im_pkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle='I.M. PrK^{ -} (GeV)', ytitle='I.M K^{ +}K^{ -} (GeV)', title='(FD), I.M K^{ +}K^{ -} vs I.M. PrK^{ -}, Pass All Cuts')

    plot_2d_sidebyside(can, histos_inb, histos_outb,'im_pkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle='I.M. PrK^{ -} (GeV)', ytitle='counts', title='(FD), I.M. PrK^{ -}, Pass All Cuts', shades=[1.5, 1.58, 1.78, 1.9] )

    plot_2d_sidebyside(can, histos_inb, histos_outb,'im_pkm_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='I.M. PrK^{ -} (GeV)', ytitle='counts', title='(FD), I.M. PrK^{ -}, PA, Pass #phi', shades=[1.5, 1.58, 1.78, 1.9] )

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phitrento_imkpkm_pass_all_pass_additional', lab, output_pdfname,
              xtitle='#phi_{trento} (deg)', ytitle='I.M. K^{ +}K^{ -} (GeV)', title='Invariant Mass K^{ +}K^{ -} vs #phi_{trento}')


    ###############
    ## general momentum resolution from reconstructed particles (not simulation resolution) of final final state.
    add_text_page(can, lab, "COMP-Calc. #Delta P PassAll #phi Peak", output_pdfname)

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'res_ekpkmX_mntm_vs_p_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{pro} (GeV)', ytitle='#Delta P_{pro} (GeV)', title='(FD),#Delta P_{pro} vs P_{pro}, Pass, #phi peak')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'res_epkpX_mntm_vs_p_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{K^{ -}} (GeV)', ytitle='#Delta P_{K^{ -}} (GeV)', title='(FD),#Delta P_{K^{ -}} vs P_{K^{ -}},Pass, #phi peak')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'res_epkmX_mntm_vs_p_pass_all_pass_additional_pass_phi_mass', lab, output_pdfname,
              xtitle='P_{K^{ +}} (GeV)', ytitle='#Delta P_{K^{ +}} (GeV)', title='(FD),#Delta P_{K^{ +}} vs P_{K^{ +}},Pass, #phi peak')

    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################
    #######################################################################################################################

    ## compare the pid types by overlapping the distributions
        #electron
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC PIDType1/PIDType2, P_{e} Final #phi (FD) Events' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC PIDType1/PIDType2, #Theta_{e} Final #phi (FD) Events' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_ele_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='REC PIDType1/PIDType2, #Phi_{e} Final #phi (FD) Events' , xtitle='#Phi_{el} (deg)', ytitle='counts', set_marker_style=True, scale=True)
    '''
    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_ele_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{el} vs P_{el}, PassAll, '+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{el} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_ele_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{el} vs #phi_{el}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{el} (deg)', ytitle='#theta (deg)')
    '''

    #proton
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, P_{p} Final #phi (FD) Events' , xtitle='P_{p} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #Theta_{p} Final #phi (FD) Events' , xtitle='#Theta_{p} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #Phi_{p} Final #phi (FD) Events' , xtitle='#Phi_{p} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    '''
    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_pro_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{p} vs P_{p}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{p} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_pro_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{p} vs #phi_{p}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{p} (deg)', ytitle='#theta (deg)')
    '''


    #kaons
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, P_{K^{ +}} Final #phi (FD) Events' , xtitle='P_{K^{ +}} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #Theta_{K^{ +}} Final #phi (FD) Events' , xtitle='#Theta_{K^{ +}} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #Phi_{K^{ +}} Final #phi (FD) Events' , xtitle='#Phi_{K^{ +}} (deg)', ytitle='counts', set_marker_style=True, scale=True)
    '''
    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_kp_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ +}} vs P_{K^{ +}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{K^{ +}} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_kp_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ +}} vs #phi_{K^{ +}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{K^{ +}} (deg)', ytitle='#theta (deg)')
    '''


    #kaonM
    plot_overlap_1d( can, histos_inb, histos_outb, 'p_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, P_{K^{ -}} Final #phi (FD) Events' , xtitle='P_{K^{ -}} (GeV)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'theta_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #Theta_{K^{ -}} Final #phi (FD) Events' , xtitle='#Theta_{K^{ -}} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phi_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #Phi_{K^{ -}} Final #phi (FD) Events' , xtitle='#Phi_{K^{ -}} (deg)', ytitle='counts', set_marker_style=True, scale=True)
    '''
    plot_2d_sidebyside(can, histos_inb, histos_outb, 'p_km_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ -}} vs P_{K^{ -}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='P_{K^{ -}} (GeV)', ytitle='#theta (deg)')

    plot_2d_sidebyside(can, histos_inb, histos_outb, 'phi_km_theta_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                       title='#theta_{K^{ -}} vs #phi_{K^{ -}}, PassAll,'+str(kpkm_mass_limit[0])+'<Mass(K^{ +}K^{ -})<'+str(kpkm_mass_limit[1]) , xtitle='#phi_{K^{ -}} (deg)', ytitle='#theta (deg)')
    '''

    # variables that are used for events selection
    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_pro_pass_all', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #theta_{Copl} Pr. Pass Me, MM2' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True, vlines =[9,9])

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_kp_pass_all', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #theta_{Copl} K^{ +}. Pass Me, MM2' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True, vlines =[9,9])

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_km_pass_all', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #theta_{Copl} K^{ -} Pass Me, MM2' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True, vlines =[9,9])


    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_pro_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #theta_{Copl} Pr. Final #phi (FD) Events' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_kp_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #theta_{Copl} K^{ +}. Final #phi (FD) Events' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'cpl_km_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #theta_{Copl} K^{ -} Final #phi (FD) Events' , xtitle='#theta_{Copl} (deg)', ytitle='counts', set_marker_style=True, scale=True)



    ##### now look at q2, xb, -t    
    ## dont scale to compare overall reduction if there is one
    plot_overlap_1d( can, histos_inb, histos_outb, 'w_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, W Final #phi (FD) Events' , xtitle='W (GeV)', ytitle='counts', set_marker_style=True, vlines=[2])

    plot_overlap_1d( can, histos_inb, histos_outb, 'q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV^{2})', ytitle='counts', set_marker_style=True, vlines=[1])

    plot_overlap_1d( can, histos_inb, histos_outb, 'xb_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, X_{b} Final #phi (FD) Events' , xtitle='X_{b}', ytitle='counts', set_marker_style=True)
        
    plot_overlap_1d( can, histos_inb, histos_outb, 't_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, -t Final #phi (FD) Events' , xtitle='-t (GeV^{2})', ytitle='counts', set_marker_style=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phitrento_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #phi_{trento} Final #phi (FD) Events' , xtitle='#phi_{trento} (deg)', ytitle='counts', set_marker_style=True)



    plot_overlap_1d( can, histos_inb, histos_outb, 'w_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, W Final #phi (FD) Events' , xtitle='W (GeV)', ytitle='counts', set_marker_style=True, scale=True, vlines=[2])

    plot_overlap_1d( can, histos_inb, histos_outb, 'q2_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True, vlines=[1])

    plot_overlap_1d( can, histos_inb, histos_outb, 'xb_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, X_{b} Final #phi (FD) Events' , xtitle='X_{b}', ytitle='counts', set_marker_style=True, scale=True)
        
    plot_overlap_1d( can, histos_inb, histos_outb, 't_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, -t Final #phi (FD) Events' , xtitle='-t (GeV^{2})', ytitle='counts', set_marker_style=True, scale=True)

    plot_overlap_1d( can, histos_inb, histos_outb, 'phitrento_pass_all_pass_additional_pass_phi_mass', output_pdfname, lab, 
                     title='PIDType1/PIDType2, #phi_{trento} Final #phi (FD) Events' , xtitle='#phi_{trento} (deg)', ytitle='counts', set_marker_style=True, scale=True)




    can.Print('{}]'.format(output_pdfname))
