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


def plot_acceptance( canvas, histos, title_formatter, hnames, save_name, label, 
                     title=None, xtitle=None, ytitle=None , x_range=None):

    canvas.Clear()
    canvas.Divide(1,1)

    root_is_dumb = []
    
    if isinstance(histos.get(title_formatter.format(hnames[0]), default_histo), TH1F) and isinstance(histos.get(title_formatter.format(hnames[1]), default_histo), TH1F):
        h1 = histos.get(title_formatter.format(hnames[0])).Clone(title_formatter.format(hnames[0] + '_clone1'))
        h2 = histos.get(title_formatter.format(hnames[1])).Clone(title_formatter.format(hnames[0] + '_clone1'))
        #h1.RebinX(10)
        #h2.RebinX(10)

        opt=''
        bin_center = array('d')
        bin_count = array('d')
        bin_err = array('d')
        zeros = array('d')

        for bb in range(1, h1.GetNbinsX()):
             
            h1_bc = h1.GetBinCenter(bb)
            h1_count = h1.GetBinContent(bb)
            h2_count = h2.GetBinContent(bb)
            if x_range and h1_bc > x_range[1]:
                continue
                
            bin_center.append(h1_bc)
            zeros.append(0.0)
            if h1_count > 0:
                bin_count.append( h2_count / h1_count)                            
                bin_error = 1.0/h1_count * np.sqrt( h2_count*h2_count / h1_count + h2_count)
                bin_err.append(bin_error)
            else:
                bin_count.append( 0.0 )
                bin_err.append(0)

        graph = TGraphErrors(len(bin_center), bin_center, bin_count, zeros, bin_err)
        graph.SetName('accp'+title_formatter.format(hnames[0]))
        graph.SetTitle('')
        graph.SetMarkerStyle(8)
        graph.SetMarkerColor(2)
        graph.SetMarkerSize(1)
        
        can.cd(1)
        graph.Draw("AP")
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)
        

    canvas.Print(save_name)
        
        
        

def plot_overlap_1d( canvas, histos, title_formatter, hnames, save_name, label, 
                     title=None, xtitle=None, ytitle=None, 
                     hline=None, vlines=None, logy=None, set_marker_style=None, scale=None, name=None):

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
            histos[title_formatter.format(hnames[0])].SetLineColor(kRed)
            histos[title_formatter.format(hnames[1])].SetLineColor(kBlack)
            

        if logy:
            print(' setting log Y ')
            gPad.SetLogy() 
            histos.get(title_formatter.format(hnames[0]), default_histo).SetMinimum(0.1)

        if scale:            
            histos[title_formatter.format(hnames[0])].Scale(  histos[title_formatter.format(hnames[1])].GetMaximum()/histos[title_formatter.format(hnames[0])].GetMaximum())


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

        leg = TLegend(0.7, 0.7, 0.9, 0.9)
        leg.AddEntry(histos[title_formatter.format(hnames[0])],name[0],'p')
        leg.AddEntry(histos[title_formatter.format(hnames[1])],name[1],'p')
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

    can = TCanvas('can', 'can', 1000, 1000)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    
    ## first compare the final state selected phi reconstructed phi events in GEMC to generated distribution
    # electrons
    
    
    #theta_pro_sim_rec_all->Divide(theta_pro_sim_gen)
    plot_overlap_1d( can, histos, 'p_ele_{}',['sim_gen','sim_rec_all'], output_pdfname, lab,
                     title='Gen. vs Rec, P_{e} All Ele Rec.' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, logy=True, name=['GEN','All REC'])

    plot_acceptance( can, histos, 'p_ele_{}',['sim_gen','sim_rec_all'], output_pdfname, lab,
                     title='Gen vs Rec, P_{e} All Ele Rec.' , xtitle='P_{el} (GeV)', ytitle='Acceptance')

    plot_overlap_1d( can, histos, 'p_ele_{}',['sim_gen_ele_present','sim_rec_all'], output_pdfname, lab,
                     title='Gen&Rec vs Rec, P_{e} All Ele Rec.' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, logy=True,name=['REC&GEN','All REC'])

    plot_overlap_1d( can, histos, 'theta_ele_{}',['sim_gen','sim_rec_all'], output_pdfname, lab,
                     title='Gen. vs Rec, #Theta_{e} All Ele Rec.' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, logy=True, name=['GEN','All REC'])

    plot_acceptance( can, histos,'theta_ele_{}',['sim_gen','sim_rec_all'], output_pdfname, lab,
                     title='Accep, #Theta_{e} All Ele. Rec' , xtitle='#Theta_{el} (deg)', ytitle='Acceptance')

    plot_overlap_1d( can, histos, 'theta_ele_{}',['sim_gen_ele_present','sim_rec_all'], output_pdfname, lab,
                     title='Gen&Rec vs Rec, #Theta_{e}' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, logy=True, name=['REC&GEN','All REC'])

    plot_overlap_1d( can, histos, 'theta_ele_{}',['sim_gen','sim_gen_ele_present'], output_pdfname, lab,
                     title='Gen&Rec, #Theta_{e}' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, logy=True, name=['REC&GEN','All REC'])


    plot_overlap_1d( can, histos, 'phi_ele_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Theta_{e}' , xtitle='#Phi_{el} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    plot_acceptance( can, histos, 'phi_ele_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #phi_{e}  All Ele Rec' , xtitle='#phi_{el} (deg)', ytitle='Acceptance', x_range=[-180,180])

    plot_overlap_1d( can, histos, 'phi_ele_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Theta_{e} Final #phi (FD) Events' , xtitle='#Phi_{el} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','REC'])



    #plot_acceptance( can, histos, 'phi_ele_{}', ['sim_gen_ele_present','sim_rec_all'], output_pdfname, lab, 
    #                 title='Accep, #phi_{e} Gen and Rec' , xtitle='#phi_{el} (deg)', ytitle='Acceptance', x_range=[-180,180])

    #####################
    # proton

    
    plot_overlap_1d( can, histos, 'p_pro_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, P_{pr},Trig Ele, All Pro Rec' , xtitle='P_{pr} (GeV)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    plot_acceptance( can, histos, 'p_pro_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, P_{pr}, Trig Ele, All Pro Rec' , xtitle='P_{Pr} (GeV)', ytitle='Acceptance', x_range=[0,3])

    
    plot_overlap_1d( can, histos, 'theta_pro_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Theta_{pr},Trig Ele, All Pro Rec' , xtitle='#Theta_{pr} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    
    plot_overlap_1d( can, histos, 'theta_pro_{}', ['sim_gen_pro_present','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. and Rec, #Theta_{pr},Trig Ele, All Pro Rec' , xtitle='#Theta_{pr} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['REC&GEN','All REC'])

    plot_acceptance( can, histos, 'theta_pro_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #Theta_{pr}, Trig Ele, All Pro Rec' , xtitle='#Theta_{Pr} (deg)', ytitle='Acceptance', x_range=[0,35])


    plot_overlap_1d( can, histos, 'phi_pro_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Phi_{pr},Trig Ele, All Pro Rec' , xtitle='#Phi_{pr} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    plot_acceptance( can, histos, 'phi_pro_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #Phi_{pr},Trig Ele, All Pro Rec' , xtitle='#Phi_{pr} (deg)', ytitle='counts',x_range=[-180,180])
    

    
    # kaonP
    
    plot_overlap_1d( can, histos, 'p_kp_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, P_{kp}, Trig Ele, All K^{+} Rec' , xtitle='P_{kp} (GeV)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN', 'All Rec'])
    
    plot_overlap_1d( can, histos, 'theta_kp_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab,     
                     title='Gen. vs Rec, #Theta_{kp} Trig Ele, All K^{+} Rec' , xtitle='#Theta_{kp} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])
    
    
    plot_acceptance( can, histos, 'theta_kp_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #Theta_{kp},Trig Ele, All K^{+} Rec' , xtitle='#Theta_{kp} (deg)', ytitle='Acceptance', x_range=[0,35])

    plot_overlap_1d( can, histos, 'phi_kp_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs All Rec, #Theta_{kp},Trig Ele, All K^{+} Rec' , xtitle='#Phi_{kp} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    plot_acceptance( can, histos, 'phi_kp_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #Phi_{kp}, Trig Ele, All K^{ +} Rec' , xtitle='#Phi_{kp} (deg)', ytitle='Acceptance')

    # kaonM
    plot_overlap_1d( can, histos, 'p_km_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs All Rec, P_{km},Trig Ele, All K^{ -} Rec' , xtitle='P_{km} (GeV)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    plot_overlap_1d( can, histos, 'theta_km_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Theta_{km},Trig Ele, All K^{ -} Rec' , xtitle='#Theta_{km} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])
    
    plot_acceptance( can, histos, 'theta_km_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #Theta_{km}, Trig Ele, All K^{ -} Rec' , xtitle='#Theta_{km} (deg)', ytitle='Acceptance', x_range=[0,35])

    plot_overlap_1d( can, histos, 'phi_km_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Theta_{km} Final, Trig Ele, All K^{ -} Rec' , xtitle='#Phi_{km} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','All REC'])

    plot_acceptance( can, histos, 'phi_km_{}', ['sim_gen','sim_rec_all'], output_pdfname, lab, 
                     title='Accep, #Phi_{km}, Trig Ele, All K^{ -} Rec' , xtitle='#Phi_{km} (deg)', ytitle='Acceptance')
    
    ##############################################################
    ##############################################################
    ## check the distributions for when the phi is presented - 
    
    # lets look at the electron theta and phi distributions and acceptance
    plot_overlap_1d( can, histos, 'p_ele_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, P_{e} Final #phi (FD) Events' , xtitle='P_{el} (GeV)', ytitle='counts', set_marker_style=True, logy=True,  name=['GEN','#phi REC'])
    
    plot_overlap_1d( can, histos, 'theta_ele_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, #Theta_{e} Final #phi (FD) Events' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True, logy=True, name=['GEN','#phi REC'])

    plot_acceptance( can, histos, 'theta_ele_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Accep, #Theta_{e} Final #phi (FD) Events' , xtitle='#Theta_{el} (deg)', ytitle='Acceptance')
    
    
    
    ###### now check the sim rec and sim gen per rec for q2, -t and xb 
    plot_overlap_1d( can, histos, 'q2_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, Q^{2}, Final #phi' , xtitle='Q^{2} (GeV^{2})', ytitle='counts',  set_marker_style=True, name=['GEN','#phi REC'])

    plot_overlap_1d( can, histos, 'xb_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, X_{b}, Final #phi' , xtitle='X_{b}', ytitle='counts', set_marker_style=True, name=['GEN','#phi REC'])
    
    plot_overlap_1d( can, histos, 't_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, -t, Final #phi' , xtitle='-t (GeV^{2})', ytitle='counts',  set_marker_style=True, name=['GEN','#phi REC'])

    plot_overlap_1d( can, histos, 'phitrento_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, #phi_{trento}, Final #phi' , xtitle='#phi_{trento} (deg)', ytitle='counts',  set_marker_style=True, name=['GEN','#phi REC'])

    ###### get the acceptance for q2, -t , xb, phi trento
    # q2
    plot_overlap_1d( can, histos, 'q2_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV^{2})', ytitle='counts', logy=True, set_marker_style=True,  name=['GEN','#phi REC'])

    plot_acceptance( can, histos, 'q2_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Accep, Q^{2} Final #phi (FD) Events' , xtitle='Q^{2} (GeV^{2})', ytitle='Acceptance', x_range=[0,7])

    # xb
    plot_overlap_1d( can, histos, 'xb_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, X_{b} Final #phi (FD) Events' , xtitle='X_{b}', ytitle='counts', logy=True, set_marker_style=True,  name=['GEN','#phi REC'])

    plot_acceptance( can, histos, 'xb_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Accp, X_{b} Final #phi (FD) Events' , xtitle='X_{b}', ytitle='counts', x_range=[0.0, 0.6])

    # -t 
    plot_overlap_1d( can, histos, 't_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, -t Final #phi (FD) Events' , xtitle='-t (GeV^{2})', ytitle='counts', logy=True, set_marker_style=True,  name=['GEN','#phi REC'])

    plot_acceptance( can, histos, 't_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Accep, -t Final #phi (FD) Events' , xtitle='-t (GeV^{2})', ytitle='Acceptance', x_range=[0,2])
    # phi trento
    plot_overlap_1d( can, histos, 'phitrento_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Gen. vs Rec, #phi_{trento} Final #phi (FD) Events' , xtitle='#phi_{trento} (deg)', ytitle='counts', logy=True, set_marker_style=True, name=['GEN','#phi REC'])

    plot_acceptance( can, histos, 'phitrento_{}', ['sim_gen','sim_rec'], output_pdfname, lab, 
                     title='Accep, #phi_{trento} Final #phi (FD) Events' , xtitle='#phi_{trento} (deg)', ytitle='Acceptance', x_range = [-180,180])


    '''


    ## compare generated and reconstructed for the same event. 

    
    can.Clear()
    can.Update()
    plot_overlap_1d( can, histos, 'p_ele_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, P_{e} Final #phi' , xtitle='p_{el} (GeV)', ytitle='counts', set_marker_style=True)

    plot_overlap_1d( can, histos, 'theta_ele_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, #Theta_{e} Final #phi' , xtitle='#Theta_{el} (deg)', ytitle='counts', set_marker_style=True)
    
    ## now for protons 
    plot_overlap_1d( can, histos, 'p_pro_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, P_{pr} Final #phi' , xtitle='p_{pr} (GeV)', ytitle='counts', set_marker_style=True)

    plot_overlap_1d( can, histos, 'theta_pro_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
                     title='FD:Rec and Gen/Event, #Theta_{pr} Final #phi' , xtitle='#Theta_{pr} (deg)', ytitle='counts', set_marker_style=True)
    

    ## now for KaonP
    #plot_overlap_1d( can, histos, 'p_kp_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
    #                 title='FD:Rec and Gen/Event, P_{K^{ +}} Final #phi' , xtitle='p_{K^{ +}} (GeV)', ytitle='counts', set_marker_style=True)

    #plot_overlap_1d( can, histos, 'theta_kp_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
    #                 title='FD:Rec and Gen/Event, #Theta_{K^{ +}} Final #phi' , xtitle='#Theta_{K^{ +}} (deg)', ytitle='counts',set_marker_style=True) 

    ## now for KaonM
    #plot_overlap_1d( can, histos, 'p_km_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
    #                 title='FD:Rec and Gen/Event, P_{K^{ -}} Final #phi' , xtitle='p_{K^{ -}} (GeV)', ytitle='counts',  set_marker_style=True)

    #plot_overlap_1d( can, histos, 'theta_km_{}', ['sim_rec','sim_gen_per_rec'], output_pdfname, lab, 
    #                 title='FD:Rec and Gen/Event, #Theta_{K^{ -}} Final #phi' , xtitle='#Theta_{K^{ -}} (deg)', ytitle='counts', set_marker_style=True) 
 

    '''
    
    
    
    
    

    
    

    can.Print('{}]'.format(output_pdfname))
        
