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
                  gPad, gStyle, TLatex, TGraphErrors, TMultiGraph, TLegend,
                  TLine, 
                  kRed, kBlack, kBlue, kWhite, kGreen, kMagenta, kYellow, kPink, kSpring)

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



def add_text_page(can, label, text, save_name):
    """ Write some text onto a page. """
    can.Clear()
    can.cd(1)
    label.DrawLatex(0.1, 0.5, text)
    can.Print(save_name)
    

def formatAsyTxt(fin, sig_range ):
    x_cnt = array('d')
    y_val = array('d')
    x_err = array('d')
    y_err = array('d')

    x_cnt_sig = array('d')
    y_val_sig = array('d')
    x_err_sig = array('d')
    y_err_sig = array('d')

    cc=0
    weighted_avg = 0
    weighted_avg_den = 0
    x_cnt_sig_avg=0
    print(' loading parameters')
    for ll in fin:
        rr = ll.split(' ')
        print rr
        a= float(rr[0])
        b=float(rr[1])
        c=float(rr[2])
        

        ## signal
        if cc > sig_range[0] and cc < sig_range[1]: # typically 0, 5
            print(' loading signal cc {} '.format(cc))
            x_cnt_sig.append(a)
            y_val_sig.append(b)
            x_err_sig.append(0.0)
            y_err_sig.append(c)
            ## calculate most prob weighted average
            ## assuming sampling from IID , gaus
            weighted_avg+= ( b*c)
            weighted_avg_den+= (c)
            x_cnt_sig_avg+=a 
            print('weighted avg {} '.format(weighted_avg/weighted_avg_den))
                            

        else:
            x_cnt.append(a)
            y_val.append(b)
            x_err.append(0.0)
            y_err.append(c)
        cc+=1

    return {'x_cnt':x_cnt,'y_val': y_val, 'x_err':x_err, 'y_err':y_err, 'x_cnt_sig':x_cnt_sig, 'y_val_sig':y_val_sig,'x_err_sig': x_err_sig, 
            'y_err_sig':y_err_sig,'weighted_avg': weighted_avg, 'weighted_avg_den':weighted_avg_den, 'x_cnt_sig_avg':x_cnt_sig_avg}
    
if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    output_pdfname = args.output_prefix + '.pdf'
    setup_global_options() 

    can = TCanvas('can', 'can', 1000, 600)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    can.Clear()
    can.Update()

    asy_per_mass = open('asy_vs_imkpkm_pidtype1.txt','r')
    asy_per_mass_remove_resonant1 = open('asy_vs_imkpkm_pidtype1_remove_resonant1.txt','r')
    asy_per_mass_remove_resonant2 = open('asy_vs_imkpkm_pidtype1_remove_resonant2.txt','r')

    asy_res_nom = formatAsyTxt(asy_per_mass, [0,5])
    asy_res_rr1 = formatAsyTxt(asy_per_mass_remove_resonant1, [0,5])
    asy_res_rr2 = formatAsyTxt(asy_per_mass_remove_resonant2, [0,5])


    # nominal
    x_cnt = asy_res_nom['x_cnt'] 
    y_val = asy_res_nom['y_val']
    x_err = asy_res_nom['x_err']
    y_err = asy_res_nom['y_err']
    x_cnt_sig = asy_res_nom['x_cnt_sig']
    y_val_sig = asy_res_nom['y_val_sig']
    x_err_sig = asy_res_nom['x_err_sig']
    y_err_sig = asy_res_nom['y_err_sig']
    weighted_avg = asy_res_nom['weighted_avg']
    weighted_avg_den = asy_res_nom['weighted_avg_den']
    x_cnt_sig_avg = asy_res_nom['x_cnt_sig_avg']


    # remove first 1520 resonant
    x_cnt_rr1 = asy_res_rr1['x_cnt'] 
    y_val_rr1 = asy_res_rr1['y_val']
    x_err_rr1 = asy_res_rr1['x_err']
    y_err_rr1 = asy_res_rr1['y_err']
    x_cnt_sig_rr1 = asy_res_rr1['x_cnt_sig']
    y_val_sig_rr1 = asy_res_rr1['y_val_sig']
    x_err_sig_rr1 = asy_res_rr1['x_err_sig']
    y_err_sig_rr1 = asy_res_rr1['y_err_sig']
    weighted_avg_rr1 = asy_res_rr1['weighted_avg']
    weighted_avg_den_rr1 = asy_res_rr1['weighted_avg_den']
    x_cnt_sig_avg_rr1 = asy_res_rr1['x_cnt_sig_avg']

    # remove second reson
    x_cnt_rr2 = asy_res_rr2['x_cnt'] 
    y_val_rr2 = asy_res_rr2['y_val']
    x_err_rr2 = asy_res_rr2['x_err']
    y_err_rr2 = asy_res_rr2['y_err']
    x_cnt_sig_rr2 = asy_res_rr2['x_cnt_sig']
    y_val_sig_rr2 = asy_res_rr2['y_val_sig']
    x_err_sig_rr2 = asy_res_rr2['x_err_sig']
    y_err_sig_rr2 = asy_res_rr2['y_err_sig']
    weighted_avg_rr2 = asy_res_rr2['weighted_avg']
    weighted_avg_den_rr2 = asy_res_rr2['weighted_avg_den']
    x_cnt_sig_avg_rr2 = asy_res_rr2['x_cnt_sig_avg']

    x_cnt_weight =  array('d')
    y_weight = array('d')
    x_cnt_weight_er = array('d')
    y_weight_err = array('d')

    x_cnt_weight_rr1 =  array('d')
    y_weight_rr1 = array('d')
    x_cnt_weight_er_rr1 = array('d')
    y_weight_err_rr1 = array('d')

    x_cnt_weight_rr2 =  array('d')
    y_weight_rr2 = array('d')
    x_cnt_weight_er_rr2 = array('d')
    y_weight_err_rr2 = array('d')

    # get weighted variance
    weight_err = 0
    weight_err_den = 0
    print('--> weighted average is {}'.format(weighted_avg/weighted_avg_den))
    print('--> weighted bin center is {}'.format(x_cnt_sig_avg/len(x_cnt_sig)))
          
    for mm in range(0,len(y_val_sig)):
        weight_err+=(y_err_sig[mm] * (y_val_sig[mm] - (weighted_avg/weighted_avg_den))**2)
        

    weight_err= weight_err/weighted_avg_den * (len(y_err_sig)/(len(y_err_sig)-1))
    #print('--> weighted error {}'.format(weight_err))
    std_weight_err = np.sqrt(weight_err * (1.0/len(y_val_sig)))
    x_cnt_weight.append(x_cnt_sig_avg/len(x_cnt_sig))
    y_weight.append(weighted_avg/weighted_avg_den)
    x_cnt_weight_er.append(0.0)
    y_weight_err.append(weighted_avg_den/(len(y_err_sig)))#std_weight_err)
    print('--> weighted error {}'.format((weighted_avg_den/(len(y_err_sig))))) # average of the error bar values


    ## get the weighted info for removing first resonace 
    weight_err_rr1 = 0
    weight_err_den_rr1 = 0
    print('-->_rr1 weighted average is {}'.format(weighted_avg_rr1/weighted_avg_den_rr1))
    print('-->_rr1 weighted bin center is {}'.format(x_cnt_sig_avg_rr1/len(x_cnt_sig_rr1)))
          
    for mm in range(0,len(y_val_sig_rr1)):
        weight_err_rr1+=(y_err_sig_rr1[mm] * (y_val_sig_rr1[mm] - (weighted_avg_rr1/weighted_avg_den_rr1))**2)
        

    weight_err_rr1= weight_err_rr1/weighted_avg_den_rr1 * (len(y_err_sig_rr1)/(len(y_err_sig_rr1)-1))
    #print('-->_rr1 weighted error {}'.format(weight_err_rr1))
    std_weight_err = np.sqrt(weight_err_rr1 * (1.0/len(y_val_sig_rr1)))
    x_cnt_weight_rr1.append(x_cnt_sig_avg_rr1/len(x_cnt_sig_rr1))
    y_weight_rr1.append(weighted_avg_rr1/weighted_avg_den_rr1)
    x_cnt_weight_er_rr1.append(0.0)
    y_weight_err_rr1.append(weighted_avg_den_rr1/(len(y_err_sig_rr1)))#std_weight_err)
    print('-->_rr1 weighted error {}'.format((weighted_avg_den_rr1/(len(y_err_sig_rr1)))))


    ## get the weighted info for removing second resonace 
    weight_err_rr2 = 0
    weight_err_den_rr2 = 0
    print('-->_rr2 weighted average is {}'.format(weighted_avg_rr2/weighted_avg_den_rr2))
    print('-->_rr2 weighted bin center is {}'.format(x_cnt_sig_avg_rr2/len(x_cnt_sig_rr2)))
          
    for mm in range(0,len(y_val_sig_rr2)):
        weight_err_rr2+=(y_err_sig_rr2[mm] * (y_val_sig_rr2[mm] - (weighted_avg_rr2/weighted_avg_den_rr2))**2)
        

    weight_err_rr2 = weight_err_rr2/weighted_avg_den_rr2 * (len(y_err_sig_rr2)/(len(y_err_sig_rr2)-1))
    #print('-->_rr2 weighted error {}'.format(weight_err_rr2)) error on data point itself
    std_weight_err = np.sqrt(weight_err_rr2 * (1.0/len(y_val_sig_rr2)))
    x_cnt_weight_rr2.append(x_cnt_sig_avg_rr2/len(x_cnt_sig_rr2))
    y_weight_rr2.append(weighted_avg_rr2/weighted_avg_den_rr2)
    x_cnt_weight_er_rr2.append(0.0)
    y_weight_err_rr2.append(weighted_avg_den_rr2/(len(y_err_sig_rr2)))#std_weight_err)
    print('-->_rr2 weighted error {}'.format((weighted_avg_den_rr2/(len(y_err_sig_rr2)))))

    # plot the nominal weighted bsa without removing any resonant
    avg_asy_graph = TGraphErrors(len(x_cnt_weight), x_cnt_weight, y_weight, x_cnt_weight_er , y_weight_err)
    avg_asy_graph.SetTitle('')
    avg_asy_graph.SetMarkerStyle(29)
    avg_asy_graph.SetMarkerSize(1)
    avg_asy_graph.SetLineColor(kGreen)
    avg_asy_graph.SetMarkerColor(kGreen)

    ## remove the first 1520 resonance
    avg_asy_graph_rr1 = TGraphErrors(len(x_cnt_weight_rr1), x_cnt_weight_rr1, y_weight_rr1, x_cnt_weight_er_rr1 , y_weight_err_rr1)
    avg_asy_graph_rr1.SetTitle('')
    avg_asy_graph_rr1.SetMarkerStyle(29)
    avg_asy_graph_rr1.SetMarkerSize(0.75)
    avg_asy_graph_rr1.SetLineColor(kGreen+4)
    avg_asy_graph_rr1.SetMarkerColor(kGreen+4)

    ## remove the second lambda resonance
    avg_asy_graph_rr2 = TGraphErrors(len(x_cnt_weight_rr2), x_cnt_weight_rr2, y_weight_rr2, x_cnt_weight_er_rr2 , y_weight_err_rr2)
    avg_asy_graph_rr2.SetTitle('')
    avg_asy_graph_rr2.SetMarkerStyle(29)
    avg_asy_graph_rr2.SetMarkerSize(0.5)
    avg_asy_graph_rr2.SetLineColor(kSpring)
    avg_asy_graph_rr2.SetMarkerColor(kSpring)
    
    
    asy_per_mass_graph = TGraphErrors(len(x_cnt), x_cnt, y_val, x_err, y_err)
    asy_per_mass_graph.SetTitle('')
    asy_per_mass_graph.SetMarkerStyle(21)
    asy_per_mass_graph.SetMarkerSize(1.5)
    asy_per_mass_graph.SetLineColor(kBlack)
    asy_per_mass_graph.SetMarkerColor(kBlack)

    # first resonance
    asy_per_mass_graph_rr1 = TGraphErrors(len(x_cnt_rr1), x_cnt_rr1, y_val_rr1, x_err_rr1, y_err_rr1)
    asy_per_mass_graph_rr1.SetTitle('')
    asy_per_mass_graph_rr1.SetMarkerStyle(21)
    asy_per_mass_graph_rr1.SetMarkerSize(0.7)
    asy_per_mass_graph_rr1.SetLineColor(kMagenta+4)
    asy_per_mass_graph_rr1.SetMarkerColor(kMagenta+4)

    # second resonance
    asy_per_mass_graph_rr2 = TGraphErrors(len(x_cnt_rr2), x_cnt_rr2, y_val_rr2, x_err_rr2, y_err_rr2)
    asy_per_mass_graph_rr2.SetTitle('')
    asy_per_mass_graph_rr2.SetMarkerStyle(21)
    asy_per_mass_graph_rr2.SetMarkerSize(0.4)
    asy_per_mass_graph_rr2.SetLineColor(kMagenta-3)
    asy_per_mass_graph_rr2.SetMarkerColor(kMagenta-3)


    # nominal , no removal 
    sig_asy_per_mass_graph = TGraphErrors(len(x_cnt_sig), x_cnt_sig, y_val_sig, x_err_sig, y_err_sig)
    sig_asy_per_mass_graph.SetTitle('')
    sig_asy_per_mass_graph.SetMarkerStyle(20)
    sig_asy_per_mass_graph.SetMarkerSize(2)    
    sig_asy_per_mass_graph.SetLineColor(kRed)
    sig_asy_per_mass_graph.SetMarkerColor(kRed)

    # remove first resonsance
    sig_asy_per_mass_graph_rr1 = TGraphErrors(len(x_cnt_sig_rr1), x_cnt_sig_rr1, y_val_sig_rr1, x_err_sig_rr1, y_err_sig_rr1)
    sig_asy_per_mass_graph_rr1.SetTitle('')
    sig_asy_per_mass_graph_rr1.SetMarkerStyle(20)
    sig_asy_per_mass_graph_rr1.SetMarkerSize(2)    
    sig_asy_per_mass_graph_rr1.SetLineColorAlpha(kRed+3, 0.75)
    sig_asy_per_mass_graph_rr1.SetMarkerColorAlpha(kRed+3, 0.75)

    sig_asy_per_mass_graph_rr2 = TGraphErrors(len(x_cnt_sig_rr2), x_cnt_sig_rr2, y_val_sig_rr2, x_err_sig_rr2, y_err_sig_rr2)
    sig_asy_per_mass_graph_rr2.SetTitle('')
    sig_asy_per_mass_graph_rr2.SetMarkerStyle(20)
    sig_asy_per_mass_graph_rr2.SetMarkerSize(2)    
    sig_asy_per_mass_graph_rr2.SetLineColorAlpha(kRed-6, 0.50)
    sig_asy_per_mass_graph_rr2.SetMarkerColorAlpha(kRed-6, 0.50)

    mg_asy_per_mass = TMultiGraph()
    mg_asy_per_mass.SetTitle('')
    mg_asy_per_mass.Add(asy_per_mass_graph)
    mg_asy_per_mass.Add(sig_asy_per_mass_graph)
    mg_asy_per_mass.Add(avg_asy_graph)
    mg_asy_per_mass.Add(asy_per_mass_graph_rr1)
    mg_asy_per_mass.Add(sig_asy_per_mass_graph_rr1)
    mg_asy_per_mass.Add(avg_asy_graph_rr1)

    mg_asy_per_mass.Add(asy_per_mass_graph_rr2)
    mg_asy_per_mass.Add(sig_asy_per_mass_graph_rr2)
    mg_asy_per_mass.Add(avg_asy_graph_rr2)

    mg_asy_per_mass.GetHistogram().SetMaximum(0.5)
    mg_asy_per_mass.GetHistogram().SetMinimum(-0.5)
    mg_asy_per_mass.GetXaxis().SetLimits(0.95, 1.4)
    
    # draw grid
    h = gPad.DrawFrame(1.4, -0.5, 0.95, 0.5)
    h.GetXaxis().SetNdivisions(-9)
    h.GetYaxis().SetNdivisions(-9)
    #h.GetXaxis().SetLabelSize(0.025)
    #h.GetYaxis().SetLabelSize(0.025)
    gPad.SetGrid()
    

    mg_asy_per_mass.Draw('AP')

    
    background_fit = TF1('bck_asy',"[0]",1.031, 1.23)
    asy_per_mass_graph.Fit("bck_asy","R")
    background_fit.SetLineColor(kBlue)
    background_fit.SetLineWidth(3)
    background_fit.Draw('same')

    ## fit background from data without resonant 1, but with resonant 2
    background_fit_rr1 = TF1('bck_asy_rr1',"[0]",1.031, 1.23)
    asy_per_mass_graph_rr1.Fit("bck_asy_rr1","R")
    background_fit_rr1.SetLineWidth(2)
    background_fit_rr1.SetLineStyle(4)    
    background_fit_rr1.SetLineColorAlpha(kYellow, 0.75)
    background_fit_rr1.Draw('same')

    ## fit background from data without resonant 2, but with resonant 1
    background_fit_rr2 = TF1('bck_asy_rr2',"[0]",1.031, 1.23)
    asy_per_mass_graph_rr2.Fit("bck_asy_rr2","R")
    background_fit_rr2.SetLineWidth(1)
    background_fit_rr2.SetLineStyle(4)    
    background_fit_rr2.SetLineColorAlpha(kPink, 0.35 )
    background_fit_rr2.Draw('same')
    
    start_range =1.03
    end_range = 1.23
    delta_shift=0.05
    for bg in range(1,5):
        start_range+=delta_shift
        end_range+=delta_shift
        print('start background fit range {} - end range {} '.format(start_range,end_range))

    lab.DrawLatex(0.1, 0.925, 'BSA Per I.M. of K^{ +}K^{ -} Slice')        
    lab.DrawLatex(0.5, 0.015,'I.M. of K^{ +}K^{ -} (GeV)')            
    lab.SetTextAngle(90)
    lab.DrawLatex(0.035, 0.5, 'Beam Spin Asy.')
    lab.SetTextAngle(0)

    leg = TLegend(0.5, 0.6, 0.9, 0.9)
    leg.AddEntry(asy_per_mass_graph,'Sideband Asy','p')
    leg.AddEntry(sig_asy_per_mass_graph,'Sig. Asy','p')
    leg.AddEntry(avg_asy_graph,'Weighted Asy','p')
    leg.AddEntry(asy_per_mass_graph_rr1,'Sideband Asy, Remove 1520','p')
    leg.AddEntry(sig_asy_per_mass_graph_rr1,'Sig. Asy, Remove 1520','p')
    leg.AddEntry(avg_asy_graph_rr1,'Weighted Asy, Remove 1520','p')
    leg.AddEntry(asy_per_mass_graph_rr2,'Sideband Asy, Remove 1820','p')
    leg.AddEntry(sig_asy_per_mass_graph_rr2,'Sig. Asy, Remove 1820','p')
    leg.AddEntry(avg_asy_graph_rr2,'Weighted Asy, Remove 1820','p')
    #leg.SetFillStyle(0)
    leg.Draw('same')

    can.Print(output_pdfname)


    can.Print('{}]'.format(output_pdfname))
