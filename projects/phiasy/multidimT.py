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
#  
#if you have a canvas, TCanvas *canvas, do
#  canvas->SetBatch(kTRUE);


from array import array 
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors, TLegend,
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
    
    # original binning
    #q2_bins = [1, 2, 2.5, 3, 3.5]
    #xb_bins = [0.1, 0.20, 0.275, 0.35, 0.7]
    q2_bins = [1, 3, 8]
    xb_bins = [0.1, 0.4, 0.8]

    can.Divide(len(xb_bins)-1, len(q2_bins)-1)#, 0, 0)
    
    cc=1
    for qq in range(1, len(q2_bins)):
        for xx in range(1, len(xb_bins)):
            can.cd(cc)
            t_in_q2_xb = histos['t_q2b{}_xbb{}'.format(qq,xx)]
            t_in_q2_xb.Draw()
            lab.DrawLatex(0.1, 0.925, 'Final #phi, -t, 1.010<Mass(K^{ +}K^{ -})<1.050')
            lab.DrawLatex(0.5, 0.8, 'Q2B{}:{}-{}'.format(qq, q2_bins[qq-1], q2_bins[qq]))
            lab.DrawLatex(0.5, 0.75, 'XbB{}:{}-{}'.format(xx, xb_bins[xx-1], xb_bins[xx]))
            lab.DrawLatex(0.6, 0.025, '-t (GeV^{2})')
            
            cc+=1

            
    can.Print(output_pdfname)
    can.Print('{}]'.format(output_pdfname))

