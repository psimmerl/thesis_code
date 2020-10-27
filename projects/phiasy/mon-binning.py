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
                  TLine, kRed, kBlack, kBlue, kWhite, kTRUE, gROOT)

from ROOT import( RooRealVar, RooArgSet, RooArgList, RooDataHist,
                  RooGaussian, RooFit, RooAddPdf, RooPolynomial, RooFormulaVar, RooLinkedList)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
global_dump = []

gROOT.SetBatch(kTRUE)

def readInBinningInfo(fname):
    fin = open(fname,'r')
    temp = []
    for ll in fin:
        print ll
        temp.append(float(ll))
    return temp    

def readInOverlapBinningInfo(fname):
    ovrlp_list = []
    fin = open(fname,'r')
    for ll in fin:
        temp_ll = ll.split(' ')
        print(temp_ll)
        ovrlp_list.append([float(temp_ll[0]), float(temp_ll[1])])
    return ovrlp_list


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

def plot_slices( canvas, histos, title_formatter, save_name, label,
                      title=None, xtitle=None, ytitle=None , plot_details=None, field=None):

    canvas.Clear()
    canvas.Divide(1,1)

    root_is_dumb = []
    canvas.cd(1)
    page_max = 0
    signal_sum=0
    
    ## define max of y axis here
    ymax = histos.get(title_formatter.format(2), default_histo).GetMaximum()
    
    start_bin = 1
    if 'sliding_bin_list' in plot_details.keys():
        start_bin=0


    for bb in range(start_bin, plot_details['max_bins']+1):
        if isinstance(histos.get(title_formatter.format(bb), default_histo), TH1F):
            h_temp = histos.get(title_formatter.format(bb), default_histo)
            print('--> Number of Entries: {}'.format(h_temp.GetEntries()))
            #if bb > 0 and bb < 5: signal_sum+=h_temp.GetEntries()
            if 'colors' in plot_details.keys():
                if bb%2 == 0:
                    print(' alternating color of bin groups - RED' )
                    h_temp.SetLineColor(kRed)
                else:
                    print(' alternating color of bin groups - BLUE' )
                    h_temp.SetLineColor(kBlue)

            
            h_temp.Draw('same')
            root_is_dumb.append(h_temp)
            if( bb == 1 and not 'sliding_bin_list' in plot_details.keys() and 'max_xaxis' in plot_details.keys()):
                h_temp.SetMaximum(ymax*1.2)
                h_temp.GetXaxis().SetRangeUser(0.95, plot_details['max_xaxis'])#1.0861)#134)
                #label.DrawLatex(0.55,0.7,'Total Entries In Peak Range'+' = {}'.format(h_temp.GetEntries()))
                canvas.Update()
            if( bb == 0 and 'sliding_bin_list' in plot_details.keys()):
                h_temp.SetMaximum(1575)
                h_temp.GetXaxis().SetRangeUser(0.95, plot_details['max_xaxis'])#1.0861)#134)
                #label.DrawLatex(0.55,0.7,'Total Entries In Peak Range'+' = {}'.format(h_temp.GetEntries()))
                canvas.Update()
    print('--> Total number of entries in signal region is {}'.format(signal_sum))

    if title:
        label.DrawLatex(0.1, 0.925, title)

    if xtitle:
        label.DrawLatex(0.4, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.035, 0.5, ytitle)
        label.SetTextAngle(0)
    
    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)
        

    canvas.Print(save_name)


def plot_page_bins( canvas, histos, title_formatter, save_name, label,
                      title=None, xtitle=None, ytitle=None , plot_details=None, field=None, vlines=None):


    canvas.Clear()

    root_is_dumb = []    
    pol=1
    if 'bp' in plot_details.keys():
        pol= plot_details['bp']

    g_list = []
    ncols=3
    nrows=int(math.ceil( (plot_details['max_bins']-1)/ncols) ) + 1
    print(' number of cols {} number of rows {}'.format(ncols,nrows))
    canvas.Divide(2,2)#ncols, nrows)

    start_bin = 1
    can_shift=0
    if 'sliding_bin_list' in plot_details.keys():
        start_bin=0
        can_shift=1

    for bb in range(start_bin, plot_details['max_bins']):
        canvas.cd(bb+can_shift)
        if isinstance(histos.get(title_formatter.format('pos',bb), default_histo), TH1F):
            print(' --> mass bin {}'.format(bb) )
            h_temp_pos = histos.get(title_formatter.format('pos',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('pos',bb)).GetTitle())
            h_temp_neg = histos.get(title_formatter.format('neg',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('neg',bb)).GetTitle())

            h_temp_pos.SetLineColor(kRed)
            h_temp_neg.SetLineColor(kBlue)

            leg=TLegend(0.5, 0.7, 0.7, 0.9)
            leg.AddEntry(h_temp_pos,'pos','l')
            leg.AddEntry(h_temp_neg,'neg','l')
            
            h_temp_pos.SetMinimum(0)
            h_temp_pos.SetMaximum(h_temp_pos.GetMaximum()*1.1)

            h_temp_pos.Draw('HE')
            h_temp_neg.Draw('HE+same')
            leg.Draw('same')

            root_is_dumb.append(h_temp_pos)
            root_is_dumb.append(h_temp_neg)
            root_is_dumb.append(leg)

            
            if 'bin_list' in plot_details.keys():
                label.DrawLatex(0.22, 0.925, '{0:.4f}<{2}<{1:.4f}'.format(plot_details['bin_list'][bb-1], plot_details['bin_list'][bb], 'K^{ +}K^{ -}'))            

            if 'sliding_bin_list' in plot_details.keys():
                print(plot_details['sliding_bin_list'][bb])
                label.DrawLatex(0.22, 0.925, '{0:.4f}<{2}<{1:.4f}'.format(plot_details['sliding_bin_list'][bb][0], plot_details['sliding_bin_list'][bb][1], 'K^{ +}K^{ -}'))            

            if( h_temp_pos.GetMaximum() >  h_temp_neg.GetMaximum() ): 
                h_temp_pos.SetMaximum( h_temp_pos.GetMaximum()*1.1)
            else:                    
                h_temp_neg.SetMaximum( h_temp_neg.GetMaximum()*1.1 )
                
            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.5, 0.025, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.035, 0.5, ytitle)
                label.SetTextAngle(0)

            if field:
                flabel=TLatex()
                flabel.SetTextColorAlpha(1, 0.376)
                flabel.SetNDC()
                flabel.SetTextFont(42)
                flabel.SetTextSize(0.08)
                flabel.DrawLatex(0.125, 0.847, field)


    canvas.Print(save_name)            
    
            

def plot_asy( canvas, histos, title_formatter, save_name, label,
                      title=None, xtitle=None, ytitle=None , plot_details=None, f_out=None, f_out2=None, field=None):

    canvas.Clear()
    canvas.Update()

    root_is_dumb = []    
    pol=1
    if 'bp' in plot_details.keys():
        pol= plot_details['bp']

    g_list = []
    ncols=3
    nrows=int(math.ceil(plot_details['max_bins']-1)/ncols) + 1
    print(' number of cols {} number of rows {}'.format(ncols,nrows))
    canvas.Divide(2,2)#ncols, nrows)

    fout=None
    if f_out:
        fout = open(f_out,'w')

    start_bin = 1
    can_shift=0
    if 'sliding_bin_list' in plot_details.keys():
        start_bin=0
        can_shift=1
    asy_per_mass = []
    for bb in range(start_bin, plot_details['max_bins']):
        #if bb != 2: continue ## CHANGE BACK # SO REMOVE ME
        canvas.cd(bb+can_shift)
        if isinstance(histos.get(title_formatter.format('pos',bb), default_histo), TH1F) and isinstance(histos.get(title_formatter.format('neg',bb), default_histo), TH1F) :
            h_temp_pos = histos.get(title_formatter.format('pos',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('pos',bb)).GetTitle())
            h_temp_neg = histos.get(title_formatter.format('neg',bb), default_histo).Clone("clone"+histos.get(title_formatter.format('neg',bb)).GetTitle())

            bc=0
            if 'bin_list' in plot_details.keys():
                bc=plot_details['bin_list'][bb-1] + (plot_details['bin_list'][bb] - plot_details['bin_list'][bb-1])/2.0
            elif 'sliding_bin_list' in plot_details.keys():
                bc=plot_details['sliding_bin_list'][bb][0] + (plot_details['sliding_bin_list'][bb][1] - plot_details['sliding_bin_list'][bb][0])/2.0
            x_center = array('d')
            y_val = array('d')
            x_err = array('d')
            y_err = array('d')
            
             
            fout_phi_dist = None
            if f_out2:
                fout_phi_dist = open('asy_vs_trento_method1_imkpkm_bin{}_'.format(bb)+f_out2,'w')
            for ab in range(1,h_temp_pos.GetXaxis().GetNbins()+1):
                P = h_temp_pos.GetBinContent(ab) 
                N = h_temp_neg.GetBinContent(ab)
                x_cent = h_temp_pos.GetBinCenter(ab)
                nom = P-N
                den = P+N
                y=0
                yerr = 0
                print(" number of positive values %d, number of negative values %d" % (P, N))        
                if(den==0):
                    y=0
                    yerr=0
                else:
                    y=nom/den
                    yerr = math.sqrt( (1 - y**2)/(P+N))

                print(" bin center {}, asy {}, bin center error {}, asy error {} ".format(x_cent, y, 0, yerr))
                x_center.append(x_cent)
                y_val.append(y*(1.0/pol))
                x_err.append(0)
                y_err.append(yerr)
                fout_phi_dist.write(str(x_cent) + ' ' + str(y*(1.0/pol)) + ' ' + str(yerr) +' \n')
            if f_out2:
                fout_phi_dist.close()
            ##############
            graph = TGraphErrors(len(x_center), x_center, y_val, x_err, y_err)
            graph.SetMarkerStyle(22)
            graph.SetMarkerSize(2)
            graph.GetHistogram().SetMaximum(0.65)
            graph.GetHistogram().SetMinimum(-0.65)
            graph.GetXaxis().SetRangeUser(-180.0, 180.0)
            if bb%2 == 0: graph.SetMarkerColor(kRed)
            else: graph.SetMarkerColor(kBlue)
            

            # draw grid
            h = gPad.DrawFrame(-180, -0.65, -180, 0.65)
            h.GetXaxis().SetNdivisions(-7)
            h.GetYaxis().SetNdivisions(-7)
            h.GetXaxis().SetLabelSize(0.025)
            h.GetYaxis().SetLabelSize(0.025)
            gPad.SetGrid()
            
            
            g_list.append(graph)
            graph.Draw('AP')
            gPad.RedrawAxis("g")

            fit_asy = TF1('fasy{}'.format(bb),"[0]*sin(x*TMath::Pi()/180.0)", -180, 180)#, 361.0)
            fit_asy.SetParameter(0,1)
            graph.Fit('fasy{}'.format(bb))
            fit_asy.SetLineColor(kRed)
            fit_asy_chi2 = fit_asy.GetChisquare()
            fit_asy_ndf = fit_asy.GetNDF()
            root_is_dumb.append(fit_asy)
            if f_out:
                fout.write(str(bc) + ' ' + str(fit_asy.GetParameter(0)) + ' ' + str(fit_asy.GetParError(0)) + ' ' + str(fit_asy_chi2) + ' ' + str(fit_asy_ndf) + ' \n')
            
            
            # draw bin related info
            if 'bin_list' in plot_details.keys():
                label.DrawLatex(0.35, 0.925, '             ,{0:.4f}<{2}<{1:.4f}'.format(plot_details['bin_list'][bb-1], plot_details['bin_list'][bb],'K^{ +}K^{ -}'))
            elif 'sliding_bin_list' in plot_details.keys():
                label.DrawLatex(0.35, 0.925, '             ,{0:.4f}<{2}<{1:.4f}'.format(plot_details['sliding_bin_list'][bb][0], plot_details['sliding_bin_list'][bb][1],'K^{ +}K^{ -}'))
            label.DrawLatex(0.55, 0.835, 'BC:  {:.3f} GeV'.format(bc) )
            label.DrawLatex(0.55, 0.785, 'Asy_{m}'+': {:.3f}'.format(fit_asy.GetParameter(0)))
            label.DrawLatex(0.55, 0.728, '#Delta Asy_{m}'+': {:.3f}'.format(fit_asy.GetParError(0)))
            ## was 0.55 beefore 
            label.DrawLatex(0.55, 0.675, 'Chi2/NDF:{0:.1f}/{1}'.format(fit_asy_chi2,fit_asy_ndf))
                        
            
            if title:
                label.DrawLatex(0.1, 0.925, title)

            if xtitle:
                label.DrawLatex(0.5, 0.025, xtitle)

            if ytitle:
                label.SetTextAngle(90)
                label.DrawLatex(0.045, 0.5, ytitle) # was 0.035
                label.SetTextAngle(0)
        
            if field:
                flabel=TLatex()
                flabel.SetTextColorAlpha(1, 0.376)
                flabel.SetNDC()
                flabel.SetTextFont(42)
                flabel.SetTextSize(0.08)
                flabel.DrawLatex(0.125, 0.847, field)


            asy_per_mass.append([bc, fit_asy.GetParameter(0), fit_asy.GetParError(0)])
            

    

    if f_out:
        fout.close()
    canvas.Print(save_name)


def plot_page(canvas, histos, title_formatter, label, save_name,
                        xtitle=None, ytitle=None, title=None, log=False,
                        y_fit_range=None, landscape=False, x_range=None, vline=None, vline2=None, hlines = None, shades = None,
                        hline=None, vlines=None, field=None):

    root_garbage_can = []
    root_is_dumb = []
    canvas.Clear() 
    if landscape:
        canvas.Divide(1,1)
    else:
        canvas.Divide(1,1)
    canvas.cd(1)

    if isinstance(histos.get(title_formatter, default_histo), TH1F):
        
        if x_range:
            histos.get(title_formatter, default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
        if log:
            histos.get(title_formatter, default_histo).SetMinimum(0.1)
        histos.get(title_formatter, default_histo).Draw('same')
            
    if isinstance(histos.get(title_formatter, default_histo), TH2F):
        histos.get(title_formatter, default_histo).Draw('colz')
        
    if title:
        label.DrawLatex(0.1, 0.925, title)
        
    if xtitle:
        label.DrawLatex(0.5, 0.025, xtitle)
        
    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    if hlines:
        for hl in hlines:
            xmin = histos.get(title_formatter).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter).GetXaxis().GetXmax() 
            line = TLine(xmin, hl, xmax, hl)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
        
    if vline:
        line = TLine(vline, histos.get(title_formatter, default_histo).GetYaxis().GetXmin(),
                     vline, histos.get(title_formatter, default_histo).GetYaxis().GetXmax())
        line.SetLineColor(1)
        line.SetLineStyle(1)
        line.Draw('same')
        root_garbage_can.append(line)

    if vlines:
        for vv in vlines:
            gPad.Update()
            ymax=gPad.GetUymax()                
            line = TLine(vv, 0,
                         vv, ymax)
    
            line.SetLineColor(kRed)
            line.SetLineStyle(1)
            line.Draw('same')
            global_dump.append(line)

    if shades:
        ii=0
        while ii < len(shades):
            gPad.Update()
            ymax=gPad.GetUymax()                
            box = TBox(shades[ii], 0.0, shades[ii+1], ymax)
            box.SetFillColorAlpha(4+ii,0.5)
            ii+=2
            box.Draw("same")
            global_dump.append(box)


    if field:
        flabel=TLatex()
        flabel.SetTextColorAlpha(1, 0.376)
        flabel.SetNDC()
        flabel.SetTextFont(42)
        flabel.SetTextSize(0.08)
        flabel.DrawLatex(0.125, 0.847, field)


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

    can = TCanvas('can', 'can', 1800, 1600)#3000, 3000)
    can.SetBatch(kTRUE)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.047)
    lab.SetTextColor(1)


    field_setting = 'inbNoutb'
    bp = 0
    max_xaxis = 0
    binvar='binVar11'
    print('--> FIELD {} , variation {}, BP {}'.format(field_setting,binvar,bp))

    can.Print('{}['.format(output_pdfname))

    mass_bins = readInBinningInfo('bin_mass_rangesVar11_'+field_setting+'.txt')
    overlap_mass_bins = readInOverlapBinningInfo('bin_mass_slidingVar4_'+field_setting+'.txt')
    print(overlap_mass_bins)

    q2_bins=readInBinningInfo('bin_limits_q2_'+field_setting+'.txt')
    xb_bins=readInBinningInfo('bin_limits_xb_'+field_setting+'.txt')
    t_bins=readInBinningInfo('bin_limits_t_'+field_setting+'.txt')
    w_bins=readInBinningInfo('bin_limits_w_'+field_setting+'.txt')
    
    q2_centers = readInBinningInfo('bin_centers_q2_'+field_setting+'.txt')
    xb_centers = readInBinningInfo('bin_centers_xb_'+field_setting+'.txt')
    t_centers = readInBinningInfo('bin_centers_t_'+field_setting+'.txt')
    w_centers = readInBinningInfo('bin_centers_w_'+field_setting+'.txt')
    
    mass_min = mass_bins[1]
    mass_max = mass_bins[2]

    #########################
    ## official plots for the thesis for binned asymmetry. These are more accurate 
    ## to the procedure. the ones in mon-phi event are different because of the cut on the 
    ## kpkm mass depends on field setting. here it is just one value

    plot_page(can, histos, 'q2_final', lab, output_pdfname,
              xtitle='Q^{2} (GeV^{2})', ytitle='Counts', title='(FD), Q^{2}, Final Sample, '+str(mass_min)+'<K^{ +}K^{ -}<'+str(mass_max),
              field=field_setting, vlines = q2_bins)

    plot_page(can, histos, 't_final', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='Counts', title='(FD), -t, Final Sample, '+str(mass_min)+'<K^{ +}K^{ -}<'+str(mass_max),
              field=field_setting, vlines = t_bins)

    plot_page(can, histos, 'xb_final', lab, output_pdfname,
              xtitle='X_{b}', ytitle='Counts', title='(FD), X_{b}, Final Sample, '+str(mass_min)+'<K^{ +}K^{ -}<'+str(mass_max),
              field=field_setting, vlines = xb_bins)

    plot_page(can, histos, 'w_final', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='Counts', title='(FD), W, Final Sample, '+str(mass_min)+'<K^{ +}K^{ -}<'+str(mass_max),
              field=field_setting, vlines = w_bins)

    ##############3

    plot_page(can, histos, 'xbt_lowQ2_final', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='x_{b}', title='-t vs x_{b} - Final Events: Q^{2}<2.5')

    plot_page(can, histos, 'xbt_highQ2_final', lab, output_pdfname,
              xtitle='-t (GeV^{2})', ytitle='x_{b}', title='-t vs x_{b} - Final Events: Q^{2}>2.5')
              #x_range=None, vline=None, vline2=None, hlines = None, shades = None,

    plot_page(can, histos, 'wt_final', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='-t (GeV^{2})', title='W vs -t - Final Events')
    
    plot_page(can, histos, 'wt_lowQ2_final', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='-t (GeV^{2})', title='W vs -t - Final Events: Q^{2}<2.5')

    plot_page(can, histos, 'wt_highQ2_final', lab, output_pdfname,
              xtitle='W (GeV)', ytitle='-t (GeV^{2})', title='W vs -t - Final Events: Q^{2}>2.5')
    
    plot_page(can, histos, 'mass_small_bin_wbin1', lab, output_pdfname, 
              xtitle='K^{ +}K^{ -} (GeV)', ytitle='counts', title='W < 2.1')       
    plot_page(can, histos, 'mass_small_bin_wbin2', lab, output_pdfname,
              xtitle='K^{ +}K^{ -} (GeV)', ytitle='counts', title='W > 2.1 and W < 2.5')
    plot_page(can, histos, 'mass_small_bin_wbin3', lab, output_pdfname,
              xtitle='K^{ +}K^{ -} (GeV)', ytitle='counts', title='W > 2.5 and W < 3.1')


    can.Print('{}]'.format(output_pdfname))
