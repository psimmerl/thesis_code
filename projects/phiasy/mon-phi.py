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
                  gPad, gStyle, TLatex, TGraphErrors, TLegend, TEllipse,
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
                      title=None, xtitle=None, ytitle=None , plot_details=None, field=None):


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

            if bb-1 not in plot_details['bins_to_draw'] : continue
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
                

            if title :
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
                #if(den==0):
                if( P < 6 or N < 6 ):
                    continue
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
            if bb-1 in plot_details['bins_to_draw'] :

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

            fit_asy = None
            if 'special_bsa_fit' in plot_details.keys():
                fit_asy = TF1('fasy{}'.format(bb),plot_details['special_bsa_fit'][bb-1], -180, 180)
            else:
                fit_asy = TF1('fasy{}'.format(bb),"[0]*sin(x*TMath::Pi()/180.0)", -180, 180)

            fit_asy.SetParameter(0,1)
            graph.Fit('fasy{}'.format(bb))

            fit_asy.SetLineColor(kRed)
            fit_asy_chi2 = fit_asy.GetChisquare()
            fit_asy_ndf = fit_asy.GetNDF()
            root_is_dumb.append(fit_asy)
            if f_out:
                fout.write(str(bc) + ' ' + str(fit_asy.GetParameter(0)) + ' ' + str(fit_asy.GetParError(0)) + ' ' + str(fit_asy_chi2) + ' ' + str(fit_asy_ndf) + ' \n')
            
            
            # draw bin related info
            if bb-1 not in plot_details['bins_to_draw'] : continue

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

            fit_ptr = graph.Fit('fasy{}'.format(bb),'S')
            correlation_matrix = fit_ptr.GetCovarianceMatrix()
            correlation_matrix.Print()
            ## draw the correlation matrix elemenst
            if correlation_matrix.GetNcols() > 1:
                can_cor = TCanvas('can_cor'+f_out2,'can_cor'+f_out2,800,800)
                can_cor.Divide(1,1)
                can_cor.cd(1)
                a = correlation_matrix(0,0)
                b = correlation_matrix(0,1)
                c = correlation_matrix(1,1)
                eigen_val1 = (a+c)/2 + math.sqrt( ((a-c)/2)**2 + b**2 )
                eigen_val2 = (a+c)/2 - math.sqrt( ((a-c)/2)**2 + b**2 )
                ellipse_angle = math.atan2(eigen_val1-correlation_matrix(0,0),correlation_matrix(0,1))
                print('------------> lambda1 {}  and lambda2 {} and ellipse angle {}'.format(eigen_val1,eigen_val2, ellipse_angle*180.0/3.1415926))
                correlation_ellipse = TEllipse(0,0, math.sqrt(eigen_val1), math.sqrt(eigen_val2), 0, 360, ellipse_angle*180.0/3.1415926)
                #print('{}'.format( correlation_matrix(1,1)))
                h1 = gPad.DrawFrame(-0.03, -math.sqrt(eigen_val1)*1.05, 0.03, math.sqrt(eigen_val1)*1.05)
                h1.GetXaxis().SetNdivisions(-10)
                h1.GetYaxis().SetNdivisions(-10)
                h1.GetXaxis().SetLabelSize(0.025)
                h1.GetYaxis().SetLabelSize(0.025)
                gPad.SetGrid()
                correlation_ellipse.SetLineWidth(3)
                correlation_ellipse.Draw("same")
                gPad.RedrawAxis("g")
                label.DrawLatex(0.1, 0.925, title + ',Signal, Correlation of #lambda_{1},#lambda_{2}')
                label.DrawLatex(0.5, 0.025,'#lambda_{2}')
                label.SetTextAngle(90)
                label.DrawLatex(0.045, 0.5,'#lambda_{1}') # was 0.035
                label.SetTextAngle(0)
                can_cor.SaveAs('asy_vs_trento_method1_imkpkm_bin{}_correlation'.format(bb)+f_out2+'.pdf')
                root_is_dumb.append(can_cor)
                root_is_dumb.append(label)
                root_is_dumb.append(correlation_ellipse)


            asy_per_mass.append([bc, fit_asy.GetParameter(0), fit_asy.GetParError(0)])
            

    

    if f_out:
        fout.close()
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

    field_setting = 'inbNoutb'
    if 'inbNout' in output_pdfname:
        field_setting="inbNoutb"
    elif 'inb' in output_pdfname:
        field_setting="inb"
    elif 'outb' in output_pdfname:
        field_setting="outb"
    else:
        field_setting=""

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


    bp = 0
    max_xaxis = 0
    binvar='binVar11'
    if field_setting == 'inb':
        bp = 0.8592 #0.869
        max_xaxis = 1.0908
    elif field_setting == 'outb':
        bp = 0.892 # outb should have same bp as the end of inbending data run
        max_xaxis = 1.1049
    elif field_setting == 'inbNoutb':
        bp = 0.874
        max_xaxis = 1.09782

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

    
    phi_mass_sliced = { 'max_bins': 5, #len(mass_bins),
                        'colors': True,
                        'max_xaxis': max_xaxis
                    }

    phi_mass_overlap = { 'max_bins': 5, #len(overlap_mass_bins),
                         'colors': True,
                         'sliding_bin_list':True,
                         'max_xaxis': max_xaxis
                     }

    
    plot_slices( can, histos, 'mass_small_bin{}' , output_pdfname, lab,
                      title = 'Final Events, Inv. Mass Bins:', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_sliced, field=field_setting)
    # put back FD later 
    
    im_asy_details={'max_bins': 5, #len(mass_bins),
                    'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                    'bin_list': mass_bins,          
                    'bins_to_draw':[0,1,2]          
                }

    
    plot_page_bins(can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=im_asy_details, field=field_setting)


    plot_asy( can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
              title='Asy: Comb. Datasets,  ', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details, 
              f_out='asy_vs_imkpkm_pidtype4a_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'V2.txt', f_out2='pidtype4a_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'V2.txt', field=field_setting)
    
    '''
    
    ######################
    ## plot the asym for im kpkm bins that overlap to the right of the phi peak
    ## this is used for systematicall checking the sideband subtraction
    plot_slices( can, histos, 'sliding_mass_small_bin{}' , output_pdfname, lab,
                 title='Final Events, Overlap Inv. Mass Bins', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_overlap, field=field_setting)

    # put back FD later     
    im_asy_details_ovrlp={'max_bins': 4, #11,#len(overlap_mass_bins),
                          'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                          'sliding_bin_list': overlap_mass_bins
                      }
                    
    plot_page_bins(can, histos, 'trento_{}_sliding_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=im_asy_details_ovrlp, field=field_setting)

    plot_asy( can, histos, 'trento_{}_sliding_mass_bin{}', output_pdfname, lab,
              title='Asy', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details_ovrlp, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_sliding_'+'Var4'+'_'+field_setting+'.txt', f_out2='pidtype1_additionalcuts_sliding_'+'Var4'+'_'+field_setting+'.txt', field=field_setting)
    '''

    '''
    ########################
    ## systematic check of fitting the asymmetry with additional bsa moments
    ## first check the result using the sin phi and cos phi momentums
    plot_slices( can, histos, 'mass_small_bin{}' , output_pdfname, lab,
                      title = 'Final Events, Inv. Mass Bins:', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_sliced, field=field_setting)
    # put back FD later 
    
    im_asy_details_cosmoment={'max_bins': 5, #len(mass_bins),
                              'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                              'bin_list': mass_bins,                    
                              'special_bsa_fit': ["[0]*sin(x*TMath::Pi()/180.0)","([0]*sin(x*TMath::Pi()/180.0))/(1+[1]*cos(x*TMath::Pi()/180.0))",
                                                  "[0]*sin(x*TMath::Pi()/180.0)","[0]*sin(x*TMath::Pi()/180.0)",
                                                  "[0]*sin(x*TMath::Pi()/180.0)"],
                              'bins_to_draw':[0,1,2]
                          }

    
    plot_page_bins(can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=im_asy_details_cosmoment, field=field_setting)

    plot_asy( can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
              title='Asy: sin(#phi)/1+cos(#phi)', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details_cosmoment, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'_systmodbsa_V2.txt', 
              f_out2='pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'_systmodbsa_V2.txt', field=field_setting)

    #########################
    ## check with 1+cos(2phi) moment
    plot_slices( can, histos, 'mass_small_bin{}' , output_pdfname, lab,
                      title = 'Final Events, Inv. Mass Bins:', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_sliced, field=field_setting)
    # put back FD later 
    
    im_asy_details_cosmoment={'max_bins': 5, #len(mass_bins),
                              'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                              'bin_list': mass_bins,                    
                              'special_bsa_fit': ["[0]*sin(x*TMath::Pi()/180.0)","([0]*sin(x*TMath::Pi()/180.0))/(1+[1]*cos(2*x*TMath::Pi()/180.0))",
                                                  "[0]*sin(x*TMath::Pi()/180.0)","[0]*sin(x*TMath::Pi()/180.0)",
                                                  "[0]*sin(x*TMath::Pi()/180.0)"],
                              'bins_to_draw':[0,1,2]
                          }

    
    plot_page_bins(can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=im_asy_details_cosmoment, field=field_setting)

    plot_asy( can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
              title='Asy: sin(#phi)/1+cos(2#phi)', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details_cosmoment, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'_systmodbsa2_V2.txt', 
              f_out2='pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'_systmodbsa2_V2.txt', field=field_setting)


    ##########################
    ## fit with all moments
    im_asy_details_cosmoment3={'max_bins': 5, #len(mass_bins),
                              'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                              'bin_list': mass_bins,                    
                              'special_bsa_fit': ["[0]*sin(x*TMath::Pi()/180.0)","([0]*sin(x*TMath::Pi()/180.0))/([1]*cos(x*TMath::Pi()/180.0)+[2]*cos(2*x*TMath::Pi()/180.0))",
                                                  "[0]*sin(x*TMath::Pi()/180.0)","[0]*sin(x*TMath::Pi()/180.0)",
                                                  "[0]*sin(x*TMath::Pi()/180.0)"],
                              'bins_to_draw':[0,1,2]
                          }
    
    
    plot_page_bins(can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=im_asy_details_cosmoment3, field=field_setting)

    plot_asy( can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
              title='Asy: sin(#phi)/cos(#phi)+cos(2#phi)', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details_cosmoment3, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'_systmodbsa3_V2.txt', 
              f_out2='pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'_systmodbsa3_V2.txt', field=field_setting)
    '''

    '''

    
    
    ##### plot the asymettry for different bins in -t, xb, q2, w
    ### the results are plotted for each bin in mass 
    ### not many bins in t, q2, and xb so reasonable to do here 
    ## old binning
    
    t_im_asy_details={'max_bins': 5, #len(mass_bins),
                      'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                      'bin_list': mass_bins, 
                      'bins': t_bins,
                      'bin_centers':t_centers,
                      'bins_to_draw':[0,1,2]
                }

    for tt in range(1,len(t_bins)):
        print('------------------------------------------------------------> creating asymettry plots for t bin {}'.format(tt))
        plot_slices( can, histos, 't_bin'+str(tt)+'_mass_small_bin{}' , output_pdfname, lab,
                      title='FD:Final, Inv.MassBins. {0:.3f}<-t<{1:.3f}'.format(t_bins[tt-1],t_bins[tt]), xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , 
                     plot_details=phi_mass_sliced, field=field_setting)

        plot_asy( can, histos, 'trento_{}_t_bin'+str(tt)+'_mass_bin{}', output_pdfname, lab,
                  title='Asy,{0:.3f}<-t<{1:.3f}'.format(t_bins[tt-1],t_bins[tt]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=t_im_asy_details, 
                  f_out='asy_vs_imkpkm_tbin{0}_additionalcuts_{2}_{1}V2.txt'.format(tt,field_setting,binvar), 
                  f_out2='pidtype1_additionalcuts_'+binvar+'_'+field_setting+'_tbin{}V2.txt'.format(tt), field=field_setting)
    
    #######################################
    
    q2_im_asy_details={'max_bins': 5,#len(mass_bins),
                       'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                       'bin_list': mass_bins,
                       'bins': q2_bins,
                       'bin_centers':q2_centers,
                       'bins_to_draw':[0,1,2]                       
                }

    ## q2 bins
    for qq in range(1,len(q2_bins)):
        print('------------------------------------------------------------> creating asymettry plots for Q2 bin {}'.format(qq))
        plot_slices( can, histos, 'q2_bin'+str(qq)+'_mass_small_bin{}' , output_pdfname, lab,
                      title='FD:Final, Inv.MassBins. {0:.3f}<Q2<{1:.3f}'.format(q2_bins[qq-1], q2_bins[qq]), xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , 
                     plot_details=phi_mass_sliced, field=field_setting)

        plot_asy( can, histos, 'trento_{}_q2_bin'+str(qq)+'_mass_bin{}', output_pdfname, lab,
                  title='Asy,{0:.3f}<Q2<{1:.3f}'.format(q2_bins[qq-1], q2_bins[qq]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=q2_im_asy_details,
                  f_out='asy_vs_imkpkm_q2bin{0}_additionalcuts_{2}_{1}V2.txt'.format(qq,field_setting,binvar), 
                  f_out2='pidtype1_additionalcuts_'+binvar+'_'+field_setting+'_q2bin{}V2.txt'.format(tt), field=field_setting)

    #######################################
    xb_im_asy_details={'max_bins': 5,#len(mass_bins),
                       'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                       'bin_list': mass_bins,
                       'bins': xb_bins,
                       'bin_centers': xb_centers,
                       'bins_to_draw':[0,1,2]
                       }
    
    ## xb bins 
    for xx in range(1, len(xb_bins)):
        print('------------------------------------------------------------> creating asymettry plots for Xb bin {}'.format(xx))
        plot_slices( can, histos, 'xb_bin'+str(xx)+'_mass_small_bin{}' , output_pdfname, lab,
                     title='FD:Final, Inv.MassBins. {0:.3f}<Xb<{1:.3f}'.format(xb_bins[xx-1], xb_bins[xx]), xtitle='K^{ +}  K^{ -}  (GeV)', ytitle='counts' , 
                     plot_details=phi_mass_sliced, field=field_setting)

        plot_asy( can, histos, 'trento_{}_xb_bin'+str(xx)+'_mass_bin{}', output_pdfname, lab,
                  title='Asy,{0:.3f}<Xb<{1:.3f}'.format(xb_bins[xx-1], xb_bins[xx]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=xb_im_asy_details,
                  f_out='asy_vs_imkpkm_xbbin{0}_additionalcuts_{2}_{1}V2.txt'.format(xx,field_setting,binvar), 
                  f_out2='pidtype1_additionalcuts_'+binvar+'_'+field_setting+'_xbbin{}V2.txt'.format(tt), field=field_setting)
    


    ####################
    w_im_asy_details={'max_bins': 5, #len(mass_bins),
                      'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                      'bin_list': mass_bins, 
                      'bins': w_bins,
                      'bin_centers': w_centers,
                      #'ignore_bad_bids': True,
                      'bins_to_draw':[0,1,2]
                }

    for ww in range(1,len(w_bins)):
        print('------------------------------------------------------------> creating asymettry plots for w bin {}'.format(ww))
        plot_slices( can, histos, 'w_bin'+str(ww)+'_mass_small_bin{}' , output_pdfname, lab,
                      title='FD:Final, Inv.MassBins. {0:.3f}<W<{1:.3f}'.format(w_bins[ww-1],w_bins[ww]), xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , 
                     plot_details=phi_mass_sliced, field=field_setting)

        plot_page_bins(can, histos, 'trento_{}_w_bin'+str(ww)+'_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=w_im_asy_details, field=field_setting)

        plot_asy( can, histos, 'trento_{}_w_bin'+str(ww)+'_mass_bin{}', output_pdfname, lab,
                  title='Asy,{0:.3f}<W<{1:.3f}'.format(w_bins[ww-1],w_bins[ww]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=w_im_asy_details, 
                  f_out='asy_vs_imkpkm_wbin{0}_additionalcuts_{2}_{1}V2.txt'.format(ww,field_setting,binvar), 
                  f_out2='pidtype1_additionalcuts_'+binvar+'_'+field_setting+'_wbin{}V2.txt'.format(ww), field=field_setting)
    '''

    '''
    # attemping to increase bins
    #######################################
    ## bin the asymmetry into xb and q2 bins 
    xb_q2_im_asy_details={'max_bins': 5,#len(mass_bins),
                          'bp': bp, #0.869, # run weighted asymettry is ~2/3 have 85.9 and the last ~ 1/3  have 89.2 %. weighted average is ~ 86.9 %.
                          'bin_list': mass_bins,
                          'bins': xb_bins,
                          #'bin_centers': xb_centers
                       }
    plot_slices( can, histos, 'xb_bin1_q2_bin1_mass_small_bin{}' , output_pdfname, lab,
                      title = 'Final Events, Inv. Mass Bins, xb vs Q^{2} Bin 1:', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_sliced, field=field_setting)
    
    plot_page_bins(can, histos, 'trento_{}_xb_bin1_q2_bin1_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=xb_q2_im_asy_details, field=field_setting)

    plot_asy( can, histos, 'trento_{}_xb_bin1_q2_bin1_mass_bin{}', output_pdfname, lab,
              title='Asy: Inb/Outb, Xb1, Q^{2}1, ', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=xb_q2_im_asy_details, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_xb_bin1_q2_bin1_'+binvar+'_'+field_setting+'.txt', 
              f_out2='pidtype1_additionalcuts_xb_bin1_q2_bin1_'+binvar+'_'+field_setting+'.txt', field=field_setting)


    ###############3
    ## high xb, q2 bins
    plot_slices( can, histos, 'xb_bin2_q2_bin2_mass_small_bin{}' , output_pdfname, lab,
                      title = 'Final Events, Inv. Mass Bins, xb vs Q^{2} Bin 2:', xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , plot_details=phi_mass_sliced, field=field_setting)
    
    plot_page_bins(can, histos, 'trento_{}_xb_bin2_q2_bin2_mass_bin{}', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=xb_q2_im_asy_details, field=field_setting)

    plot_asy( can, histos, 'trento_{}_xb_bin2_q2_bin2_mass_bin{}', output_pdfname, lab,
              title='Asy: Inb/Outb, Xb2, Q^{2}2,   ', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=xb_q2_im_asy_details, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_xb_bin2_q2_bin2_'+binvar+'_'+field_setting+'.txt', 
              f_out2='pidtype1_additionalcuts_xb_bin2_q2_bin2_'+binvar+'_'+field_setting+'.txt', field=field_setting)
    

    
    ## low q2 so q2 < 2.75
    for ww in range(1,len(w_bins)):
        print('------------------------------------------------------------> creating asymettry plots for w bin {}'.format(ww))
        plot_slices( can, histos, 'w_bin'+str(ww)+'_mass_small_bin{}_lowq2' , output_pdfname, lab,
                      title='FD:Final, Inv.MassBins. {0:.3f}<W<{1:.3f}'.format(w_bins[ww-1],w_bins[ww]), xtitle='K^{ +}  K^{ -}   (GeV)', ytitle='counts' , 
                     plot_details=phi_mass_sliced, field=field_setting)

        plot_page_bins(can, histos, 'trento_{}_w_bin'+str(ww)+'_mass_bin{}_lowq2', output_pdfname, lab,
                   title='#phi_{trento}', xtitle='#phi_{trento} (deg)', ytitle='counts' , plot_details=w_im_asy_details, field=field_setting)

        plot_asy( can, histos, 'trento_{}_w_bin'+str(ww)+'_mass_bin{}_lowq2', output_pdfname, lab,
                  title='Asy,{0:.3f}<W<{1:.3f}'.format(w_bins[ww-1],w_bins[ww]), xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=w_im_asy_details, 
                  f_out='asy_vs_imkpkm_wbin{0}_lowq2_additionalcuts_{2}_{1}.txt'.format(ww,field_setting,binvar), 
                  f_out2='pidtype1_additionalcuts_'+binvar+'_'+field_setting+'_wbin{}_lowq2.txt'.format(ww), field=field_setting)

    '''
    '''
    #######################
    ## determine bsa using BP calculated differently
    im_asy_details={'max_bins': 5, #len(mass_bins),
                    'bp': 0.879, # method 1 to calc bp in thesis
                    'bin_list': mass_bins,          
                    'bins_to_draw':[0,1,2]          
                }

    plot_asy( can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
              title='Asy: Comb. Datasets,  ', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'V2_bpm1.txt', 
              f_out2='pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'V2_bpm1.txt', field=field_setting)


    ## using 88.5% bp from m2 in thesis 
    im_asy_details={'max_bins': 5, #len(mass_bins),
                    'bp': 0.885, # method 1 to calc bp in thesis
                    'bin_list': mass_bins,          
                    'bins_to_draw':[0,1,2]          
                }

    plot_asy( can, histos, 'trento_{}_mass_bin{}', output_pdfname, lab,
              title='Asy: Comb. Datasets,  ', xtitle='#phi_{trento} (deg)', ytitle='Asy' , plot_details=im_asy_details, 
              f_out='asy_vs_imkpkm_pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'V2_bpm2.txt', 
              f_out2='pidtype1_additionalcuts_rmvres12_'+binvar+'_'+field_setting+'V2_bpm2.txt', field=field_setting)


    '''
    can.Print('{}]'.format(output_pdfname))
