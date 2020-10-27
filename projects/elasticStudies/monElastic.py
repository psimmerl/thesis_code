from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLatex, TLine
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kBlack

from ROOT import Fit
from ROOT import Math

import numpy as np

from array import array
import math
import sys,os

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)
                     

acc_beam_charge = { 5306: 417438.34,
                    5418: 51939.55,
                    5419: 41387.21,                    
                    5443: 56469.41,
                    5444: 118864.08,
                    5453: 236320.25
                    }

print(' List of accum beam charge')




def rebinHist(title, rebinX=1, rebinY=1, start_range = 1, end_range = 7):
    for ii in range(start_range,end_range):        
        if title.format(ii) in hhs:
            hhs[title.format(ii)].RebinX(rebinX)
            hhs[title.format(ii)].RebinY(rebinY)


def plot_page(canvas, histos, histo_title, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False):
    
    canvas.Clear() 
        
    if isinstance(histos[histo_title], TH1F):
        histos[histo_title].Draw()
    elif isinstance(histos[histo_title], TH2F):
        histos[histo_title].Draw('colz')
        if log:
            gPad.SetLogz() 
    else:
        #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
        pass
    
    if title:
        histos[histo_title].SetTitle('')
        label.DrawLatex(0.1, 0.925, title)

    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    canvas.Print(save_name)


def make_plots(can, canx, cany, fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[], lines_to_draw=[]):
    can.Clear()
    can.Divide(canx,cany)
    cc=1
    root_is_dumb = []
    for hh in hists:
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.Draw('colz')
        gPad.SetLogz()


        if len(lines_to_draw) > 0 :
            for ll in lines_to_draw:
                ll.SetLineColor(kRed)
                ll.Draw('same')
                root_is_dumb.append(ll)

        root_is_dumb.append(hh)
        
        cc+=1
    can.Print(fname)

def plots_overlap(can, histos, title_formatter, label, save_name, title=None, xtitle=None, ytitle=None,
                  x_range=None, y_range=None, rebinX=None,  alternative_counter=None, min_range=0, max_range=6, 
                  y_fit_range=None, scale=None, vline=None):
    print(' making plots in plot_overlap')
    can.Clear()
    can.Divide(1,1)
    can.cd(1)
    root_is_dumb = []
    leg = TLegend(0.1,0.7,0.48,0.9)
    for ii in range(min_range,max_range):
        i = ii
        if alternative_counter:
            if ii < len(alternative_counter):
                i = alternative_counter[ii]                

        if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
            if y_fit_range:
                fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])
                
            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
            if y_range:
                histos.get(title_formatter.format(i), default_histo).GetYaxis().SetRangeUser(y_range[0], y_range[1])

                                
            #histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
            histos.get(title_formatter.format(i), default_histo).SetLineColor(ii+1)            
            leg.AddEntry(histos.get(title_formatter.format(i), default_histo),title+' '+str(i),"l");
            histos.get(title_formatter.format(i), default_histo).Draw('same')

            if y_fit_range:
                label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
            if rebinX:
                histos.get(title_formatter.format(i), default_histo).RebinX(rebinX)

            if scale:
                histos.get(title_formatter.format(i), default_histo).Scale(scale)

            if vline:
                can.Update()
                line = TLine(vline, float(y_range[0]),
                             vline, float(y_range[1]))
                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line.Draw('same')
                root_is_dumb.append(line)



        if title:
            histos.get(title_formatter.format(i), default_histo).SetTitle('')            
            label.DrawLatex(0.1, 0.925, title )

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    leg.Draw()
    can.Print(save_name)



def fit_slices(histo, x_range, x_bin_step, y_fit_range, y_fit_special):
    print(' Status of special fit is ')
    print(y_fit_special)
    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    root_trash_dump=[]
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)
        fit = None
        fit_gaus = None
        fit_pol2 = None
        #fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus')
        #fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
        #fit.SetName(histo.GetTitle() + '_fit{}'.format(i))
        
        if y_fit_range and not y_fit_special:
            print(y_fit_range)
            print(y_fit_special)
            fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus(0)+pol2(3)',  y_fit_range[0], y_fit_range[1])
            fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
            fit.SetName(histo.GetTitle() + '_fit{}'.format(i))        
            fit.SetParameter(1, 1.0 )#0.5 * (y_fit_range[0] + y_fit_range[1]))
            fit.SetParameter(0, histo.GetMaximum())
            fit.SetParameter(2, 0.0869)
            fit.SetParameter(3,10)
            fit.SetParameter(4,-10)
            fit.SetParameter(5,10)            
            projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])

            fit_gaus = TF1('gaus_'+histo.GetTitle() + '_fit{}'.format(i), 'gaus', y_fit_range[0], y_fit_range[1])
            fit_gaus.SetParameter(0, fit.GetParameter(0))
            fit_gaus.SetParameter(1, fit.GetParameter(1))
            fit_gaus.SetParameter(2, fit.GetParameter(2))
            fit_gaus.SetLineColor(kBlue)
            #fit_gaus.Draw('same')

            fit_pol2 = TF1('pol2_'+histo.GetTitle() + '_fit{}'.format(i), 'pol2', y_fit_range[0], y_fit_range[1])
            fit_pol2.SetParameter(0, fit.GetParameter(3))
            fit_pol2.SetParameter(1, fit.GetParameter(4))
            fit_pol2.SetParameter(2, fit.GetParameter(5))
            fit_pol2.SetLineColor(kGreen)
            #fit_pol2.Draw('same')
            root_trash_dump.append(fit_gaus)
            root_trash_dump.append(fit_pol2)
            

        if y_fit_special:
            print(" OMG USING SPECIAL FIT FUNCTION " )
            fit = getFitFractionalHeight(projec,histo.GetTitle() + 'special_fit{}'.format(i), 0.6, fits)

        #projec.Fit(fit, 'R', '',  y_fit_range[0], y_fit_range[1])

        slices.append(projec)
        root_trash_dump.append(fit)
        fits.append([fit, fit_gaus, fit_pol2])

        x_low = histo.GetXaxis().GetBinCenter(x_bin)
        x_high = histo.GetXaxis().GetBinCenter(x_bin + x_bin_step)
        x_values.append(0.5 * (x_high + x_low))

    means = array('d')
    means_err = array('d')
    stds = array('d')
    stds_err = array('d')
    zeros = array('d')

    for f in fits:
        means.append( round(f[0].GetParameter(1),3) )
        means_err.append(f[0].GetParError(1))
        stds.append(f[0].GetParameter(2))
        stds_err.append(f[0].GetParError(1)) ### change if want the sigma of dist.
        zeros.append(0.0)

    graph = TGraphErrors(len(x_values), x_values, means, zeros, stds_err)
    graph.SetName('g_' + histo.GetName())
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(kRed)

    htitle = histo.GetName()
    # write out to txt files here
    #f_txt_out = open(htitle+'_bcsim.txt','w')
    #for ii in range(0, len(x_values) ):
    #    f_txt_out.write(str(x_values[ii]) + ' ' + str(means[ii]) + ' ' + str(stds[ii]) + '\n' )
    #f_txt_out.close()
        

    ## find the x-min and x-max values to use for fitting range limits 
    ## of the delta p or theta vs p or theta
    min_i = next((i for i, x in enumerate(means) if x), None) # x!= 0 for strict match
    max_i = [i for i, e in enumerate(means) if e != 0]
    x_fit_range_min=0
    x_fit_range_max=0    
    if False:# min_i != None and max_i != 0 and len(x_values)>0 :        
        x_fit_range_min=x_values[min_i]
        x_fit_range_max=x_values[max_i[-1]]
            
    return graph, slices, fits, x_fit_range_min, x_fit_range_max


def plot_fits(canvas, histos, x_range, x_bin_step, title_formatter,
              save_name, label, y_fit_range, y_range=None,
              title=None, xtitle=None, ytitle=None, hline=None, 
              fit_special=None, output_txt=None, start_range=1, end_range=7, manual_fit_par=None, special_titles=None):


    print( ' plot fits function status of fit special is ' )
    print(fit_special)
    canvas.Clear()
    canvas.Divide(3,2)

    root_is_dumb = []

    #add can shift if starting at non traditional value of 1.
    can_shift=0
    if start_range < 0: can_shift=-start_range+1
    if start_range == 0: can_shift=1
    print(' need to shift canvas by {}'.format(can_shift))

    for i in range(start_range,end_range):#len(histos)):
        canvas.cd(i+can_shift)

        title_temp = title_formatter.format(i)
        graph, slices, fits, x_min, x_max = fit_slices(histos.get(title_temp, default_histo2d), x_range, x_bin_step, y_fit_range, fit_special)
        graph.SetMarkerStyle(8)
        graph.SetMarkerSize(1)
        
        if y_range:
            graph.GetHistogram().SetMinimum(y_range[0])
            graph.GetHistogram().SetMaximum(y_range[1])
            graph.Draw('AP')
            root_is_dumb.append(graph)
        else:
            graph.Draw('AP')
            root_is_dumb.append(graph)
            
        #if hline:
        line = TLine(x_range[0], hline, x_range[1], hline)
        #line.SetLineStyle(8)
        line.SetLineWidth(2)
        line.SetLineColor(kRed)
        line.Draw()
        root_is_dumb.append(line)
            
        if title:
            graph.SetTitle("")
            label.DrawLatex(0.1, 0.925, title + str(i))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

        '''
        if output_txt:
            print('>>> trying to fit graph')
            #            if manual_fit_par: 
                x_min, x_max = manual_fit_par[i][0], manual_fit_par[i][1]
                print(' using manual fit parameters {} {}'.format(x_min, x_max))
                poly_fit = TF1(title + str(i), "gaus*pol2",x_min, x_max)                    
            graph.Fit(poly_fit,'R')
            poly_fit.SetLineColor(kRed)
            poly_fit.Draw('same')
        '''
            
            

    canvas.Print(save_name)

    # For the slices 
    slice_can = TCanvas('slice_can', 'slice_can', 1200, 1600)
    #slice_pdfname = title_formatter.split('_{}')[0] + '_slices.pdf'
    slice_pdfname = title_formatter + '_slices.pdf'
    slice_can.Print(slice_pdfname + '[')

    out_items = {}
    for i in range(start_range,end_range): #usually 1,7
        title = title_formatter.format(i)
        graph, slices, fits, x_min, x_max = fit_slices(histos.get(title, default_histo2d), x_range, x_bin_step, y_fit_range, fit_special)

        out_items[i]=[graph,slices,fits]
        # Size of slices page
        nrows = 5
        ncols = int(np.ceil(len(slices) / nrows) + 1)
        slice_can.Clear() 
        slice_can.Divide(ncols, nrows)
        for j, (s,f) in enumerate(zip(slices, fits)):
            slice_can.cd(j+1)
            s.Draw()
            f[0].Draw('same')
            f[1].Draw('same')
            f[2].Draw('same')
            integral_min = f[1].GetParameter(1) - 3.0*f[1].GetParameter(2)
            integral_max = f[1].GetParameter(1) + 3.0*f[1].GetParameter(2)
            integrated_events = f[1].Integral(integral_min, integral_max)            
            lab.DrawLatex(0.15, 0.84, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(f[0].GetParameter(1), f[0].GetParameter(2)))
            lab.DrawLatex(0.15, 0.74, 'Signal: {0:6.4f}'.format((integrated_events)))
            if special_titles:
                lab.DrawLatex(0.15, 0.69, 'BC: {}'.format((special_titles)))
                lab.DrawLatex(0.15, 0.64, 'Sig. to BC: {0:6.10f}'.format((integrated_events/special_titles)))

            
        slice_can.Print(slice_pdfname)
    slice_can.Print(slice_pdfname + ']')
    return out_items



def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False,
                     y_fit_range=None, landscape=False, x_range=None, y_range=None, vline=None,
                     hline=None, rebinX=None, rebinY=None, scale=None, alternative_counter=None):

    root_garbage_can = []
    
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    
    ii=0
    for ii in range(0,6):
        i = ii
        
        if alternative_counter:
            if i <= len(alternative_counter):
                i = alternative_counter[ii]
                #print i
        canvas.cd(ii+1)

        #if isinstance(histos[title_formatter.format(i)], TH1F):
        if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
            if y_fit_range:
                fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])
                
            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                                
            #histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
            histos.get(title_formatter.format(i), default_histo).Draw()

            if y_fit_range:
                label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
            if rebinX:
                histos.get(title_formatter.format(i), default_histo).RebinX(rebinX)
            if rebinY:
                histos.get(title_formatter.format(i), default_histo).RebinX(rebinY)

            if scale:
                if len(scale) > 0 and ii < len(scale):
                    histos.get(title_formatter.format(i), default_histo).Scale(scale[ii])
                else:
                    histos.get(title_formatter.format(i), default_histo).Scale(scale[0])
            
        #elif isinstance(histos[title_formatter.format(i)], TH2F):
        elif isinstance(histos.get(title_formatter.format(i), default_histo), TH2F):
            # histos[title_formatter.format(i)].Draw('colz')
            histos.get(title_formatter.format(i), default_histo).Draw('colz')
            if log:
                gPad.SetLogz() 

            if y_range:
                histos.get(title_formatter.format(i), default_histo).GetYaxis().SetRangeUser(y_range[0], y_range[1])

        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass

        if vline:
            line = TLine(vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmin(),
                         vline, histos.get(title_formatter.format(i), default_histo).GetYaxis().GetXmax())
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)

        if hline:
            xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
            xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 
            line = TLine(xmin, hline, xmax, hline)
            line.SetLineColor(1)
            line.SetLineStyle(1)
            line.Draw('same')
            root_garbage_can.append(line)
            
        if title:
            histos.get(title_formatter.format(i), default_histo).SetTitle('')            
            label.DrawLatex(0.1, 0.925, title + ' ' + str(i))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)


def plot_ratio( canvas, histos, title_formatter, label, save_name,
                xtitle=None, ytitle=None, title=None, log=False, min_range=0, max_range=6,
                y_fit_range=None, landscape=False, x_range=None, y_range=None, vline=None,
                hline=None, rebinX=None, rebinY=None, scale=None, alternative_counter=None, plot_ratio_against=None):

    can.Clear()
    can.Divide(1,1)
    can.cd(1)
    root_is_dumb = []
    leg = TLegend(0.1,0.7,0.48,0.9)
    
    for ii in range(min_range,max_range):        
        i = ii
        if alternative_counter:
            if ii < len(alternative_counter):
                i = alternative_counter[ii]                

        if plot_ratio_against:
            if i in plot_ratio_against: continue
            print(' plotting {}'.format(i) )
            if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
                if y_fit_range:
                    fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                    histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])
                
                if scale:
                    print(' Scaling histogram by {}'.format(scale[i]))
                    histos.get(title_formatter.format(i), default_histo).Scale(scale[i])
                    histos.get(title_formatter.format(plot_ratio_against[0], default_histo)).Scale(scale[plot_ratio_against[0]])
            

                temp_xmin = histos.get(title_formatter.format(i)).GetXaxis().GetXmin()
                temp_xmax = histos.get(title_formatter.format(i)).GetXaxis().GetXmax() 

                h_temp = TH1F('temp_'+histos.get(title_formatter.format(i), default_histo).GetTitle(),'temp_'+histos.get(title_formatter.format(i), default_histo).GetTitle(), histos.get(title_formatter.format(i), default_histo).GetXaxis().GetNbins(), temp_xmin, temp_xmax)
                h_temp.Add(histos.get(title_formatter.format(i), default_histo),1)
                            
                h_temp.Divide(histos.get(title_formatter.format(plot_ratio_against[0]), default_histo))
                h_temp.SetLineColor(ii+1)            
                leg.AddEntry(histos.get(title_formatter.format(i), default_histo),title+' '+str(i),"l");
                h_temp.SetTitle('')
                h_temp.Draw('same')
                root_is_dumb.append(h_temp)

                if x_range:
                    histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])

                if y_range:
                    #histos.get(title_formatter.format(i), default_histo).GetYaxis().SetRangeUser(y_range[0], y_range[1])
                    h_temp.GetYaxis().SetRangeUser(y_range[0], y_range[1])

                if y_fit_range:
                    label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
                if rebinX:
                    histos.get(title_formatter.format(i), default_histo).RebinX(rebinX)
                    

                if vline:
                    can.Update()
                    line = TLine(vline, float(y_range[0]),
                             vline, float(y_range[1]))
                    line.SetLineColor(kRed)
                    line.SetLineStyle(1)
                    line.Draw('same')
                    root_is_dumb.append(line)



                if title:
                    histos.get(title_formatter.format(i), default_histo).SetTitle('')            
                    label.DrawLatex(0.1, 0.925, title )
                    
                if xtitle:
                    label.DrawLatex(0.5, 0.015, xtitle)

                if ytitle:
                    label.SetTextAngle(90)
                    label.DrawLatex(0.04, 0.5, ytitle)
                    label.SetTextAngle(0)

    leg.Draw()
    root_is_dumb.append(leg)
    can.Print(save_name)


    
print('opening file')

datatype = sys.argv[2]
hhs={}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj
    

base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_elasticXS_'+datatype+'.pdf'
can = TCanvas('can','can',1200,600)
can.Print('{}['.format(mon_out_file_name))
can.Divide(2,1)

lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)


gStyle.SetOptStat(00000)

#########################################################
#########################################################
# code to compare reconstructed DATA to reconstructed SIM
#########################################################
#########################################################

# compare electron kinematics - p, theta,

### runs 
runs = [5418, 5419, 5306, 5443, 5444, 5453]
n_entries = {}
for ii in runs:
    n_entries[ii] = hhs['hw_el_pass_all_r'+str(ii)].GetEntries()

print n_entries




'''
plot_sector_page(can, hhs,'hptheta_el_pass_all{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) P vs #theta S', xtitle='p (GeV)', ytitle='#theta (deg)', log=True)
plot_sector_page(can, hhs,'hvztheta_el_pass_all{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) Vz vs #theta S', xtitle='vz (cm)', ytitle='#theta (deg)', log=True)


plot_page(can, hhs, "hw_el_pass_all", lab, save_name=mon_out_file_name, 
          title='El (FD) W ep->eX', xtitle='W (GeV)' )

plot_sector_page(can, hhs,'hw_el_pass_all_{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W S ', xtitle='W (GeV)')#scale=1541.8455)

plot_sector_page(can, hhs,'hw_small_el_pass_all_{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - zoom S', xtitle='W (GeV)', rebinX=2)

plot_sector_page(can, hhs,'hw_small_el_pass_all_{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - zoom S', xtitle='W (GeV)', rebinX=2)

plot_sector_page(can, hhs,'hw_el_pass_all_0_r{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - Sector 1 Run', xtitle='W (GeV)', alternative_counter=runs)

plot_sector_page(can, hhs,'hw_el_pass_all_1_r{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - Sector 2 Run', xtitle='W (GeV)', alternative_counter=runs, vline = 0.938)

plot_sector_page(can, hhs,'hw_el_pass_all_2_r{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - Sector 3 Run', xtitle='W (GeV)', alternative_counter=runs, vline = 0.938)

plot_sector_page(can, hhs,'hw_el_pass_all_3_r{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - Sector 4 Run', xtitle='W (GeV)', alternative_counter=runs, vline = 0.938)

plot_sector_page(can, hhs,'hw_el_pass_all_4_r{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - Sector 5 Run', xtitle='W (GeV)', alternative_counter=runs, vline = 0.938)

plot_sector_page(can, hhs,'hw_el_pass_all_5_r{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - Sector 6 Run', xtitle='W (GeV)', alternative_counter=runs, vline = 0.938)


plots_overlap(can,hhs,'hw_el_pass_all_0_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 1 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

plots_overlap(can,hhs,'hw_el_pass_all_1_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 2 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

plots_overlap(can,hhs,'hw_el_pass_all_2_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 3 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

plots_overlap(can,hhs,'hw_el_pass_all_3_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 4 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

plots_overlap(can,hhs,'hw_el_pass_all_4_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 5 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

plots_overlap(can,hhs,'hw_el_pass_all_5_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 6 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

# plot the indivual runs scaled by their beam charge current 
plots_overlap(can,hhs,'hw_small_el_pass_all_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 1 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

# plot it for different binning size
plots_overlap(can,hhs,'hw_el_pass_all_r{}', lab, save_name=mon_out_file_name,
              title='El. (FD) W - Sector 1 All Runs', xtitle='W (GeV)', ytitle='counts', alternative_counter=runs, y_range=[0,1000], vline=0.938)

'''
#plot_ratio(can,hhs,'hw_el_pass_all_r{}', lab, save_name=mon_out_file_name, title='El (FD) W - R5418 and R  ', xtitle='W (GeV)', ytitle='ratio', alternative_counter=runs, plot_ratio_against = [5418, 5306, 5444, 5453], scale=acc_beam_charge, y_range=[0,10] )

y_lims = {5418: 70,
          5419: 100,
          5306: 300,
          5443: 1000,
          5444: 1500,
          5453: 60 }
for rr in runs:
    plots_overlap(can,hhs,'hw_el_pass_all_{}_r'+str(rr), lab, save_name=mon_out_file_name,
                  title='El. (FD) W - Run '+str(rr)+' Sector ', xtitle='W (GeV)', ytitle='counts', y_range=[0,y_lims[rr]], vline=0.938)


############# now look at the 2D W vs Theta dist
plot_page(can, hhs, "hwtheta_el_pass_all", lab, save_name=mon_out_file_name, 
          title='El (FD) W vs #Theta_{el} ep->eX All Runs All Sector', xtitle='#Theta_{el} (deg)', ytitle='W (GeV)' )

for rr in runs:
    plot_sector_page(can, hhs, "hwtheta_el_pass_all_{}_r"+str(rr), lab, save_name=mon_out_file_name, 
                     title='El (FD) W vs #Theta_{el} ep->eX Run ' + str(rr) + ' Sector', xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', rebinY=8 )
'''
rebinHist('hwtheta_el_pass_all_{}_r5418', rebinY=12, start_range=0, end_range=6)
fr_5418 = plot_fits(can,hhs,x_range=[7,13],x_bin_step=1, title_formatter='hwtheta_el_pass_all_{}_r5418', save_name=mon_out_file_name, label=lab, y_fit_range=[0.8,1.35], y_range=[0.60, 1.35], xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', title='Elastic W vs #Theta_{el} Run 5418 Sector ', hline = 0.938272, special_titles=acc_beam_charge[5418])


rebinHist('hwtheta_el_pass_all_{}_r5419', rebinY=12, start_range=0, end_range=6)
fr_5419 = plot_fits(can,hhs,x_range=[7,13],x_bin_step=1, title_formatter='hwtheta_el_pass_all_{}_r5419', save_name=mon_out_file_name, label=lab, y_fit_range=[0.8,1.35], y_range=[0.60, 1.35], xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', title='Elastic W vs #Theta_{el} Run 5419 Sector ', hline = 0.938272, special_titles=acc_beam_charge[5419])
'''
rebinHist('hwtheta_el_pass_all_{}_r5306', rebinY=12, start_range=0, end_range=7)
fr_5306 = plot_fits(can,hhs,x_range=[7,13],x_bin_step=1, title_formatter='hwtheta_el_pass_all_{}_r5306', save_name=mon_out_file_name, label=lab, y_fit_range=[0.8,1.35], y_range=[0.60, 1.35], xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', title='Elastic W vs #Theta_{el} Run 5306 Sector ', hline = 0.938272, special_titles=acc_beam_charge[5306])

rebinHist('hwtheta_el_pass_all_{}_r5443', rebinY=12, start_range=0, end_range=7)
fr_5443 = plot_fits(can,hhs,x_range=[5,10],x_bin_step=1, title_formatter='hwtheta_el_pass_all_{}_r5443', save_name=mon_out_file_name, label=lab, y_fit_range=[0.8,1.35], y_range=[0.60, 1.35], xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', title='Elastic W vs #Theta_{el} Run 5443 Sector ', hline = 0.938272, special_titles=acc_beam_charge[5443])

rebinHist('hwtheta_el_pass_all_{}_r5444', rebinY=12, start_range=0, end_range=7)
fr_5444 = plot_fits(can,hhs,x_range=[5,10],x_bin_step=1, title_formatter='hwtheta_el_pass_all_{}_r5444', save_name=mon_out_file_name, label=lab, y_fit_range=[0.8,1.35], y_range=[0.60, 1.35], xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', title='Elastic W vs #Theta_{el} Run 5444 Sector ', hline = 0.938272, special_titles=acc_beam_charge[5444])
'''
rebinHist('hwtheta_el_pass_all_{}_r5453', rebinY=12, start_range=0, end_range=6)
fr_5453 = plot_fits(can,hhs,x_range=[5,10],x_bin_step=1, title_formatter='hwtheta_el_pass_all_{}_r5453', save_name=mon_out_file_name, label=lab, y_fit_range=[0.8,1.35], y_range=[0.60, 1.35], xtitle='#Theta_{el} (deg)', ytitle='W (GeV)', title='Elastic W vs #Theta_{el} Run 5453 Sector ', hline = 0.938272, special_titles=acc_beam_charge[5453])
'''

can.Clear()
can.Divide(3,2)
l_of_runs=[fr_5306, fr_5444, fr_5443]
rrr = [5306, 5444,5443]
root_dump = []

delta_theta=0.5
colors = [kRed, kBlack, kBlue]
for ii in range(1,7):
    can.cd(ii)
    mg_nom_rates = TMultiGraph()
    graphs = []
    cc=0
    for fit_results in l_of_runs:
        gr, sl, ft = fit_results[ii]# = fr_5306[ii]
        theta_center=5
        cnt = array('d')
        val = array('d')
        cnt_err = array('d')
        val_err = array('d')
        
        for ff in ft:
            theta_center+=delta_theta
            charge_nom = acc_beam_charge[rrr[cc]]
            integral_min = ff[1].GetParameter(1) - 3.0*ff[1].GetParameter(2)
            integral_max = ff[1].GetParameter(1) + 3.0*ff[1].GetParameter(2)
            integrated_events = ff[1].Integral(integral_min, integral_max)            
            #print(' integral min. {}, integral max. {} - events {} '.format(integral_min, integral_max, integrated_events))
            print(rrr[cc])
            print(' theta center {} , max events {}'.format(theta_center, integrated_events))
            if integrated_events/charge_nom < 0 or integrated_events/charge_nom  > 0.001 : continue

            cnt.append(theta_center)
            val.append(integrated_events/charge_nom)
            cnt_err.append(0)
            val_err.append(0)

        g_run=TGraphErrors(len(cnt), cnt, val, cnt_err, val_err)
        g_run.SetMarkerStyle(8)
        g_run.SetMarkerSize(1)
        g_run.SetMarkerColor(colors[cc])
        g_run.SetTitle('Run {} S {}'.format(rrr[cc],ii))
        #mg_nom_rates.SetTitle('Run {} Sector {}'.format(rrr[cc],ii))
        mg_nom_rates.Add(g_run,)
        cc+=1
    mg_nom_rates.SetTitle('Integrated Events/Beam Charge Sector {}'.format(ii))
    mg_nom_rates.Draw('AP')
    gPad.SetLogy()
    root_dump.append(mg_nom_rates)
   
can.Print(mon_out_file_name)








can.Print('{}]'.format(mon_out_file_name))















    
    
    
    
