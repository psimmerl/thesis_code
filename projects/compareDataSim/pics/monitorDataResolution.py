from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile
from ROOT import TGraphErrors, TBox
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import kRed
from array import array
import math
import sys

def fit_slices(histo, x_range, x_bin_step):

    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    histo.RebinY(4)
    
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)
        fit = TF1(histo.GetTitle()+'_fit{}'.format(i), 'gaus')

        fit_max = projec.GetMean() + 4.0*projec.GetRMS()
        fit_min = projec.GetMean() - 4.0*projec.GetRMS()

        
        projec.Fit(fit,'R','',fit_min,fit_max)

        slices.append(projec)
        fits.append(fit)

        x_low = histo.GetXaxis().GetBinCenter(x_bin)
        x_high = histo.GetXaxis().GetBinCenter(x_bin + x_bin_step)
        x_values.append(0.5 * (x_high + x_low))

    means = array('d')
    means_err = array('d')
    stds = array('d')
    stds_err = array('d')
    zeros = array('d')

    for f in fits:
        means.append(f.GetParameter(1))
        means_err.append(f.GetParError(1))
        stds.append(f.GetParameter(2))
        stds_err.append(f.GetParError(2))
        zeros.append(0.0)

    graph = TGraphErrors(len(x_values), x_values, means, zeros, stds)
    graph.SetName('g_' + histo.GetName())

    return graph, slices, fits


def make_plots(can, fname, hists, x_ranges, rebinx, min_y, max_y):
    can.Clear()
    can.Divide(3,1)
    cc=1
    root_is_dumb = []
    for hh in hists:
        can.cd(cc)
        g_fit, sl_fit, f_fit  = fit_slices(hh,x_ranges[cc-1],rebinx)
        g_fit.SetTitle(hh.GetTitle()+'_slices')
        g_fit.SetMarkerStyle(21)
        g_fit.SetMarkerSize(1)
        g_fit.Draw('AP')
        g_fit.GetHistogram().SetMaximum(max_y)
        g_fit.GetHistogram().SetMinimum(min_y)
        g_fit.Draw('AP')
        root_is_dumb.append(g_fit)
        cc+=1
    can.Print(fname)

def make_slice_plots(incan, fname, hist, ranges, can_x, can_y, rebinx ):
    
    sl_container=[]
        
    grph, slices, fits = fit_slices(hist, ranges, rebinx)
    incan.Clear()
    incan.Divide(4,4)

    ss=1
    for sl in slices:
        incan.cd(ss)
        sl.SetTitle(sl.GetTitle()+'slice'+str(ss))
        sl.Draw()
        sl_container.append(sl)
        ss+=1
    incan.Print(fname)


def define2DHist(htemp, nameX, nameY):
    htemp.SetTitleX(nameX)
    htemp.SetTitleY(nameY)
    
fin = TFile(sys.argv[1])
tag = sys.argv[2]

my_pdf_name=data_type+'res_phi_'+tag+'.pdf'
my_pdf_fits=data_type+'res_fits_phi_'+tag+'.pdf'
my_pdf_slices=data_type+'res_phi_slices_'+tag+'.pdf'

hhs = {}
for kk in fin.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj


can =TCanvas('can','can',1200,600)

#open pdf file 
can.Print('{}['.format(my_pdf_name))
can.Divide(2,1)
can.cd(1)
hhs['epX_vs_theta_pr'].Draw('colz')
can.cd(2)
hhs['epX_vs_p_pr'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
hhs['epkpX_vs_theta_km'].Draw('colz')
can.cd(2)
hhs['epkpX_vs_p_km'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
hhs['epkmX_vs_theta_kp'].Draw('colz')
can.cd(2)
hhs['epkmX_vs_p_kp'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
hhs['ekpkmX_vs_theta_pr'].Draw('colz')
can.cd(2)
hhs['ekpkmX_vs_p_pr'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

#deltas
can.Divide(2,2)
can.cd(1)
hhs['delta_theta_vs_theta_pr'].Draw('colz')
can.cd(2)
hhs['delta_theta_vs_p_pr'].Draw('colz')
can.cd(3)
hhs['delta_p_vs_theta_pr'].Draw('colz')
can.cd(4)
hhs['delta_p_vs_p_pr'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,2)
can.cd(1)
hhs['delta_theta_vs_theta_kp'].Draw('colz')
can.cd(2)
hhs['delta_theta_vs_p_kp'].Draw('colz')
can.cd(3)
hhs['delta_p_vs_theta_kp'].Draw('colz')
can.cd(4)
hhs['delta_p_vs_p_kp'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,2)
can.cd(1)
hhs['delta_theta_vs_theta_km'].Draw('colz')
can.cd(2)
hhs['delta_theta_vs_p_km'].Draw('colz')
can.cd(3)
hhs['delta_p_vs_theta_km'].Draw('colz')
can.cd(4)
hhs['delta_p_vs_p_km'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()


#close the pdf file
can.Print('{}]'.format(my_pdf_name))


#######################
## plot the fit means
rebinx=8
canfits = TCanvas('canfits','canfits',1000,1200)
canfits.Print('{}['.format(my_pdf_fits))
canfits.Divide(2,1)
canfits.cd(1)
g_epXthetapr, sl_epXthetapr, f_epXthetapr = fit_slices(hhs['epX_vs_theta_pr'], [22,38], rebinx)
g_epXthetapr.SetMarkerStyle(21)
g_epXthetapr.SetMarkerSize(1)
g_epXthetapr.Draw('AP')
g_epXthetapr.GetHistogram().SetMaximum(3.00)
g_epXthetapr.GetHistogram().SetMinimum(0.00)
g_epXthetapr.Draw('AP')
canfits.cd(2)
g_epXmntmpr, sl_epXmntmpr, f_epXmntmpr = fit_slices(hhs['epX_vs_p_pr'], [0.5,1.5], rebinx)
g_epXmntmpr.SetMarkerStyle(21)
g_epXmntmpr.SetMarkerSize(1)
g_epXmntmpr.Draw('AP')
g_epXmntmpr.GetHistogram().SetMaximum(3.00)
g_epXmntmpr.GetHistogram().SetMinimum(0.00)
g_epXmntmpr.Draw('AP')
canfits.Print(my_pdf_fits)
canfits.Clear()

canfits.Divide(2,1)
canfits.cd(1)
g_epkpXthetakm, sl_epkpXthetakm, f_epkpXthetakm = fit_slices(hhs['epkpX_vs_theta_km'], [12,20], rebinx)
g_epkpXthetakm.SetMarkerStyle(21)
g_epkpXthetakm.SetMarkerSize(1)
g_epkpXthetakm.Draw('AP')
g_epkpXthetakm.GetHistogram().SetMaximum(0.40)
g_epkpXthetakm.GetHistogram().SetMinimum(0.00)
g_epkpXthetakm.Draw('AP')
canfits.cd(2)
g_epkpXmntmkm, sl_epkpXmntmkm, f_epkpXmntmkm = fit_slices(hhs['epkpX_vs_p_km'], [0.85,3.0], rebinx)
g_epkpXmntmkm.SetMarkerStyle(21)
g_epkpXmntmkm.SetMarkerSize(1)
g_epkpXmntmkm.Draw('AP')
g_epkpXmntmkm.GetHistogram().SetMaximum(0.40)
g_epkpXmntmkm.GetHistogram().SetMinimum(0.00)
g_epkpXmntmkm.Draw('AP')
canfits.Print(my_pdf_fits)
canfits.Clear()


canfits.Divide(2,1)
canfits.cd(1)
g_epkmXthetakp, sl_epkmXthetakp, f_epkmXthetakp = fit_slices(hhs['epkmX_vs_theta_kp'], [9,20], rebinx)
g_epkmXthetakp.SetMarkerStyle(21)
g_epkmXthetakp.SetMarkerSize(1)
g_epkmXthetakp.Draw('AP')
g_epkmXthetakp.GetHistogram().SetMaximum(0.40)
g_epkmXthetakp.GetHistogram().SetMinimum(0.00)
g_epkmXthetakp.Draw('AP')
canfits.cd(2)
g_epkmXmntmkp, sl_epkmXmntmkp, f_epkmXmntmkp = fit_slices(hhs['epkmX_vs_p_kp'], [0.85,3.0], rebinx)
g_epkmXmntmkp.SetMarkerStyle(21)
g_epkmXmntmkp.SetMarkerSize(1)
g_epkmXmntmkp.Draw('AP')
g_epkmXmntmkp.GetHistogram().SetMaximum(0.40)
g_epkmXmntmkp.GetHistogram().SetMinimum(0.00)
g_epkmXmntmkp.Draw('AP')
canfits.Print(my_pdf_fits)
canfits.Clear()

canfits.Divide(2,1)
canfits.cd(1)
g_ekpkmXthetapr, sl_epkmXthetakp, f_epkmXthetakp = fit_slices(hhs['ekpkmX_vs_theta_pr'], [21,38], rebinx)
g_ekpkmXthetapr.SetMarkerStyle(21)
g_ekpkmXthetapr.SetMarkerSize(1)
g_ekpkmXthetapr.Draw('AP')
g_ekpkmXthetapr.GetHistogram().SetMaximum(1.10)
g_ekpkmXthetapr.GetHistogram().SetMinimum(0.70)
g_ekpkmXthetapr.Draw('AP')
canfits.cd(2)
g_ekpkmXmntmpr, sl_ekpkmXmntmpr, f_ekpkmXmntmpt = fit_slices(hhs['ekpkmX_vs_p_pr'], [0.44,1.3], rebinx)
g_ekpkmXmntmpr.SetMarkerStyle(21)
g_ekpkmXmntmpr.SetMarkerSize(1)
g_ekpkmXmntmpr.Draw('AP')
g_ekpkmXmntmpr.GetHistogram().SetMaximum(1.10)
g_ekpkmXmntmpr.GetHistogram().SetMinimum(0.70)
g_ekpkmXmntmpr.Draw('AP')
canfits.Print(my_pdf_fits)
canfits.Clear()

delta_theta_vs_theta_plots = [hhs['delta_theta_vs_theta_pr'], 
                              hhs['delta_theta_vs_theta_kp'],
                              hhs['delta_theta_vs_theta_km']]

delta_theta_vs_p_plots = [hhs['delta_theta_vs_p_pr'], 
                          hhs['delta_theta_vs_p_kp'],
                          hhs['delta_theta_vs_p_km']]

delta_p_vs_theta_plots = [hhs['delta_p_vs_theta_pr'], 
                          hhs['delta_p_vs_theta_kp'],
                          hhs['delta_p_vs_theta_km']]

delta_p_vs_p_plots = [hhs['delta_p_vs_p_pr'], 
                      hhs['delta_p_vs_p_kp'],
                      hhs['delta_p_vs_p_km']]

delta_theta_vs_theta_ranges = [ [20,38], [9,22], [9,22]]
delta_theta_vs_p_ranges = [ [0.5,1.75],[1.02,2.75], [1.02,2.75]]
delta_p_vs_theta_ranges = [ [20,38], [9,22],[9,22] ]
delta_p_vs_p_ranges = [ [0.5,1.75], [0.8, 2.75], [0.8, 2.75]]

canfits.SetCanvasSize(900,300);
make_plots(canfits, my_pdf_fits, delta_theta_vs_theta_plots, delta_theta_vs_theta_ranges ,rebinx,-3,3)

make_plots(canfits, my_pdf_fits, delta_theta_vs_p_plots, delta_theta_vs_p_ranges ,rebinx,-3.0, 3.0)

make_plots(canfits, my_pdf_fits, delta_p_vs_theta_plots, delta_p_vs_theta_ranges ,rebinx,-1.3, 1.3)

make_plots(canfits, my_pdf_fits, delta_p_vs_p_plots, delta_p_vs_p_ranges ,rebinx, -1.3, 1.3)

can.Print('{}]'.format(my_pdf_fits))

    
canslices = TCanvas('canslices','canslices',900,900)
canslices.Print('{}['.format(my_pdf_slices))

hh_counter=0
for hh in delta_theta_vs_theta_plots:
    make_slice_plots(canslices, my_pdf_slices, hh, delta_theta_vs_theta_ranges[hh_counter], 3, 3, rebinx) 
    hh_counter+=1

hh_counter=0
for hh in delta_theta_vs_p_plots:
    make_slice_plots(canslices, my_pdf_slices, hh, delta_theta_vs_p_ranges[hh_counter], 3, 3, rebinx) 
    hh_counter+=1

hh_counter=0
for hh in delta_p_vs_theta_plots:
    make_slice_plots(canslices, my_pdf_slices, hh, delta_p_vs_theta_ranges[hh_counter], 3, 3, rebinx) 
    hh_counter+=1

hh_counter=0
for hh in delta_p_vs_p_plots:
    make_slice_plots(canslices, my_pdf_slices, hh, delta_p_vs_p_ranges[hh_counter], 3, 3, rebinx) 
    hh_counter+=1


canslices.Print('{}]'.format(my_pdf_slices))

    
