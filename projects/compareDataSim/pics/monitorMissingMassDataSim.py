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


def make_plots(can, fname, hists, x_ranges, rebinx, y_range):
    can.Clear()
    can.Divide(2,1)
    cc=1
    root_is_dumb = []
    
    for hh in hists:
        can.cd(cc)
        g_fit, sl_fit, f_fit  = fit_slices(hh,x_ranges[cc-1],rebinx)
        g_fit.SetTitle(hh.GetTitle()+'_slices')
        g_fit.SetMarkerStyle(21)
        g_fit.SetMarkerSize(1)
        g_fit.Draw('AP')
        g_fit.GetHistogram().SetMaximum(y_range[cc-1][1])
        g_fit.GetHistogram().SetMinimum(y_range[cc-1][0])
        g_fit.Draw('AP')
        root_is_dumb.append(g_fit)
        cc+=1
    can.Print(fname)


def define2DHist(htemp, nameX, nameY):
    htemp.SetTitleX(nameX)
    htemp.SetTitleY(nameY)
    

fin = TFile(sys.argv[1])
find = TFile(sys.argv[2])
signal_sigma = 5.0
rebinx=4


hhs = {}
for kk in fin.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj

hhd = {}
for kk in find.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhd[obj.GetName()] = obj

my_pdf_name='monitor_data_sim_missing_mass_inb.pdf'

can =TCanvas('can','can',1800,900)

#open pdf file 
can.Print('{}['.format(my_pdf_name))
can.Divide(2,2)
can.cd(1)
hd_epkpkmX=hhs['epkpkmX_pass_c123']
hs_epkpkmX=hhd['epkpkmX_pass_c123']
hd_epkpkmX.SetLineColor(kRed)
hs_epkpkmX.SetLineColor(kBlue)
hs_epkpkmX.Draw()
hd_epkpkmX.Draw('same')


can.cd(2)
hd_ekpkmX=hhs['ekpkmX_pass_c013']
hs_ekpkmX=hhd['ekpkmX_pass_c013']
hd_ekpkmX.SetLineColor(kRed)
hs_ekpkmX.SetLineColor(kBlue)
hs_ekpkmX.Draw()
hd_ekpkmX.Draw('same')


can.cd(3)
hd_epkpX=hhs['epkpX_pass_c023']
hs_epkpX=hhd['epkpX_pass_c023']
hd_epkpX.SetLineColor(kRed)
hs_epkpX.SetLineColor(kBlue)
hs_epkpX.Draw()
hd_epkpX.Draw('same')

can.cd(4)
hd_epkmX=hhs['epkmX_pass_c012']
hs_epkmX=hhd['epkmX_pass_c012']
hd_epkmX.SetLineColor(kRed)
hs_epkmX.SetLineColor(kBlue)
hs_epkmX.Draw()
hd_epkmX.Draw('same')



can.Print('{}]'.format(my_pdf_name))
















'''

hhs['mc_res_el_p_FD'].Draw('colz')
can.cd(2)
hhs['mc_res_el_theta_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
hhs['mc_res_pr_p_FD'].Draw('colz')
can.cd(2)
hhs['mc_res_pr_theta_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
hhs['mc_res_kp_p_FD'].Draw('colz')
can.cd(2)
hhs['mc_res_kp_theta_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
hhs['mc_res_km_p_FD'].Draw('colz')
can.cd(2)
hhs['mc_res_km_theta_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

mc_res_el_fd = [hhs['mc_res_el_p_FD'],
                hhs['mc_res_el_theta_FD']]

mc_res_pr_fd = [hhs['mc_res_pr_p_FD'],
                hhs['mc_res_pr_theta_FD']]

mc_res_kp_fd = [hhs['mc_res_kp_p_FD'],
                hhs['mc_res_kp_theta_FD']]

mc_res_km_fd = [hhs['mc_res_km_p_FD'],
                hhs['mc_res_km_theta_FD']]

el_kin_ranges=[ [4.5, 7.8], [7.5, 23] ]
pr_kin_ranges=[ [0.5, 1.6], [20,37] ]
kp_kin_ranges=[ [0.96, 3.0], [7.1, 20.0] ]
km_kin_ranges=[ [0.96, 3.2], [10.5, 22.5] ]

el_y_rang = [[-0.15, 0.15], [-2.0, 2.0]]
pr_y_rang = [[-0.15, 0.15], [-2.0, 2.0]]
kp_y_rang = [[-0.15, 0.15], [-2.0, 2.0]]
km_y_rang = [[-0.15, 0.15], [-2.0, 2.0]]

can.SetCanvasSize(1000,500);
make_plots(can, my_pdf_name, mc_res_el_fd, el_kin_ranges ,rebinx,el_y_rang)

make_plots(can, my_pdf_name, mc_res_pr_fd, pr_kin_ranges ,rebinx,pr_y_rang)

make_plots(can, my_pdf_name, mc_res_kp_fd, kp_kin_ranges ,rebinx,kp_y_rang)

make_plots(can, my_pdf_name, mc_res_km_fd, km_kin_ranges ,rebinx, km_y_rang)


#close the pdf file
can.Print('{}]'.format(my_pdf_name))

'''
