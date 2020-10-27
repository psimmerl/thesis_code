from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TLine
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
    #histo.RebinY(4)
    
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

def rebinXY(hist,rebinfactorX, rebinfactorY):
    hist.RebinY(4)
    hist.RebinX(4)
    

fin = TFile(sys.argv[1])
tag=sys.argv[2]
signal_sigma = 5.0
rebinx=4


hhs = {}
for kk in fin.GetListOfKeys():
    obj = kk.ReadObj()
    
    hhs[obj.GetName()] = obj

my_pdf_name='mc_simres_phi_'+tag+'.pdf'
my_pdf_fits='mc_simres_fits_phi_'+tag+'.pdf'

can =TCanvas('can','can',1200,600)
rebin_hist=True

zLine1a = TLine(0.0, 0.0, hhs['mc_res_el_p_FD'].GetXaxis().GetXmax(), 0.0)
zLine1a.SetLineColor(kRed)
zLine1a.SetLineWidth(2)
zLine1b = TLine(0.0, 0.0, hhs['mc_res_el_theta_FD'].GetXaxis().GetXmax(), 0.0)
zLine1b.SetLineColor(kRed)
zLine1b.SetLineWidth(2)
zLine2a = TLine(0.0, 0.0, hhs['mc_res_pr_p_FD'].GetXaxis().GetXmax(), 0.0)
zLine2a.SetLineColor(kRed)
zLine2a.SetLineWidth(2)
zLine2b = TLine(0.0, 0.0, hhs['mc_res_pr_theta_FD'].GetXaxis().GetXmax(), 0.0)
zLine2b.SetLineColor(kRed)
zLine2b.SetLineWidth(2)


#open pdf file 
can.Print('{}['.format(my_pdf_name))
can.Divide(2,1)
can.cd(1)
if rebin_hist:
    rebinXY(hhs['mc_res_el_p_FD'],4,4)
hhs['mc_res_el_p_FD'].SetTitle('REC - GEN Electron mntm; p (GeV); #Delta p')
hhs['mc_res_el_p_FD'].Draw('colz')
zLine1a.Draw('same')
can.cd(2)
if rebin_hist:
    rebinXY(hhs['mc_res_el_theta_FD'],4,4)
hhs['mc_res_el_theta_FD'].SetTitle('REC - GEN Electron #theta; #theta (deg); #Delta #theta')
hhs['mc_res_el_theta_FD'].Draw('colz')
zLine1b.Draw('same')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
if rebin_hist:
    rebinXY(hhs['mc_res_pr_p_FD'],4,4)
hhs['mc_res_pr_p_FD'].SetTitle('REC - GEN Proton mntm; p (GeV); #Delta p')
hhs['mc_res_pr_p_FD'].Draw('colz')
zLine2a.Draw('same')
can.cd(2)
if rebin_hist:
    rebinXY(hhs['mc_res_pr_theta_FD'],4,4)
hhs['mc_res_pr_theta_FD'].SetTitle('REC - GEN Proton #theta; #theta (deg); #Delta #theta')
hhs['mc_res_pr_theta_FD'].Draw('colz')
zLine2b.Draw('same')
can.Print(my_pdf_name)
can.Clear()

can.SetCanvasSize(1600,800)
can.Divide(4,2)
low_theta=10
for ii in range(1,8):
    max_theta = low_theta + ii*2.5
    min_theta = low_theta + (ii-1)*2.5
    can.cd(ii)
    hhs['mc_res_pr_p_FD_mintheta'+str(min_theta)+'_maxtheta'+str(max_theta)].SetTitle('Proton #Delta P, #theta_{rec} ('+str(min_theta)+','+str(max_theta)+')' ) 
    hhs['mc_res_pr_p_FD_mintheta'+str(min_theta)+'_maxtheta'+str(max_theta)].Draw('colz')
can.Print(my_pdf_name)

can.Clear()
can.SetCanvasSize(1200,600)

can.Divide(2,1)
can.cd(1)
if rebin_hist:
    rebinXY(hhs['mc_res_kp_p_FD'],4,4)
hhs['mc_res_kp_p_FD'].SetTitle('REC - GEN K^{+} mntm; p (GeV); #Delta p')
hhs['mc_res_kp_p_FD'].Draw('colz')
zLine2a.Draw('same')
can.cd(2)
if rebin_hist:
    rebinXY(hhs['mc_res_kp_theta_FD'],4,4)
hhs['mc_res_kp_theta_FD'].SetTitle('REC - GEN K^{+} #theta; #theta (deg); #Delta #theta')
hhs['mc_res_kp_theta_FD'].Draw('colz')
zLine2b.Draw('same')
can.Print(my_pdf_name)
can.Clear()

can.Divide(2,1)
can.cd(1)
can.cd(1)
if rebin_hist:
    rebinXY(hhs['mc_res_km_p_FD'],4,4)
hhs['mc_res_km_p_FD'].SetTitle('REC - GEN K^{-} mntm; p (GeV); #Delta p')
hhs['mc_res_km_p_FD'].Draw('colz')
zLine2a.Draw('same')
can.cd(2)
if rebin_hist:
    rebinXY(hhs['mc_res_km_theta_FD'],4,4)
hhs['mc_res_km_theta_FD'].SetTitle('REC - GEN K^{-} #theta; #theta (deg); #Delta #theta')
hhs['mc_res_km_theta_FD'].Draw('colz')
zLine2b.Draw('same')
can.Print(my_pdf_name)
can.Clear()


can.SetCanvasSize(1800,1200)
can.Divide(3,2)
can.cd(1)
hhs['xb_el_theta_prior_cuts'].Draw('colz')
can.cd(2)
hhs['xb_el_theta_prior_cuts_FD'].Draw('colz')
can.cd(3)
hhs['xb_el_theta_epkpkmXcut_FD'].Draw('colz')
can.cd(4)
hhs['xb_el_theta_epkmXcut_FD'].Draw('colz')
can.cd(5)
hhs['xb_el_theta_epkpXcut_FD'].Draw('colz')
can.cd(6)
#hhs['xb_el_theta_ekpkmXcut_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()


can.SetCanvasSize(1800,1200)
can.Divide(3,2)
can.cd(1)
hhs['xb_pr_theta_prior_cuts'].Draw('colz')
can.cd(2)
hhs['xb_pr_theta_prior_cuts_FD'].Draw('colz')
can.cd(3)
hhs['xb_pr_theta_epkpkmXcut_FD'].Draw('colz')
can.cd(4)
hhs['xb_pr_theta_epkmXcut_FD'].Draw('colz')
can.cd(5)
hhs['xb_pr_theta_epkpXcut_FD'].Draw('colz')
can.cd(6)
#hhs['xb_pr_theta_ekpkmXcut_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.SetCanvasSize(1800,1200)
can.Divide(3,2)
can.cd(1)
hhs['xb_kp_theta_prior_cuts'].Draw('colz')
can.cd(2)
hhs['xb_kp_theta_prior_cuts_FD'].Draw('colz')
can.cd(3)
hhs['xb_kp_theta_epkpkmXcut_FD'].Draw('colz')
can.cd(4)
hhs['xb_kp_theta_epkmXcut_FD'].Draw('colz')
can.cd(5)
hhs['xb_kp_theta_epkpXcut_FD'].Draw('colz')
can.cd(6)
#hhs['xb_kp_theta_ekpkmXcut_FD'].Draw('colz')
can.Print(my_pdf_name)
can.Clear()

can.SetCanvasSize(1800,1200)
can.Divide(3,2)
can.cd(1)
hhs['xb_km_theta_prior_cuts'].Draw('colz')
can.cd(2)
hhs['xb_km_theta_prior_cuts_FD'].Draw('colz')
can.cd(3)
hhs['xb_km_theta_epkpkmXcut_FD'].Draw('colz')
can.cd(4)
hhs['xb_km_theta_epkmXcut_FD'].Draw('colz')
can.cd(5)
hhs['xb_km_theta_epkpXcut_FD'].Draw('colz')
can.cd(6)
#hhs['xb_km_theta_ekpkmXcut_FD'].Draw('colz')
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

el_y_rang = [[-0.1, 0.1], [-0.5, 0.5]]
pr_y_rang = [[-0.15, 0.15], [-2.0, 2.0]]
kp_y_rang = [[-0.05, 0.05], [-0.5, 0.5]]
km_y_rang = [[-0.05, 0.05], [-0.5, 0.5]]

can.SetCanvasSize(1000,500);
make_plots(can, my_pdf_name, mc_res_el_fd, el_kin_ranges ,rebinx,el_y_rang)

make_plots(can, my_pdf_name, mc_res_pr_fd, pr_kin_ranges ,rebinx,pr_y_rang)

make_plots(can, my_pdf_name, mc_res_kp_fd, kp_kin_ranges ,rebinx,kp_y_rang)

make_plots(can, my_pdf_name, mc_res_km_fd, km_kin_ranges ,rebinx, km_y_rang)




#close the pdf file
can.Print('{}]'.format(my_pdf_name))

