from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLine
from ROOT import TGraphErrors, TMultiGraph
from ROOT import gROOT, gBenchmark, gStyle, gPad, TLatex
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen, kSpring, kBlack
from ROOT import kTRUE, gROOT
import numpy as np

from array import array
import math
import sys,os

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)

gROOT.SetBatch(kTRUE);

def getFitLimits(h_temp, sig):
    temp_mean = h_temp.GetMean()
    temp_rms = h_temp.GetRMS()
    temp_min_fit_range = temp_mean - sig*temp_rms
    temp_max_fit_range = temp_mean + sig*temp_rms
    return [ temp_min_fit_range, temp_max_fit_range, temp_mean, temp_rms ]


def getFitFractionalHeight( h_temp, h_name, percent_max, save_fits):
    # fitting routine adopted from Andrew Puckett.
    print(' special fitting routine ')
    print(h_name)
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
    
    fit_temp = TF1(h_name,'gaus(0)', xlow, xhigh )
    fit_temp.SetParameter(0, histmax)
    fit_temp.SetParameter(1, h_temp.GetBinCenter(h_temp.GetMaximumBin()))
    
                          
    
    h_temp.Fit(h_name,"R") # add Q for quiet
    save_fits.append(fit_temp)
    
    temp_mean = fit_temp.GetParameter(1)
    temp_rms = fit_temp.GetParameter(2)
    print('fit results {} {}'.format(temp_mean, temp_rms))
    return fit_temp#[temp_mean, temp_rms ]


def plot_page(canvas, histos, histo_title, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False, vline=None, fit=None, rebinX=None, rebinY=None):
    
    canvas.Clear() 
    line=None
    root_garbage_can = []


    if isinstance(histos[histo_title], TH1F):
         histos[histo_title].Draw()
    elif isinstance(histos[histo_title], TH2F):
        if rebinX:
            histos[histo_title].RebinX(rebinX)
        if rebinY:
            histos[histo_title].RebinY(rebinY)
        histos[histo_title].Draw('colz')
        if log:
            gPad.SetLogz() 
    else:
        #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
        pass
    
    if title:
        histos.get(histo_title, default_histo).SetTitle('')
        label.DrawLatex(0.1, 0.925, title)
        
    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    if vline:
        can.Update()
        ymax = gPad.GetUymax()
        for vv in vline:
            print(ymax)
            print(vv)
        
            line = TLine(vv, 0,
                         vv, ymax )
            
            line.SetLineColor(kRed)
            line.SetLineStyle(1)
            line.SetLineWidth(1)
            line.Draw('same')
            root_garbage_can.append(line)
        
    if fit:        
        fit_result = getFitFractionalHeight(histos[histo_title], histos[histo_title].GetTitle(), 0.4, root_garbage_can)
        fit_result.SetLineColor(kGreen)
        fit_result.Draw('same')
        label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.6f}'.format(fit_result.GetParameter(1), fit_result.GetParameter(2)))

        


    canvas.Print(save_name)

def plot_page_many_overlap( canvas, histos, histo_title, out_come, label, save_name,
                            titles=None, log=False, vline=None, fit=None, rebinX=None, rebinY=None, line_colors=None):
    
    canvas.Clear() 
    line=None
    root_garbage_can = []
    canvas.Divide(2,2)

    cc=0
    
    for hh in histo_title:
        canvas.cd(cc+1)
        ii=0
        print(hh)
        for oo in out_come:
            print(oo)
            if line_colors:
                histos[hh.format(oo)].SetLineColor(line_colors[ii])
                                
            if isinstance(histos[hh.format(oo)], TH1F):
                print(' drawing pic ')
                histos[hh.format(oo)].Draw("same")
                root_garbage_can.append(histos[hh.format(oo)])
            else:
                #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
                pass
            ii+=1
        if titles:
            histos.get(hh.format(oo), default_histo).SetTitle('')
            # title
            label.DrawLatex(0.1, 0.925, titles[cc][0])                   
            # x title
            label.DrawLatex(0.5, 0.015, titles[cc][1])
            # y title
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, titles[cc][2])
            label.SetTextAngle(0)
            root_garbage_can.append(label)
            #plot bin information
            #label.DrawLatex(0.12, 0.85, 'BinWidth:{0:.5f} GeV'.format(histos[hh].GetBinWidth(2)))
            #label.DrawLatex(0.12, 0.8, 'range: [{0:1.2f},{1:1.2f}]'.format(histos[hh].GetXaxis().GetXmin(), histos[hh].GetXaxis().GetXmax()))


        cc+=1
        
    canvas.Print(save_name)


def plot_page_many(canvas, histos, histo_title, label, save_name,
                   titles=None, log=False, vline=None, fit=None, rebinX=None, rebinY=None):
    
    canvas.Clear() 
    line=None
    root_garbage_can = []
    canvas.Divide(2,2)

    cc=0
    for hh in histo_title:
        canvas.cd(cc+1)
        if isinstance(histos[hh], TH1F):
            histos[hh].Draw()
        elif isinstance(histos[hh], TH2F):
            if rebinX:
                histos[hh].RebinX(rebinX)
            if rebinY:
                histos[hh].RebinY(rebinY)
            histos[hh].Draw('colz')
            if log:
                gPad.SetLogz() 
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass
    
        if titles:
            histos.get(hh, default_histo).SetTitle('')
            # title
            label.DrawLatex(0.1, 0.925, titles[cc][0])                   
            # x title
            label.DrawLatex(0.5, 0.015, titles[cc][1])
            # y title
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, titles[cc][2])
            label.SetTextAngle(0)

        # plot bin information
        label.DrawLatex(0.12, 0.85, 'BinWidth:{0:.5f} GeV'.format(histos[hh].GetBinWidth(2)))
        label.DrawLatex(0.12, 0.8, 'range: [{0:1.2f},{1:1.2f}]'.format(histos[hh].GetXaxis().GetXmin(), histos[hh].GetXaxis().GetXmax()))


        if vline:
            can.Update()
            ymax = gPad.GetUymax()
            for vv in vline:
                print(ymax)
                print(vv)
                line = TLine(vv, ymax,
                             vv, ymax )
                #histos.get(histo_title, default_histo).GetYaxis().GetXmax())
                #histos.get(histo_title, default_histo).GetYaxis().GetXmin(),
                line.SetLineColor(kRed)
                line.SetLineStyle(1)
                line.SetLineWidth(1)
                line.Draw('same')
                root_garbage_can.append(line)
        
        if fit:        
            fit_result = getFitFractionalHeight(histos[hh], histos[hh].GetTitle(), 0.4, root_garbage_can)
            fit_result.SetLineColor(kGreen)
            fit_result.Draw('same')
            label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.6f}'.format(fit_result.GetParameter(1), fit_result.GetParameter(2)))
        cc+=1
        


    canvas.Print(save_name)


def plot_page_overlap(canvas, histos, histo_titles, label, save_name,
                      xtitle=None, ytitle=None, title=None, log=False, line_colors=None, fit=None, which_to_fit=None, pass_rate=None):
    
    canvas.Clear() 
    canvas.Divide(1,1)
    canvas.cd(1)
    root_garbage_can=[]
    ii=0
    

    for histo_title in histo_titles:
        if isinstance(histos[str(histo_title)], TH1F):
            if line_colors:
                histos[histo_title].SetLineColor(line_colors[ii])
            else:
                histos[histo_title].SetLineColor(ii+1)
            if ii < 1:
                histos[histo_title].Draw()
            else:
                histos[histo_title].Draw('same')
        elif isinstance(histos[histo_title], TH2F):
            histos[histo_title].Draw('colz')
            if log:
                gPad.SetLogz() 
        else:
            #raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            pass
    
        if title:
            histos.get(histo_title, default_histo).SetTitle('')
            label.DrawLatex(0.1, 0.925, title)
            
        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)
        if fit and which_to_fit == ii:        
            fit_result = getFitFractionalHeight(histos[histo_title], histos[histo_title].GetName(), 0.4, root_garbage_can)
            label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit_result.GetParameter(1), fit_result.GetParameter(2)))


        ii+=1

    label.DrawLatex(0.15, 0.825, "All: {}".format(pass_rate[0]))
    label.DrawLatex(0.15, 0.775, "#color[2]{Pass Me}"+": {}".format(pass_rate[2]))
    label.DrawLatex(0.15, 0.725, "#color[4]{Fail Me}"+": {}".format(pass_rate[1]))
    canvas.Print(save_name)



# script to plot, view, correct electron and proton momentum
# in each theta and phi bin

# load data and file
ap = argparse.ArgumentParser()#Paul addition
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
datatype = args.output_prefix 
ff = TFile(input_rootfile)#sys.argv[1])
#datatype = sys.argv[2]

hhs={}
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj


base=os.path.basename(ff.GetName())
print('Analysing file %s' % (base) )

mon_out_file_name='monitor_phi_'+datatype+'.pdf'
can = TCanvas('can','can')
can.SetBatch(kTRUE)

can.SetCanvasSize(1200,1200)
can.Print('{}['.format(mon_out_file_name))

lab = TLatex()
lab.SetNDC()
lab.SetTextFont(42)
lab.SetTextSize(0.05)
lab.SetTextColor(1)


gStyle.SetOptStat(00000)

part_type=['el','pr','kp','km']
had_type=['pro','kp','km']
mm2_type = ['eX','epX','ekpX','epkpX','epkmX','ekpkmX','epkpkmXe']




'''

for pp in part_type:
    plot_page(can,hhs,pp+'_p_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='P (GeV)', ytitle='counts', title=pp+' P - pass all excl. cuts')

for pp in part_type:
    plot_page(can,hhs,pp+'_theta_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='#Theta (deg)', ytitle='counts', title=pp+' #Theta - pass all excl. cuts')

for pp in part_type:
    plot_page(can,hhs,pp+'_ptheta_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='P (GeV)', ytitle='#Theta (deg)', title=pp+' #Theta vs P - pass all excl. cuts')

for pp in part_type:
    plot_page(can,hhs,pp+'_phitheta_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='P (GeV)', ytitle='#Theta (deg)', title=pp+' #Theta vs P - pass all excl. cuts')

for pp in had_type:
    plot_page(can,hhs,pp+'_beta_vs_p_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='P (GeV)', ytitle='#beta', title=pp+' #beta vs P - pass all excl. cuts')


###
## look at distribution of phi events in FD
n_ev_pr_fd_kp_fd_km_fd = hhs['el_ptheta_FD_FD_FD_pass_all_testing'].GetEntries()
n_ev_pr_fd_kp_fd_km_cd = hhs['el_ptheta_FD_FD_CD_pass_all_testing'].GetEntries()
n_ev_pr_fd_kp_cd_km_cd = hhs['el_ptheta_FD_CD_CD_pass_all_testing'].GetEntries()
n_ev_pr_cd_kp_fd_km_fd = hhs['el_ptheta_CD_FD_FD_pass_all_testing'].GetEntries()
#n_ev_pr_cd_kp_cd_km_fd = hhs['el_ptheta_CD_CD_FD_pass_all_testing'].GetEntries()
n_ev_pr_cd_kp_cd_km_cd = hhs['el_ptheta_CD_CD_CD_pass_all_testing'].GetEntries()
n_ev_pr_fd_kp_cd_km_fd = hhs['el_ptheta_FD_CD_FD_pass_all_testing'].GetEntries()
n_ev_pr_cd_kp_fd_km_cd = hhs['el_ptheta_CD_FD_CD_pass_all_testing'].GetEntries()
print(n_ev_pr_fd_kp_fd_km_fd)
print(n_ev_pr_fd_kp_fd_km_cd)
print(n_ev_pr_fd_kp_cd_km_cd)
print(n_ev_pr_cd_kp_fd_km_fd)
#print(n_ev_pr_cd_kp_cd_km_fd)
print(n_ev_pr_cd_kp_cd_km_cd)
print(n_ev_pr_fd_kp_cd_km_fd)
print(n_ev_pr_cd_kp_fd_km_cd)


for mm2 in mm2_type:
    plot_page(can,hhs,mm2+'_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle=mm2, ytitle='counts', title='ep -> ' + mm2 + ' - pass all excl. cuts')


plot_page(can,hhs,'pt_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='P_{t}', ytitle='counts', title=' P_{t} - pass all excl. cuts')
plot_page(can,hhs,'xq2_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='X_{b}', ytitle='Q^{2} (GeV^{2})', title=' Q^{2} vs X_{b} - pass all excl. cuts')
plot_page(can,hhs,'q2t_pass_all_exclusivity',lab,save_name=mon_out_file_name, xtitle='Q^{2} (GeV^{2})', ytitle='-t (GeV^{2})', title='-t vs Q^{2} - pass all excl. cuts')
plot_page(can,hhs,'txb_pass_all_exclusivity',lab,save_name=mon_out_file_name, ytitle='X_{b}', xtitle='-t (GeV^{2})', title='X_{b} vs -t - pass all excl. cuts')


plot_page(can,hhs,'xq2_pass_all_exclusivity_select_phi',lab,save_name=mon_out_file_name, xtitle='X_{b}', ytitle='Q^{2} (GeV^{2})', title=' Q^{2} vs X_{b} - Final #phi Events')
plot_page(can,hhs,'q2t_pass_all_exclusivity_select_phi',lab,save_name=mon_out_file_name, xtitle='Q^{2} (GeV^{2})', ytitle='-t (GeV^{2})', title='-t vs Q^{2} - Final #phi Events', rebinX=4, rebinY=4)
plot_page(can,hhs,'txb_pass_all_exclusivity_select_phi',lab,save_name=mon_out_file_name, ytitle='X_{b}', xtitle='-t (GeV^{2})', title='X_{b} vs -t - Final #phi Events')


sig=4
epkpkmxe_min = 0.08 - 0.0398*sig
epkpkmxe_max = 0.08 + 0.0398*sig

epkpX_min = 0.248 - 0.059*sig
epkpX_max = 0.248 + 0.059*sig

epkmX_min = 0.248 - 0.055*sig
epkmX_max = 0.248 + 0.055*sig

ekpkmX_min = 0.904 - 0.0719*sig
ekpkmX_max = 0.904 + 0.0719*sig

plot_page(can,hhs,'epkpkmXe_raw',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='M. Energy (GeV)', title='epK^{+}K^{-}X Missing Energy - Raw', vline=[epkpkmxe_min, epkpkmxe_max])
plot_page(can,hhs,'ekpkmX_raw',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='MM2 (GeV^2)', title='eK^{+}K^{-}X Missing Pro. - Raw')
plot_page(can,hhs,'epkpX_raw',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='MM2 (GeV^2)', title='epK^{+}X Missing K^{-} - Raw')
plot_page(can,hhs,'epkmX_raw',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='MM2 (GeV^2)', title='epK^{-}X Missing K_{+} - Raw')


#### now make plots from the makeFinalPlots.py file - still good for when needing a png, but this puts all in one spot
'''

## what I cut on for selecting the events below
cut_limits = [0.0295104760535 - 4.0*0.0676640007682, 0.0295104760535 + 4.0*0.0676640007682]
print( cut_limits)
plot_page(can, hhs, 'missing_energy_pass_all_but_missing_energy', lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Missing Energy (GeV)', title='epK^{+}K^{-}X Missing Energy To Cut On', vline=[cut_limits[0], cut_limits[1]])

mm_plots = ['mm_epkpkmX_mod_combssize1',
            'mm_ekpkmX_mod_combssize1',
            'mm_epkpX_mod_combssize1',
            'mm_epkmX_mod_combssize1']

mm_pass_me_plots = ['mm_epkpkmX_mod_pass_only_missing_energy',
                    'mm_ekpkmX_mod_pass_only_missing_energy',
                    'mm_epkpX_mod_pass_only_missing_energy',
                    'mm_epkmX_mod_pass_only_missing_energy']

mm_fail_me_plots = ['mm_epkpkmX_mod_fail_only_missing_energy',
                    'mm_ekpkmX_mod_fail_only_missing_energy',
                    'mm_epkpX_mod_fail_only_missing_energy',
                    'mm_epkmX_mod_fail_only_missing_energy']

mm_plots_title = [ ['epK^{ +}K^{ -}X Missing Mass^{2}, Before M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                   ['eK^{ +}K^{ -}X Missing Mass^{2}, Before M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                   ['epK^{ +}X Missing Mass^{2}, Before M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                   ['epK^{ -}X Missing Mass^{2}, Before M_{e} Cut','Missing Mass^{2} (GeV^{2})','counts'] ]

mm_pass_me_plots_title = [ ['epK^{ +}K^{ -}X Missing Mass^{2}, Pass M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                           ['eK^{ +}K^{ -}X Missing Mass^{2}, Pass M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                           ['epK^{ +}X Missing Mass^{2}, Pass M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                           ['epK^{ -}X Missing Mass^{2}, Pass M_{e} Cut','Missing Mass^{2} (GeV^{2})','counts'] ]

mm_fail_me_plots_title = [ ['epK^{ +}K^{ -}X Missing Mass^{2}, Fail M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                           ['eK^{ +}K^{ -}X Missing Mass^{2}, Fail M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                           ['epK^{ +}X Missing Mass^{2}, Fail M_{e} Cut', 'Missing Mass^{2} (GeV^{2})','counts'],
                           ['epK^{ -}X Missing Mass^{2}, Fail M_{e} Cut','Missing Mass^{2} (GeV^{2})','counts'] ]




# final sample before cuts on missing energy
plot_page_many(can, hhs, mm_plots, lab, save_name=mon_out_file_name, titles = mm_plots_title)

# events passing me cuts
plot_page_many(can, hhs, mm_pass_me_plots, lab, save_name=mon_out_file_name, titles = mm_pass_me_plots_title)

# events failing me cuts
plot_page_many(can, hhs, mm_fail_me_plots, lab, save_name=mon_out_file_name, titles = mm_fail_me_plots_title)

    

all_pass_fail_mm = ['mm_epkpkmX_mod_{}','mm_ekpkmX_mod_{}','mm_epkpX_mod_{}','mm_epkmX_mod_{}']
out_comes=['combssize1', 'pass_only_missing_energy', 'fail_only_missing_energy']
all_pass_fail_titles = [ ['epK^{ +}K^{ -}X Missing Mass^{2}', 'Missing Mass^{2} (GeV^{2})','counts'], 
                         [ 'eK^{ +}K^{ -}X Missing Mass^{2}', 'Missing Mass^{2} (GeV^{2})','counts'],
                         ['epK^{ +}X Missing Mass^{2}', 'Missing Mass^{2} (GeV^{2})','counts' ],
                         ['epK^{ -}X Missing Mass^{2}','Missing Mass^{2} (GeV^{2})','counts'] ]

plot_page_many_overlap( can, hhs, all_pass_fail_mm, out_comes, lab, save_name=mon_out_file_name,
                        titles=all_pass_fail_titles, line_colors=[kBlack, kRed, kBlue])
    

#### need the below ones - but not critical in their current format, will keep nonetheless
absoulte_pass_fail_values = [239555, 176333, 63222]
plot_page_overlap(can,hhs,['mm_epkpkmX_mod_combssize1','mm_epkpkmX_mod_fail_only_missing_energy','mm_epkpkmX_mod_pass_only_missing_energy',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Missing Mass^{2} (GeV^{2})', title='epK^{ +}K^{ -}X Missing Mass^{2}', line_colors=[kBlack,kBlue,kRed], pass_rate = absoulte_pass_fail_values)#, fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['mm_ekpkmX_mod_combssize1','mm_ekpkmX_mod_fail_only_missing_energy','mm_ekpkmX_mod_pass_only_missing_energy',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Missing Mass^{2} (GeV^{2})', title='eK^{ +}K^{ -}X Missing Mass^{2}', line_colors=[kBlack,kBlue,kRed], pass_rate = absoulte_pass_fail_values)#, fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['mm_epkpX_mod_combssize1','mm_epkpX_mod_fail_only_missing_energy','mm_epkpX_mod_pass_only_missing_energy',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Missing Mass^{2} (GeV^{2})', title='epK^{ +}X Missing Mass^{2}', line_colors=[kBlack,kBlue,kRed], pass_rate = absoulte_pass_fail_values)#, fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['mm_epkmX_mod_combssize1','mm_epkmX_mod_fail_only_missing_energy','mm_epkmX_mod_pass_only_missing_energy',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Missing Mass^{2} (GeV^{2})', title='epK^{ -}X Missing Mass^{2}', line_colors=[kBlack,kBlue,kRed], pass_rate = absoulte_pass_fail_values)#, fit=True, which_to_fit=2)




'''
plot_page_overlap(can,hhs,['epkpkmXe_raw','epkpkmXe_fail_c123','epkpkmXe_pass_c123',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Missing Energy (GeV)', title='epK^{+}K^{-}X Missing Energy', line_colors=[kBlack,kBlue,kRed], fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['ekpkmX_raw','ekpkmX_fail_c013','ekpkmX_pass_c013',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='MM2 (GeV^2)', title='eK^{+}K^{-}X Missing Pro. ', line_colors=[kBlack,kBlue,kRed],fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['epkpX_raw','epkpX_fail_c023','epkpX_pass_c023',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='MM2 (GeV^2)', title='epK^{+}X Missing K^{-}', line_colors=[kBlack,kBlue,kRed], fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['epkmX_raw','epkmX_fail_c012','epkmX_pass_c012',],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='MM2 (GeV)', title='epK^{-}X Missing K^{+}', line_colors=[kBlack,kBlue,kRed], fit=True, which_to_fit=2)

plot_page_overlap(can,hhs,['cplpro_raw','cplpro_pass_all_exclusivity'],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='#Theta (deg)', title='Coplanarity Angle Pr', line_colors=[kBlack,kRed])
plot_page_overlap(can,hhs,['cplkp_raw','cplpro_pass_all_exclusivity'],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='#Theta (deg)', title='Coplanarity Angle Pr', line_colors=[kBlack,kRed])
plot_page_overlap(can,hhs,['cplkm_raw','cplpro_pass_all_exclusivity'],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='#Theta (deg)', title='Coplanarity Angle Pr', line_colors=[kBlack,kRed])

plot_page_overlap(can,hhs,['vphimass_pos_helicity_pass_all_exclusivity','vphimass_neg_helicity_pass_all_exclusivity'],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Inv. Mass K^{+}K^{-} (GeV)', title='Inv. Mass K^{+}K^{-} - Pass All Excl. Cuts', line_colors=[kBlack,kRed])

plot_page(can,hhs,'vphimass_pass_all_exclusivity',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Inv. Mass K^{+}K^{-} (GeV)', title='Inv. Mass K^{+}K^{-} - Pass All Excl. Cuts')


#pass_all_exclusivity_select_phi_rmv_lambda
#pass_all_exclusivity_phi_rmv_lambda
#pass_all_exclusivity_rmv_both_lambdas
#pass_all_exclusivity_select_phi_rmv_lambda_pro_mntm_cut
plot_page_overlap(can,hhs,['lambdamass_raw','lambdamass_pass_all_exclusivity'],lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Inv. Mass Pr. K^{-} (GeV)', title='Inv. Mass Pr K^{-} - Pass All Excl. Cuts', line_colors=[kBlack,kRed])

plot_page(can,hhs,'lambdaphimass_pass_all_exclusivity',lab,save_name=mon_out_file_name, ytitle='K^{+}K^{-} (GeV)',xtitle='Mass Pr. K^{-} (GeV)', title='InvMass Pr. K^{-} (GeV) vs K^{+}K^{-} (GeV)' )

plot_page(can,hhs,'vphimass_pass_all_exclusivity',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Inv. Mass K^{+}K^{-} (GeV)', title='Inv. Mass K^{+}K^{-} - Pass All Excl. Cuts')
plot_page(can,hhs,'vphimass_pass_all_exclusivity_phi_rmv_lambda',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Inv. Mass K^{+}K^{-} (GeV)', title='Inv. Mass K^{+}K^{-} - Remove 1520')
plot_page(can,hhs,'vphimass_pass_all_exclusivity_rmv_both_lambdas',lab,save_name=mon_out_file_name, ytitle='counts', xtitle='Inv. Mass K^{+}K^{-} (GeV)', title='Inv. Mass K^{+}K^{-} - Remove both #Lambda s')

for pp in had_type:
    plot_page(can,hhs,pp+'_beta_vs_p_pass_all_exclusivity_remove_lambda',lab,save_name=mon_out_file_name, xtitle='P (GeV)', ytitle='#beta', title=pp+' #beta vs P - pass all excl. cuts')


plot_page(can,hhs,'pos_asy_pass_all_exclusivity_remove_lmbd_vphim',lab,save_name=mon_out_file_name, ytitle='Counts', xtitle='#Phi_{trento} (deg)', title='Select #Phi, Remove #Lambda s #Phi_{trento} (deg) POS. Hel.')

plot_page(can,hhs,'neg_asy_pass_all_exclusivity_remove_lmda_vphim',lab,save_name=mon_out_file_name, ytitle='Counts', xtitle='#Phi_{trento} (deg)', title='Select #Phi, Remove #Lambda s #Phi_{trento} (deg) NEG. Hel.')

'''



can.Print('{}]'.format(mon_out_file_name))
