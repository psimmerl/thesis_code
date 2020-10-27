from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TH2F, TFile, TMath, TLegend, TLatex
from ROOT import TGraphErrors
from ROOT import gROOT, gBenchmark, gStyle, gPad
from ROOT import TLorentzVector, TVector3
from ROOT import kRed, kBlue, kGreen
import numpy as np

from array import array
import math
import sys,os

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)


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

def make_plots_overlap(can, canx, cany, fname, hists, hh_titles, cut_before, set_axis_range=[] ):
    can.Clear()
    can.Divide(canx,cany)
    cc=1
    root_is_dumb = []
    for hh in hists:
        can.cd(cc)
        cut_before.SetTitle('')#hh_titles[cc-1])
        if cc!=1:
            cut_before.Draw()
            root_is_dumb.append(cut_before)
        hh.Draw('colz+same')
        gPad.SetLogz()

        root_is_dumb.append(hh)
        cc+=1
    can.Print(fname)


def make_plots_per_sector(can, fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[]):
    can.Clear()
    can.Divide(2,3)
    cc=1
    root_is_dumb = []
    ss=1
    for hh in hists:
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.Draw('colz')
        gPad.SetLogz()

        #print('draw cut %d' %(draw_cut ))
        #print('cc mod 3 %d' %( cc%2 ))
        if draw_cut :
            can.Update()
            pad=can.GetPad(cc)
            ymin=0
            ymax=gPad.GetUymax()
            #print(ymax)
            if len(cut_pos[ss-1]) >= 2:
               for cp in cut_pos[ss-1]:               
                   if( min_y > 0 and max_y > 0 ):
                       ymin=min_y              
                       ymax=max_y                
                   print(' sector {} : cut position at {} '.format(ss,cp))
                   cutline = TLine(cp, ymin, cp,ymax)
                   cutline.SetLineColor(kRed)
                   cutline.Draw('same')
                   root_is_dumb.append(cutline)
            else:
                cutline = TLine(cut_pos[ss-1][0], 0.0, cut_pos[ss-1][0],ymax)
                print(' sector {}:  cut position at {} '.format(ss,cut_pos[ss-1][0]))

                cutline.SetLineColor(kRed)
                cutline.Draw('same')
                root_is_dumb.append(cutline)

        root_is_dumb.append(hh)

        if cc%3 == 0:             
            ss+=1  
            if ss == 7:
                ss=0
        cc+=1


    can.Print(fname)


def make_plots_per_sector_withfit(can, fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[], fits =[], fits_to_draw=[]):
    can.Clear()
    can.Divide(2,3)
    cc=1
    root_is_dumb = []
    ss=1
    for hh in hists:
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.RebinX(rebinx)
        hh.Draw('colz')
        # fits has fit name and TF1
        if( len(fits) > 0 ):
            hh.Fit(fits[ss-1][0],'R+')
            
        gPad.SetLogz()

        if( len(fits_to_draw) > 0 ):
            for ii in range(0,len(fits_to_draw[ss-1])):
                fits_to_draw[ss-1][ii].Draw("same")

        #print('draw cut %d' %(draw_cut ))
        #print('cc mod 3 %d' %( cc%2 ))
        if draw_cut :
            can.Update()
            pad=can.GetPad(cc)
            ymin=0
            ymax=gPad.GetUymax()
            #print(ymax)
            if len(cut_pos[ss-1]) >= 2:
               for cp in cut_pos[ss-1]:               
                   if( min_y > 0 and max_y > 0 ):
                       ymin=min_y              
                       ymax=max_y                
                   print(' sector {} : cut position at {} '.format(ss,cp))
                   cutline = TLine(cp, ymin, cp,ymax)
                   cutline.SetLineColor(kRed)
                   cutline.Draw('same')
                   root_is_dumb.append(cutline)
            else:
                cutline = TLine(cut_pos[ss-1][0], 0.0, cut_pos[ss-1][0],ymax)
                print(' sector {}:  cut position at {} '.format(ss,cut_pos[ss-1][0]))

                cutline.SetLineColor(kRed)
                cutline.Draw('same')
                root_is_dumb.append(cutline)

        root_is_dumb.append(hh)

        if cc%3 == 0:             
            ss+=1  
            if ss == 7:
                ss=0
        cc+=1


    can.Print(fname)

def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False,
                     y_fit_range=None, landscape=False, x_range=None, y_range=None, vline=None,
                     hline=None, rebinX=None, scale=None):

    root_garbage_can = []
    
    canvas.Clear() 
    if landscape:
        canvas.Divide(3,2)
    else:
        canvas.Divide(2,3)

    for i in range(0,7):
        canvas.cd(i+1)
        
        #if isinstance(histos[title_formatter.format(i)], TH1F):
        if isinstance(histos.get(title_formatter.format(i), default_histo), TH1F):
            if y_fit_range:
                fit = TF1(title_formatter.format(i) + '_fit', 'gaus')
                histos.get(title_formatter.format(i), default_histo).Fit(fit, '', 'R', y_fit_range[0], y_fit_range[1])
                
            if x_range:
                histos.get(title_formatter.format(i), default_histo).GetXaxis().SetRangeUser(x_range[0], x_range[1])
                                
            histos.get(title_formatter.format(i), default_histo).SetFillColorAlpha(55, 0.65)
            histos.get(title_formatter.format(i), default_histo).Draw()

            if y_fit_range:
                label.DrawLatex(0.15, 0.86, '#mu = {0:6.4f}, #sigma = {1:6.4f}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            
            if rebinX:
                histos.get(title_formatter.format(i), default_histo).RebinX(rebinX)

            if scale:
                histos.get(title_formatter.format(i), default_histo).Scale(scale)
            
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
            label.DrawLatex(0.1, 0.925, title + ' S' + str(i))

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

    


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

h_el_ptheta = []
h_el_ptheta_name = []
h_el_vzphi = []
h_el_vzphi_name = []
h_el_vztheta = []
h_el_vztheta_name = []
h_el_phiw = []
h_el_phiw_name = []


h_el_w_small =[]
h_el_w_small_names =[]
for ss in range(0,6):
    if 'hw_small_el_pass_all_' + str(ss) in hhs:
        h_el_w_small.append(hhs['hw_small_el_pass_all_' + str(ss)])
        h_el_w_small[ss].RebinX(2)
        h_el_w_small_names.append(" W Sector " + str(ss+1) + "; W (GeV); counts")

    if 'hptheta_el_pass_all'+str(ss) in hhs:
        h_el_ptheta.append(hhs['hptheta_el_pass_all'+str(ss)])
        h_el_ptheta_name.append('El P vs #theta PCAL Sector '+str(ss+1)+'; P (GeV); #theta (deg)')

    if 'hvztheta_el_pass_all'+str(ss) in hhs:
        h_el_vztheta.append(hhs['hvztheta_el_pass_all'+str(ss)])
        h_el_vztheta_name.append('El Vz vs #theta PCAL Sector ' + str(ss+1)+'; Vz (cm); #theta (deg)')
    


h_el_vzphi.append(hhs['hphivz_el_pass_all'])
h_el_vzphi_name.append('El Vz vs #phi; vz (cm); #phi (deg)')

h_el_phiw.append(hhs['hphiw_el_pass_all'])
h_el_phiw_name.append('El #phi vs W (GeV); #phi (deg); W (GeV)')


plot_sector_page(can, hhs,'hptheta_el_pass_all{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) P vs #theta', xtitle='p (GeV)', ytitle='#theta (deg)', log=True)
plot_sector_page(can, hhs,'hvztheta_el_pass_all{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) Vz vs #theta', xtitle='vz (cm)', ytitle='#theta (deg)', log=True)

make_plots(can,1,1,mon_out_file_name,h_el_vzphi,h_el_vzphi_name)
make_plots(can,1,1,mon_out_file_name,h_el_phiw,h_el_phiw_name)

plot_sector_page(can, hhs,'hw_el_pass_all_{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W ', xtitle='W (GeV)', scale=1541.8455)




plot_sector_page(can, hhs,'hw_small_el_pass_all_{}', lab, save_name=mon_out_file_name,
                 title='El. (FD) W - zoom ', xtitle='W (GeV)', rebinX=2)

'''

can.Clear()
can.Divide(2,3)
w_small_fits=[]
temp_f=[]
for ss in range(0,6):
    can.cd(ss+1)
    w_tf1=[]
    w_tf1.append('w_small_fit_sector'+str(ss))
    w_tf1.append(TF1('w_small_fit_sector'+str(ss), 'gaus(0)*pol3(3)', 0.8, 1.25))
    w_tf1[1].SetParameter(0,50)
    w_tf1[1].SetParameter(1,1.04)
    w_tf1[1].SetParameter(2,0.1)    
    #w_tf1[1].SetParameter(3,1)    
    #w_tf1[1].SetParameter(4,1)    
    #w_tf1[1].SetParameter(5,1)    
    
    w_small_fits.append(w_tf1)
    h_el_w_small[ss].Draw('same')
    h_el_w_small[ss].Fit(w_small_fits[ss][0],'R')
    
    temp_tf1_gaus = TF1("temp_gaus_s"+str(ss),"gaus",0.8,1.25)
    temp_tf1_pol3 = TF1("temp_pol3_s"+str(ss),"pol3",0.8,1.25)
    print('sector ' + str(ss))
    print(w_small_fits[ss][1].GetParameter(0))
    print(w_small_fits[ss][1].GetParameter(1))
    print(w_small_fits[ss][1].GetParameter(2))
    print(w_small_fits[ss][1].GetParameter(3))
    print(w_small_fits[ss][1].GetParameter(4))
    print(w_small_fits[ss][1].GetParameter(5))
    temp_tf1_gaus.SetParameter(0,w_small_fits[ss][1].GetParameter(0))
    temp_tf1_gaus.SetParameter(1,w_small_fits[ss][1].GetParameter(1))
    temp_tf1_gaus.SetParameter(2,w_small_fits[ss][1].GetParameter(2))

    temp_tf1_pol3.SetParameter(0,w_small_fits[ss][1].GetParameter(3))
    temp_tf1_pol3.SetParameter(1,w_small_fits[ss][1].GetParameter(4))
    temp_tf1_pol3.SetParameter(2,w_small_fits[ss][1].GetParameter(5))
    temp_tf1_gaus.SetLineColor(kBlue)
    temp_tf1_pol3.SetLineColor(kGreen)
    temp_tf1_gaus.Draw('same')
    temp_tf1_pol3.Draw('same')
    temp_f.append(temp_tf1_gaus)
    temp_f.append(temp_tf1_pol3)

    #estimate area under fitted function to data
    w_tf1[1]
    
can.Print(mon_out_file_name)

#make_plots_per_sector_withfit(can,mon_out_file_name,h_el_w_small,h_el_w_small_names,fits=w_small_fits, rebinx=2)

#set TF1 parameters based on fit results so we can draw them

w_small_fits_to_draw={}
for ss in range(0,6):
    temp_f = []
    temp_tf1_gaus = TF1("temp_gaus_s"+str(ss),"gaus",0.8,1.)
    temp_tf1_pol3 = TF1("temp_pol3_s"+str(ss),"pol3",0.8,1.)
    print('sector ' + str(ss))
    print(w_small_fits[ss][1].GetParameter(0))
    print(w_small_fits[ss][1].GetParameter(1))
    print(w_small_fits[ss][1].GetParameter(2))
    print(w_small_fits[ss][1].GetParameter(3))
    print(w_small_fits[ss][1].GetParameter(4))
    print(w_small_fits[ss][1].GetParameter(5))
    temp_tf1_gaus.SetParameter(0,w_small_fits[ss][1].GetParameter(0))
    temp_tf1_gaus.SetParameter(1,w_small_fits[ss][1].GetParameter(1))
    temp_tf1_gaus.SetParameter(2,w_small_fits[ss][1].GetParameter(2))

    temp_tf1_pol3.SetParameter(3,w_small_fits[ss][1].GetParameter(3))
    temp_tf1_pol3.SetParameter(4,w_small_fits[ss][1].GetParameter(4))
    temp_tf1_pol3.SetParameter(5,w_small_fits[ss][1].GetParameter(5))
    temp_tf1_gaus.SetLineColor(kBlue)
    temp_tf1_pol3.SetLineColor(kGreen)
    temp_f.append(temp_tf1_gaus)
    temp_f.append(temp_tf1_pol3)
    w_small_fits_to_draw[ss]=temp_f
'''
make_plots_per_sector_withfit(can,mon_out_file_name,h_el_w_small,h_el_w_small_names,fits_to_draw=w_small_fits_to_draw)
'''

can.Print('{}]'.format(mon_out_file_name))















    
    
    
    
