from ROOT import TFile, TH1F, TH2F, TF1, TLine, TMultiGraph, TGraph, TGraphErrors
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kSpring, kBlue, kAzure
from array import array

import sys,os
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
    can.Divide(3,6)
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
gStyle.SetPalette(55)

beamE = 10.604
hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj

mon_out_file_name='monitor_elpion.pdf'
can = TCanvas('can','can',4000,4000)
can.Print('{}['.format(mon_out_file_name))

elpim_ptheta_pass = [hhs['el_ptheta_passall'],hhs['pim_ptheta_passall'],hhs['el_ptheta_fail_at_least1'],hhs['pim_ptheta_fail_at_least1']]
elpim_ptheta_pass_names = ['El Pass EB #theta vs p ; p (GeV); #theta (deg)','PIM Pass EB #theta vs p ; p (GeV); #theta (deg)', 'El Fail EB #theta vs p ; p (GeV); #theta (deg)', ' PIM Fail EB #theta vs p ; p (GeV); #theta (deg)'] 
make_plots(can,2,2,mon_out_file_name,elpim_ptheta_pass,elpim_ptheta_pass_names)

elpim_phitheta_pass = [hhs['el_phitheta_passall'],hhs['pim_phitheta_passall'],hhs['el_phitheta_fail_at_least1'],hhs['pim_phitheta_fail_at_least1']]
elpim_phitheta_pass_names = ['El Pass EB #phi vs #theta  ;  #theta (deg); #phi (deg)','PIM Pass EB #phi vs #theta  ;  #theta (deg); #phi (deg)', 'El Fail EB #phi vs #theta  ;  #theta (deg); #phi (deg)', ' PIM Fail EB #phi vs #theta  ;  #theta (deg); #phi (deg)'] 
make_plots(can,2,2,mon_out_file_name,elpim_phitheta_pass,elpim_phitheta_pass_names)

elpim_vzphi_pass = [hhs['el_vzphi_passall'],hhs['pim_vzphi_passall'],hhs['el_vzphi_fail_at_least1'],hhs['pim_vzphi_fail_at_least1']]
elpim_vzphi_pass_names = ['El Pass EB Vz vs #phi ;  vz (cm) ; #phi (deg)','PIM Pass EB Vz vs #phi ;  vz (cm) ; #phi (deg)', 'El Fail EB Vz vs #phi ;  vz (cm) ; #phi (deg)', ' PIM Fail EB Vz vs #phi ;  vz (cm) ; #phi (deg)']
make_plots(can,2,2,mon_out_file_name,elpim_vzphi_pass,elpim_vzphi_pass_names)

elpim_resp_pass = [hhs['el_res_p'],hhs['pim_res_p'],hhs['el_res_p_fail_at_least1'],hhs['pim_res_p_fail_at_least1']]
elpim_resp_pass_names = ['El Pass EB \Delta p vs GEN P;  p GEN (GeV) ; \Delta p','PIM Pass EB \Delta p vs GEN P;  p GEN (GeV) ; \Delta p', 'El Fail \Delta p vs GEN P;  p GEN (GeV) ; \Delta p', ' PIM Fail EB \Delta p vs GEN P;  p GEN (GeV) ; \Delta p']
make_plots(can,2,2,mon_out_file_name,elpim_resp_pass,elpim_resp_pass_names)

elpim_restheta_pass = [hhs['el_res_theta'],hhs['pim_res_theta'],hhs['el_res_theta_fail_at_least1'],hhs['pim_res_theta_fail_at_least1']]
elpim_restheta_pass_names = ['El Pass EB \Delta \theta vs GEN \theta;  \theta GEN (deg) ; \Delta \theta','PIM Pass EB \Delta \theta vs GEN \theta;  \theta GEN (deg) ; \Delta \theta', 'El Fail EB \Delta \theta vs GEN \theta;  \theta GEN (deg) ; \Delta \theta', ' PIM Fail EB \Delta \theta vs GEN \theta;  \theta GEN (deg) ; \Delta \theta']
make_plots(can,2,2,mon_out_file_name,elpim_restheta_pass,elpim_restheta_pass_names)

elpim_resphi_pass = [hhs['el_res_phi'],hhs['pim_res_phi'],hhs['el_res_phi_fail_at_least1'],hhs['pim_res_phi_fail_at_least1']]
elpim_resphi_pass_names = ['El Pass EB \Delta  \phi vs GEN \phi;  \phi GEN (deg) ; \Delta \phi','PIM Pass EB \Delta \phi vs GEN \theta; \phi GEN (deg) ; \Delta \phi', 'El Fail EB \Delta \phi vs GEN \phi;  \phi GEN (deg) ; \Delta \phi', ' PIM Fail EB \Delta \phi vs GEN \phi;  \phi GEN (deg) ; \Delta \phi']
make_plots(can,2,2,mon_out_file_name,elpim_resphi_pass,elpim_resphi_pass_names)





elpim_nphe_pass=[hhs['el_nphe_passall'],hhs['pim_nphe_passall'],hhs['el_nphe_fail_at_least1'],hhs['pim_nphe_fail_at_least1']]
elpim_nphe_pass_names=['El Pass EB, Charge, Vz; nphe; counts', 'PIM Pass EB, Charge, Vz; nphe; counts', 'El Fail EB; nphe; counts', 'PIM Fail EB; nphe; counts']

make_plots(can,2,2,mon_out_file_name,elpim_nphe_pass,elpim_nphe_pass_names)

elpim_calsf_pass = [hhs['el_calsf_passall'],hhs['pim_calsf_passall'], hhs['el_calsf_fail_at_least1'],hhs['pim_calsf_fail_at_least1']]
elpim_calsf_pass_names = [' El pass EB; p; SF',' PIM pass EB; p; SF',' El fail EB; p; SF',' PIM fail EB; p; SF']

make_plots(can,2,2,mon_out_file_name,elpim_calsf_pass,elpim_calsf_pass_names)

elpim_pcalvsei_pass = [hhs['el_calsf_pcalvsei_passall'],hhs['pim_calsf_pcalvsei_passall'], hhs['el_calsf_pcalvsei_fail_at_least1'],hhs['pim_calsf_pcalvsei_fail_at_least1']]
elpim_pcalvsei_pass_names = [' El pass EB;PCAL/p; EI/p',' PIM pass EB;PCAL/p; EI/p',' El fail EB;PCAL/p; EI/p',' PIM fail EB; PCAL/p; EI/p']

line_eieo = TLine(0.0, 0.2, 0.2, 0.0)
line_eieo2 = TLine(0.0, 0.28, 0.28, 0.0)
make_plots(can,2,2,mon_out_file_name,elpim_pcalvsei_pass,elpim_pcalvsei_pass_names,x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[], lines_to_draw=[line_eieo])

elpim_eivseo_pass = [hhs['el_calsf_eivseo_passall'],hhs['pim_calsf_eivseo_passall'], hhs['el_calsf_eivseo_fail_at_least1'],hhs['pim_calsf_eivseo_fail_at_least1']]
elpim_eivseo_pass_names = [' El pass EB; EI/p; EO/p',' PIM pass EB; EI/p; EO/p',' El fail EB; EI/p; EO/p',' PIM fail EB; EI/p; EO/p']

make_plots(can,2,2,mon_out_file_name,elpim_eivseo_pass,elpim_eivseo_pass_names)

elpim_pcalvsecal_pass = [ hhs['el_calsf_pcalvsecal_passall'], hhs['pim_calsf_pcalvsecal_passall'], hhs['el_calsf_pcalvsecal_fail_at_least1'], hhs['pim_calsf_pcalvsecal_fail_at_least1'] ]
elpim_pcalvsecal_pass_names  = [ ' El Pass EB; PCAL/P; (EI+EO)/P','PIM Pass EB; PCAL/P; (EI+EO)/P', 'El Fail EB; PCAL/P; (EI+EO)/P', 'PIM Fail EB; PCAL/P; (EI+EO)/P' ]

make_plots(can,2,2,mon_out_file_name,elpim_pcalvsecal_pass,elpim_pcalvsecal_pass_names,x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[], lines_to_draw=[line_eieo,line_eieo2])

#not really chi2, it is a pull dist.

elpim_eivseo_passchi2 = []
elpim_eivseo_passchi2_names = []
for si in range(1,6):
    elpim_eivseo_passchi2.append(hhs['el_calsf_eivseo_passallpidchi2_'+str(si)+'sig'])
    elpim_eivseo_passchi2.append(hhs['pim_calsf_eivseo_passallpidchi2_'+str(si)+'sig'])
    #elpim_eivseo_passchi2.append(hhs['el_calsf_eivseo_passallpidchi2'+str(si)+'_fail_at_least1'])
    #elpim_eivseo_passchi2.append(hhs['pim_calsf_eivseo_passallpidchi2'+str(si)+'_fail_at_least1'])

    elpim_eivseo_passchi2_names.append('El Pass EB  - Chi Range '+str(si)+'; EI/p; EO/p')
    elpim_eivseo_passchi2_names.append('PIM Pass EB - Chi Range '+str(si)+'; EI/p; EO/p')
    #elpim_eivseo_passchi2_names.append('El fail EB - Chi Range '+str(si)+'; EI/p; EO/p')
    #elpim_eivseo_passchi2_names.append('PIM fail EB - Chi Range '+str(si)+'; EI/p; EO/p')

make_plots(can,2,5,mon_out_file_name,elpim_eivseo_passchi2,elpim_eivseo_passchi2_names)

elpim_cal_m2_pass = [hhs['el_pcal_m2u_passall'],hhs['el_eical_m2u_passall'],hhs['el_eocal_m2u_passall'], hhs['pim_pcal_m2u_passall'], hhs['pim_eical_m2u_passall'], hhs['pim_eocal_m2u_passall']]
elpim_cal_m2_pass_names = [' El Pass EB PCAL m2; M2PCAL; Counts', 'El Pass EB ECIN m2 ; M2ECIN; Counts', 'El Pass EB ECOUT m2; M2ECOUT; Counts', 'PIM Pass EB PCAL m2; M2PCAL; Counts', 'PIM Pass EB ECIN m2 ; M2ECIN; Counts', 'PIM Pass EB ECOUT m2; M2ECOUT; Counts']
make_plots(can,3,2,mon_out_file_name,elpim_cal_m2_pass,elpim_cal_m2_pass_names)



'''
elpim_eivseo_passp = []
elpim_eivseo_passp_names = []
for pp in range(10,100):
    if 'el_calsf_eivseo_prange'+str(pp)+'_passall' in hhs:
        elpim_eivseo_passp.append(hhs['el_calsf_eivseo_prange'+str(pp)+'_passall'])
        #elpim_eivseo_passp.append(hhs['pim_calsf_eivseo_prange'+str(pp)+'_passall'])
        #elpim_eivseo_passp.append(hhs['el_calsf_eivseo_prange'+str(pp)+'_fail_at_least1'])
        #elpim_eivseo_passp.append(hhs['pim_calsf_eivseo_prange'+str(pp)+'_fail_at_least1'])
        
        elpim_eivseo_passp_names.append('El Pass EB  - Mntm Range '+str(pp)+'; EI/p; EO/p')
        #elpim_eivseo_passp_names.append('PIM Pass EB - Mntm Range '+str(pp)+'; EI/p; EO/p')
        #elpim_eivseo_passp_names.append('El fail EB - Mntm Range '+str(pp)+'; EI/p; EO/p')
        #elpim_eivseo_passp_names.append('PIM fail EB - Mntm Range '+str(pp)+'; EI/p; EO/p')

make_plots(can,10,12,mon_out_file_name,elpim_eivseo_passp,elpim_eivseo_passp_names)
'''


elpim_pidchi2_pass = [hhs['el_calsfchi_passall'], hhs['pim_calsfchi_passall'], hhs['el_calsfchi_fail_at_least1'], hhs['pim_calsfchi_fail_at_least1']]
elpim_pidchi2_pass_names = [' El pass EB - Pull; pull; counts', 'PIM pass EB - Pull; pull ; counts ',' El fail EB - Pull; pull; counts', 'PIM fail EB - Pull; pull ; counts ' ] 

make_plots(can,2,2,mon_out_file_name,elpim_pidchi2_pass,elpim_pidchi2_pass_names)


total_number_of_electron_per_p = []# hhs['el_calsfchi_passall'].GetEntries()
total_number_of_pim_per_p =[] # hhs['pim_calsfchi_passall'].GetEntries()
total_number_of_elpim_per_p = [] #total_number_of_electron + total_number_of_pim

'''
n_el_per_p=0
n_pim_per_p=0
if "el_p_passall" in hhs:
    n_bins_el_p=hhs["el_p_passall"].GetXaxis().GetNbins()
    for bb in n_bins_el_p:
        n_el_bin_p = hhs["el_p_passall"].GetXaxis().GetBinContent()
        print(n_el_bin_p)

if "pim_p_passall" in hhs:
    n_bins_pim_p=hhs["pim_p_passall"].GetXaxis().GetNbins()
    for bb in n_bins_pim_p:
        n_pim_bin_p = hhs["pim_p_passall"].GetXaxis().GetBinContent()
        print(n_pim_bin_p)
'''     

hhs["pim_p_passall"].Add(hhs["el_p_passall"])
h_total_n_el_pim = hhs["pim_p_passall"]

el_temp_histo_p = []
can.Clear()
can.Divide(1,1)
can.cd(1)
for ss in reversed(range(1,6)):
    if 'el_p_passall_pidchi2'+str(ss) in hhs:
        hhs['el_p_passall_pidchi2'+str(ss)].Divide(h_total_n_el_pim)
        el_temp_histo_p.append(hhs['el_p_passall_pidchi2'+str(ss)])
        hhs['el_p_passall_pidchi2'+str(ss)].SetLineColor(kAzure-ss)
        hhs['el_p_passall_pidchi2'+str(ss)].Draw("same")
        
can.Print(mon_out_file_name)

pim_temp_histo_p = []
can.Clear()
can.Divide(1,1)
can.cd(1)

for ss in reversed(range(1,6)):
    if 'pim_p_passall_pidchi2'+str(ss) in hhs:
        hhs['pim_p_passall_pidchi2'+str(ss)].Divide(h_total_n_el_pim)
        el_temp_histo_p.append(hhs['pim_p_passall_pidchi2'+str(ss)])
        hhs['pim_p_passall_pidchi2'+str(ss)].SetLineColor(kAzure-ss)
        hhs['pim_p_passall_pidchi2'+str(ss)].Draw("same")
        
can.Print(mon_out_file_name)
 
# now determine the contamination of pions in eletron sample
temp_contamination = []
fin = open("eToPRatioPythia.txt","r")
fin_elpiratio=fin.readlines()
el_pion_rates=[]
for ll in fin_elpiratio:
    print ll.split()
    ll_split = ll.split()
    pions_to_el = 0
    if float(ll_split[1]) != 0.0 :
        pions_to_el = float(ll_split[2])/float(ll_split[1])
    print(pions_to_el)
    el_pion_rates.append(pions_to_el)
    
################################################################
## raw contamination without any additional pion rejection cuts
hh_el_cont_sect = []
hh_el_cont_sect_name = []
for ss in range(1,6):
    hh_temp = hhs['el_p_passall_pidchi2'+str(ss)]
    hh_cont = TH1F('Electron Contamination CHI' + str(ss), 'Electron Contamination ' + str(ss), 100, 0.0, 10.0);
    
    print('sector %d ' % ss )
    for bb in range(1, hh_temp.GetXaxis().GetNbins()):
        pion_as_el_ratio = hh_temp.GetBinContent(bb)
        print(' bin center %f ' % hh_temp.GetXaxis().GetBinCenter(bb))
        print(' pion as el ratio is %f ' % pion_as_el_ratio)
        print(' pion / el rates  %f ' % el_pion_rates[bb])        
        print(' bin center %f - pion to el ratio %f times pion/el %f ' % (hh_temp.GetXaxis().GetBinCenter(bb), pion_as_el_ratio, el_pion_rates[bb]) )
        hh_cont.SetBinContent(int(bb), pion_as_el_ratio*el_pion_rates[bb])
    hh_el_cont_sect.append(hh_cont)
    hh_el_cont_sect_name.append('Pion As El. Contamination PIDChi2 ' + str(ss) +'; p (GeV); Contamination %')

make_plots(can,2,3,mon_out_file_name, hh_el_cont_sect ,hh_el_cont_sect_name)


###################################################################
## contamination plots with pion rejection cut 


el_pion_rates_withcut=[hhs['el_calsf_pcalvsei_passall_rej_pions'], hhs['pim_calsf_pcalvsei_passall_rej_pions'],hhs['el_calsf_pcalvsei_passall_pass_pions'],hhs['pim_calsf_pcalvsei_passall_pass_pions']]
el_pion_rates_withcut_names=['El with New Cut; PCAL/P; EI/P','Pion with New Cut; PCAL/P; EI/P','El Fail New Cut;PCAL/P; EI/P','Pion Fail New Cut; PCAL/P; EI/P']
make_plots(can, 2, 2, mon_out_file_name, el_pion_rates_withcut, el_pion_rates_withcut_names)

el_pion_p_withcut=[hhs['el_p_passall_rej_pions'],hhs['pim_p_passall_rej_pions'],hhs['el_p_passall_pass_pions'],hhs['pim_p_passall_pass_pions']]
el_pion_p_withcut_names=['El Pass New Cut P; p (GeV); Counts','Pion Pass New Cut P; p (GeV); Counts','El Fail New Cut P; p (GeV); counts','Pion Fail New Cut P; p (GeV); counts']
make_plots(can,2,2,mon_out_file_name,el_pion_p_withcut,el_pion_p_withcut_names)


el_temp_histo_p_rej_pions = []
can.Clear()
can.Divide(1,1)
can.cd(1)
for ss in reversed(range(1,6)):
    if "el_p_passall_pass_pions_pidchi2_"+str(ss)+"sig" in hhs:
        hhs["el_p_passall_pass_pions_pidchi2_"+str(ss)+"sig"].Divide(h_total_n_el_pim)
        el_temp_histo_p.append(hhs["el_p_passall_pass_pions_pidchi2_"+str(ss)+"sig"])
        hhs["el_p_passall_pass_pions_pidchi2_"+str(ss)+"sig"].SetLineColor(kSpring-ss)
        hhs["el_p_passall_pass_pions_pidchi2_"+str(ss)+"sig"].Draw("same")
        
can.Print(mon_out_file_name)


el_p_withcut_pidchi2=[]
can.Clear()
can.Divide(1,1)
can.cd(1)
for ss in reversed(range(1,6)):
    if "el_p_passall_rej_pions_pidchi2_"+str(ss)+"sig" in hhs:
        hhs["el_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"].Divide(h_total_n_el_pim)
        el_p_withcut_pidchi2.append(hhs["el_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"])
        hhs["el_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"].SetLineColor(kSpring-ss)
        hhs["el_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"].Draw("same")
        
        
can.Print(mon_out_file_name)

pion_p_withcut_pidchi2=[]
can.Clear()
can.Divide(1,1)
can.cd(1)
for ss in reversed(range(1,6)):
    if "pim_p_passall_rej_pions_pidchi2_"+str(ss)+"sig" in hhs:
        hhs["pim_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"].Divide(h_total_n_el_pim)
        pion_p_withcut_pidchi2.append(hhs["pim_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"])
        hhs["pim_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"].SetLineColor(kSpring-ss)
        hhs["pim_p_passall_rej_pions_pidchi2_"+str(ss)+"sig"].Draw("same")
                
can.Print(mon_out_file_name)

################
## for electrons we want to check how many are rejected out by this cut to the total that would have been kept
can.Clear()
can.Divide(1,1)
can.cd(1)

h_temp_ratio_of_electrons_kept = TH1F('Electron Rates Passed New Cut ' , 'Electrons Rates Passed New Cut; p (GeV); N electrons good/ N pass EB ', 100, 0.0, 10.0);
h_temp_ratio_of_electrons_rej = TH1F('Electron Rates Failed New Cut ' , 'Electrons Rates Failed New Cut; p (GeV); N electrons good but fail/ N pass EB ', 100, 0.0, 10.0);
for bb in range(1, hhs["el_p_passall"].GetXaxis().GetNbins() ):
    if hhs["el_p_passall"].GetBinContent(bb) > 0 :
        h_temp_ratio_of_electrons_kept.SetBinContent(bb,  hhs["el_p_passall_pass_pions"].GetBinContent(bb)/hhs["el_p_passall"].GetBinContent(bb) )
        h_temp_ratio_of_electrons_rej.SetBinContent(bb,  hhs["el_p_passall_rej_pions"].GetBinContent(bb)/hhs["el_p_passall"].GetBinContent(bb) )

h_temp_ratio_of_electrons_kept.SetLineColor(kRed)
h_temp_ratio_of_electrons_rej.SetLineColor(kBlue)
h_temp_ratio_of_electrons_kept.Draw()
h_temp_ratio_of_electrons_rej.Draw("same")

can.Print(mon_out_file_name)

hh_el_cont_withcut = []
hh_el_cont_withcut_name = []
for ss in range(1,6):
    hh_temp = hhs["el_p_passall_pass_pions_pidchi2_"+str(ss)+"sig"]
    hh_cont = TH1F('Electron Contamination Rej. Pions CHI' + str(ss), 'Electron Rej. Pions Contamination ' + str(ss), 100, 0.0, 10.0);
    
    print('pidchi2 cut lvl %d ' % ss )
    for bb in range(1, hh_temp.GetXaxis().GetNbins()):
        pion_as_el_ratio = hh_temp.GetBinContent(bb)
        print(' bin center %f ' % hh_temp.GetXaxis().GetBinCenter(bb))
        print(' pion as el ratio is %f ' % pion_as_el_ratio)
        print(' pion / el rates  %f ' % el_pion_rates[bb])        
        print(' bin center %f - pion to el ratio % f times pion/el %f ' % (hh_temp.GetXaxis().GetBinCenter(bb), pion_as_el_ratio, el_pion_rates[bb]) )
        hh_cont.SetBinContent(int(bb), pion_as_el_ratio*el_pion_rates[bb])
    hh_el_cont_withcut.append(hh_cont)
    hh_el_cont_withcut_name.append('Pion As El. Contamination Rej Pions PIDChi2 ' + str(ss) +'; p (GeV); Contamination %')

make_plots(can,2,3,mon_out_file_name, hh_el_cont_withcut ,hh_el_cont_withcut_name)


'''
mg_counts_per_chi2 = TMultiGraph()


ii=0
list_save=[]
for cc in counts_per_sigma_over_mntm:
    g_temp = TGraph(len(cc),mntm_range_per_sigma[ii],cc)
    #print(mntm_range_per_sigma[ii])
    g_temp.SetLineColor(kRed+ii)
    g_temp.SetTitle("sigma range " + str(cc))
    mg_counts_per_chi2.Add(g_temp,'lp')
    ii+=1
    list_save.append(g_temp)

can.Clear()
can.Divide(1,1)
can.cd(1)
mg_counts_per_chi2.SetTitle("EL/(El+PIM) Per PIDChi2 vs Momentum; p (GeV); el/(el+pim)")
mg_counts_per_chi2.Draw("AP")
can.Print(mon_out_file_name)

pim_mg_counts_per_chi2 = TMultiGraph()

pim_ii=0
pim_list_save=[]
for pim_cc in pim_counts_per_sigma_over_mntm:
    
    pim_g_temp = TGraph(len(pim_cc),pim_mntm_range_per_sigma[pim_ii],pim_cc)
    
    pim_g_temp.SetLineColor(kRed+pim_ii)
    pim_g_temp.SetTitle("pim sigma range " + str(cc))
    pim_mg_counts_per_chi2.Add(pim_g_temp,'lp')
    pim_ii+=1
    pim_list_save.append(pim_g_temp)

can.Clear()
can.Divide(1,1)
can.cd(1)
pim_mg_counts_per_chi2.SetTitle("PIM/(El+PIM) Per PIDChi2 vs Momentum; p (GeV); pim/(el+pim)")
pim_mg_counts_per_chi2.Draw("AP")
can.Print(mon_out_file_name)
'''



can.Print('{}]'.format(mon_out_file_name))

