from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed
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

        if draw_cut:
            can.Update()
            pad=can.GetPad(cc)
            ymax=gPad.GetUymax()
            cutline = TLine(cut_pos[0], 0.0, cut_pos[1],ymax)
            cutline.SetLineColor(kRed)
            cutline.Draw('same')
            root_is_dumb.append(cutline)

        if len(lines_to_draw) > 0 :
            for ll in lines_to_draw:
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


def make_plots_per_sector_with_cut_lines(can, fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False,cut_lines=[]):
    can.Clear()
    can.Divide(3,6)
    cc=1
    root_is_dumb = []
    ss=0
    for hh in hists:
        print('drawing to canvas {}'.format(cc))
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.Draw('colz')
        gPad.SetLogz()
        print('drawing sector {}'.format(ss))
        print(hh.GetTitle())
        #print('draw cut %d' %(draw_cut ))
        #print('cc mod 3 %d' %( cc%2 ))
        if draw_cut :
            if len(cut_lines[ss]) >= 2:
                cp=cut_lines[ss]
                   
                cut_min=cp[0]
                cut_max=cp[1]
                
                cut_min.SetLineColor(kRed)
                cut_max.SetLineColor(kRed)
                cut_min.Draw('same')
                cut_max.Draw('same')
                #root_is_dumb.append(cut_min)
                #root_is_dumb.append(cut_max)
            #else:
            #    cutline = TLine(cut_pos[ss][0], 0.0, cut_pos[ss][0],ymax)
            #    cutline.SetLineColor(kRed)
            #    cutline.Draw('same')
            #    root_is_dumb.append(cutline)

        root_is_dumb.append(hh)
        
        if cc%3 == 0:             
            ss+=1  
            if ss == 6:
                ss=0
        cc+=1

    can.Print(fname)

###################################
## beta vs p for hadrons
###################################
pr_betap = TF1("pr_betap","x/TMath::Sqrt(x*x + 0.938*0.938)",0.0, 6.0)
kp_betap = TF1("kp_betap","x/TMath::Sqrt(x*x + 0.493*0.493)",0.0, 6.0)
pip_betap = TF1("pip_betap","x/TMath::Sqrt(x*x + 0.139*0.139)",0.0, 6.0)
du_betap = TF1("du_betap","x/TMath::Sqrt(x*x + 1.8756*1.8756)",0.0, 6.0)
pr_betap.SetLineColor(kRed)
kp_betap.SetLineColor(kRed)
pip_betap.SetLineColor(kRed)
du_betap.SetLineColor(kRed)
pr_betap.SetLineWidth(1)
kp_betap.SetLineWidth(1)
pip_betap.SetLineWidth(1)
du_betap.SetLineWidth(1)
hadron_betaps = [pr_betap, kp_betap, pip_betap, du_betap]


#gStyle.SetOptStat(00000)
gStyle.SetPalette(55)

beamE = 10.604
hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj

mon_out_file_name='monitor_timing.pdf'
can = TCanvas('can','can',4000,4000)
can.Print('{}['.format(mon_out_file_name))

#plot Iron Cross Time
hcross1=[]
hcross1_titles=[]



for ss in range(1,7):
    hcross1_titles.append('Timing Cross for FTOF and PCAL'+ str(ss) +'; #Delta t FTOF; #Delta t PCAL')
    hcross1.append(hhs['iron_cross_xscint_ypcal_s'+str(ss)+'passall'])

can.SetCanvasSize(4000,4000)
make_plots(can,3,2,mon_out_file_name,hcross1,hcross1_titles,draw_cut=False, cut_pos=[])

                        
hcross2=[]
hcross2_titles=[]

for ss in range(1,7):
    hcross2_titles.append('Timing Cross for HTCC and PCAL'+ str(ss) +'; #Delta t HTCC; #Delta t PCAL')
    hcross2.append(hhs['iron_cross_xhtcc_ypcal_s'+str(ss)+'passall'])

can.SetCanvasSize(4000,4000)
make_plots(can,1,1,mon_out_file_name,[hhs['iron_cross_xhtcc_ypcal_passall']],['Iron Cross HTCC vs PCAL; #Delta t FTOF; #Delta t PCAL'],draw_cut=False, cut_pos=[])

can.SetCanvasSize(4000,4000)
make_plots(can,3,2,mon_out_file_name,hcross2,hcross2_titles,draw_cut=False, cut_pos=[])
                                                                             

hcross3=[]
hcross3_titles=[]

for ss in range(1,7):
    hcross3_titles.append('Timing Cross for HTCC and FTOF'+ str(ss) +'; #Delta t HTCC; #Delta t FTOF')
    hcross3.append(hhs['iron_cross_xhtcc_yftof_s'+str(ss)+'passall'])

can.SetCanvasSize(4000,4000)
make_plots(can,1,1,mon_out_file_name,[hhs['iron_cross_xhtcc_yftof_passall']],['Cross HTCC vs FTOF; #Delta t HTCC; #Delta t FTOF'],draw_cut=False, cut_pos=[])

can.SetCanvasSize(4000,4000)
make_plots(can,3,2,mon_out_file_name,hcross3,hcross3_titles,draw_cut=False, cut_pos=[])


#######
## make plots for proton multiplicities with and without timing cuts
#######

ftof_pcal_multiplicities = [hhs["iron_cross_xftof_ypcal_passall_proton_matched"],
                            hhs["hadron_multiplicities_proton_ftof_pcal_timing_matched"],
                            hhs["iron_cross_xftof_ypcal_passall_proton_unmatched"],
                            hhs["hadron_multiplicities_proton_ftof_pcal_timing_unmatched"]
                            

                        ]
ftof_pcal_multiplicities_titles = [' Cross FTOF vs PCAL Selected Electrons; #Delta t FTOF (ns); #Delta t PCAL (ns)',
                                   ' Proton Multiplicities for FTOF vs PCAL Selected Electrons; #Number of protons; counts',
                                   ' Cross FTOF vs PCAL UN-selected Electrons; #Delta t FTOF (ns); #Delta t PCAL (ns)',
                                   ' Proton Multiplicities for FTOF vs PCAL UN-selected Electrons; #Number of protons; counts'
                                   ]

can.SetCanvasSize(4000,4000)
make_plots(can,2,2,mon_out_file_name,ftof_pcal_multiplicities,ftof_pcal_multiplicities_titles,draw_cut=False, cut_pos=[])


htcc_pcal_multiplicities = [hhs["iron_cross_xhtcc_ypcal_passall_proton_matched"],
                            hhs["hadron_multiplicities_proton_htcc_pcal_timing_matched"],
                            hhs["iron_cross_xhtcc_ypcal_passall_proton_unmatched"],
                            hhs["hadron_multiplicities_proton_htcc_pcal_timing_unmatched"]                            
                        ]
htcc_pcal_multiplicities_titles = [' Cross HTCC vs PCAL Selected Electrons; #Delta t HTCC (ns); #Delta t PCAL (ns)',
                                   ' Proton Multiplicities for HTCC vs PCAL Selected Electrons; #Number of protons; counts',
                                   ' Cross HTCC vs PCAL UN-selected Electrons; #Delta t HTCC (ns); #Delta t PCAL (ns)',
                                   ' Proton Multiplicities for HTCC vs PCAL UN-selected Electrons; #Number of protons; counts'
                                   ]

can.SetCanvasSize(4000,4000)
make_plots(can,2,2,mon_out_file_name,htcc_pcal_multiplicities,htcc_pcal_multiplicities_titles,draw_cut=False, cut_pos=[])


htcc_ftof_multiplicities = [hhs["iron_cross_xhtcc_yftof_passall_proton_matched"],
                            hhs["hadron_multiplicities_proton_htcc_ftof_timing_matched"],
                            hhs["iron_cross_xhtcc_yftof_passall_proton_unmatched"],
                            hhs["hadron_multiplicities_proton_htcc_ftof_timing_unmatched"]                            
                        ]
htcc_ftof_multiplicities_titles = [' Cross HTCC vs FTOF Selected Electrons; #Delta t HTCC (ns); #Delta t FTOF (ns)',
                                   ' Proton Multiplicities for FTOF vs HTCC Selected Electrons; #Number of protons/event; counts',
                                   ' Cross HTCC vs FTOF UN-selected Electrons; #Delta t HTCC (ns); #Delta t PCAL (ns)',
                                   ' Proton Multiplicities for FTOF vs HTCC UN-selected Electrons; #Number of protons/event; counts'
                                   ]

can.SetCanvasSize(4000,4000)
make_plots(can,2,2,mon_out_file_name,htcc_ftof_multiplicities,htcc_ftof_multiplicities_titles,draw_cut=False, cut_pos=[])


####################
## show rec chi2pid 

ftof_pcal_chi2pid_tmatched=[hhs['recpart_chi2pid_proton_pid_ftof_pcal_goodelectron2D'],
                            hhs['recpart_chi2pid_proton_pid_ftof_pcal_goodelectron_tmatched'],
                            hhs['recpart_chi2pid_proton_pid_ftof_pcal_goodelectron_tunmatched']
                            ]
ftof_pcal_chi2pid_tmatched_titles = ['Proton Chi2 - Electron FTOF PCAL #Delta t All; p (GeV); #chi 2',
                                     'Proton Chi2 - FTOF PCAL #Delta t Selected Electron ; #chi 2',
                                     'Proton Chi2 - FTOF PCAL #Delta t UN-selected Electron; #chi 2'
                                     ]

can.SetCanvasSize(4000,1300)
make_plots(can,3,1,mon_out_file_name,ftof_pcal_chi2pid_tmatched,ftof_pcal_chi2pid_tmatched_titles,draw_cut=False, cut_pos=[])

htcc_pcal_chi2pid_tmatched=[hhs['recpart_chi2pid_proton_pid_htcc_pcal_goodelectron2D'],
                            hhs['recpart_chi2pid_proton_pid_htcc_pcal_goodelectron_tmatched'],
                            hhs['recpart_chi2pid_proton_pid_htcc_pcal_goodelectron_tunmatched']
                            ]
htcc_pcal_chi2pid_tmatched_titles = ['Proton Chi2 - Electron HTCC PCAL #Delta t All; p (GeV); #chi 2',
                                     'Proton Chi2 - HTCC PCAL #Delta t Selected Electron ; #chi 2',
                                     'Proton Chi2 - HTCC PCAL #Delta t UN-selected Electron; #chi 2'
                                     ]

can.SetCanvasSize(4000,1300)
make_plots(can,3,1,mon_out_file_name,htcc_pcal_chi2pid_tmatched,htcc_pcal_chi2pid_tmatched_titles,draw_cut=False, cut_pos=[])

ftof_htcc_chi2pid_tmatched=[hhs['recpart_chi2pid_proton_pid_htcc_ftof_goodelectron2D'],
                            hhs['recpart_chi2pid_proton_pid_htcc_ftof_goodelectron_tmatched'],
                            hhs['recpart_chi2pid_proton_pid_htcc_ftof_goodelectron_tunmatched']
                            ]
ftof_htcc_chi2pid_tmatched_titles = ['Proton Chi2 - Electron FTOF HTCC #Delta t All; p (GeV); #chi 2',
                                     'Proton Chi2 - FTOF HTCC #Delta t Selected Electron ; #chi 2',
                                     'Proton Chi2 - FTOF HTCC #Delta t UN-selected Electron; #chi 2'
                                     ]

can.SetCanvasSize(4000,1300)
make_plots(can,3,1,mon_out_file_name,ftof_htcc_chi2pid_tmatched,ftof_htcc_chi2pid_tmatched_titles,draw_cut=False, cut_pos=[])

# matching htcc pcal ftof
pr_matching_all_det_boxcut = [hhs['hadron_multiplicities_proton_htcc_ftof_timing_box_cut_matched'],
                              hhs['iron_cross_xhtcc_yftof_passall_proton_box_cut_matched'],
                              hhs['iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_box_cut_matched'],
                              hhs['hadron_multiplicities_proton_htcc_ftof_timing_box_cut_unmatched'],
                              hhs['iron_cross_xhtcc_yftof_passall_proton_box_cut_unmatched'],
                              hhs['iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_box_cut_unmatched'],
                              ]

pr_matching_all_det_boxcut_titles = ['Proton Multiplicities Pass BoxCut; protons/ev; ',
                                     'Cross #Delta HTCC vs #Delta FTOF Pass BoxCut ; HTCC ; FTOF',
                                     'Cross #Delta FTOF - PCAL v.s. #Delta HTCC - FTOF Pass BoxCut; #Delta FTOF - PCAL; #Delta HTCC - FTOF',
                                     'Proton Multiplicities Fail BoxCut; protons/ev; ',
                                     'Cross #Delta HTCC vs #Delta FTOF Fail BoxCut ; HTCC ; FTOF',
                                     'Cross #Delta FTOF - PCAL v.s. #Delta HTCC - FTOF Fail BoxCut; #Delta FTOF - PCAL; #Delta HTCC - FTOF']

can.SetCanvasSize(4000,1300)
make_plots(can,3,2,mon_out_file_name,pr_matching_all_det_boxcut, pr_matching_all_det_boxcut_titles,draw_cut=False, cut_pos=[])


#############
## corner cut

iron_cross_corner_cut_time_diff = [hhs['iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall'],
                                   hhs['iron_cross_ftof_minus_pcal_vs_ftof_minus_htcc_passall_corner_cut_matched'],
                                   hhs['iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall_corner_cut_unmatched']
                                   ]

iron_cross_corner_cut_time_diff_titles=['Cross #Delta FTOF - PCAL v.s #Delta FTOF - HTCC All; #Delta FTOF - PCAL ; #Delta FTOF - HTCC ',
                                        'Cross #Delta FTOF - PCAL v.s #Delta FTOF - HTCC Pass CornerCut; #Delta FTOF - PCAL ; #Delta FTOF - HTCC ',
                                        'Cross #Delta FTOF - PCAL v.s #Delta FTOF - HTCC Fail CornerCut; #Delta FTOF - PCAL ; #Delta FTOF - HTCC '
                                        ]

can.SetCanvasSize(4000,1350)
make_plots(can,3,1,mon_out_file_name, iron_cross_corner_cut_time_diff, iron_cross_corner_cut_time_diff_titles,draw_cut=False, cut_pos=[])


# matching htcc pcal ftof
pr_matching_all_det_cornercut = [hhs['hadron_multiplicities_proton_htcc_ftof_timing_corner_cut_matched'],
                              hhs['iron_cross_xhtcc_yftof_passall_proton_corner_cut_matched'],
                              hhs['iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_matched'],
                              hhs['iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_matched'],
                              hhs['iron_cross_ftof_minus_pcal_vs_ftof_minus_htcc_passall_corner_cut_matched'],
                              hhs['hadron_multiplicities_proton_htcc_ftof_timing_corner_cut_unmatched'],                  
                              hhs['iron_cross_xhtcc_yftof_passall_proton_corner_cut_unmatched'],
                              hhs['iron_cross_ftof_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_unmatched'],
                              hhs['iron_cross_htcc_minus_pcal_vs_htcc_minus_ftof_passall_corner_cut_unmatched'],
                              hhs['iron_cross_ftof_minus_pcal_vs_ftof_minus_pcal_passall_corner_cut_unmatched']
                             ]

pr_matching_all_det_cornercut_titles = ['Proton Multiplicities Pass CornerCut; protons/ev; ',
                                     'Cross #Delta HTCC vs #Delta FTOF Pass CornerCut ; HTCC ; FTOF',
                                     'Cross #Delta FTOF - PCAL v.s. #Delta HTCC - FTOF Pass CornerCut; #Delta FTOF - PCAL; #Delta HTCC - FTOF',
                                     'Cross #Delta HTCC - PCAL v.s. #Delta HTCC - FTOF Pass CornerCut; #Delta HTCC - PCAL; #Delta HTCC - FTOF',
                                     'Cross #Delta FTOF - PCAL v.s. #Delta FTOF - HTCC Pass CornerCut; #Delta FTOF - PCAL; #Delta FTOF - HTCC',
                                     'Proton Multiplicities Fail CornerCut; protons/ev; ',
                                     'Cross #Delta HTCC vs #Delta FTOF Fail CornerCut ; HTCC ; FTOF',
                                     'Cross #Delta FTOF - PCAL v.s. #Delta HTCC - FTOF Fail CornerCut; #Delta FTOF - PCAL; #Delta HTCC - FTOF',
                                     'Cross #Delta HTCC - PCAL v.s. #Delta HTCC - FTOF Fail CornerCut; #Delta HTCC - PCAL; #Delta HTCC - FTOF',
                                     'Cross #Delta FTOF - PCAL v.s. #Delta FTOF - HTCC Fail CornerCut; #Delta FTOF - PCAL; #Delta FTOF - HTCC'                                        
                                    ]



can.SetCanvasSize(4000,1600)
make_plots(can,5,2,mon_out_file_name,pr_matching_all_det_cornercut, pr_matching_all_det_cornercut_titles,draw_cut=False, cut_pos=[])
              

          
# refining beta plots
pos_beta_redefined = [hhs['pos_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched'],
                      hhs['pos_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                      hhs['pos_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                      hhs['pos_beta_p_redfined_ftof_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                      hhs['pos_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                      hhs['pos_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched']
                      ]

pos_beta_redefined_names = ['Positives EB Beta vs P; p ; beta',
                            'Positives PCAL VzCorrStartTime - Beta vs p; p ; beta',
                            'Positives HTCC VzCorrStartTime - Beta vs p; p; beta',
                            'Positives FTOFL2 VzCorrStartTime - Beta vs p; p; beta',
                            'Positives EB StartTime - Beta vs p; p; beta',
                            'Positives EB Vt - Beta vs p; p; beta']
                            
can.SetCanvasSize(1000,1000)
make_plots(can,2,3,mon_out_file_name,pos_beta_redefined,pos_beta_redefined_names, draw_cut=False, cut_pos=[], lines_to_draw=hadron_betaps)


#unmatched - failed box cut 
pr_beta_redefined_set1 = [ hhs['proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched'],
                          hhs['proton_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                          hhs['proton_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                          hhs['proton_beta_p_redfined_ftof_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                          hhs['proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched'],
                          hhs['proton_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched']
                      ]
                          

pr_beta_redefined_names_set1 = ['EB Beta vs P Failed Time BoxCut; p ; beta',
                            'PCAL VzCorrStartTime - Beta vs p Failed Time BoxCut; p ; beta',
                            'HTCC VzCorrStartTime - Beta vs p Failed Time BoxCut; p; beta',
                            'FTOFL2 VzCorrStartTime - Beta vs p Failed Time BoxCut; p; beta',
                            'EB StartTime - Beta vs p Failed Time BoxCut; p; beta',
                            'EB Vt - Beta vs p Failed Time BoxCut; p; beta']

can.SetCanvasSize(1000,1000)
make_plots(can,2,3,mon_out_file_name,pr_beta_redefined_set1,pr_beta_redefined_names_set1, draw_cut=False, cut_pos=[], lines_to_draw=hadron_betaps)

#matched - passed box cut 
pr_beta_redefined_set2 = [ hhs['proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tmatched'],
                          hhs['proton_beta_p_redefined_pcal_st_htcc_ftof_pcal_goodelectron_tmatched'],
                          hhs['proton_beta_p_redefined_htcc_st_htcc_ftof_pcal_goodelectron_tmatched'],
                          hhs['proton_beta_p_redefined_ftof_st_htcc_ftof_pcal_goodelectron_tmatched'],
                          hhs['proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tmatched'],
                          hhs['proton_beta_p_redefined_vt_st_htcc_ftof_pcal_goodelectron_tmatched']
                      ]
                          

pr_beta_redefined_names_set2 = ['EB Beta vs P Passed Time BoxCut; p ; beta',
                            'PCAL VzCorrStartTime - Beta vs p Passed Time BoxCut; p ; beta',
                            'HTCC VzCorrStartTime - Beta vs p Passed Time BoxCut; p; beta',
                            'FTOFL2 VzCorrStartTime - Beta vs p Passed Time BoxCut; p; beta',
                            'EB StartTime - Beta vs p Passed Time BoxCut; p; beta',
                            'EB Vt - Beta vs p Passed Time BoxCut; p; beta']

can.SetCanvasSize(1000,1000)
make_plots(can,2,3,mon_out_file_name,pr_beta_redefined_set2,pr_beta_redefined_names_set2, draw_cut=False, cut_pos=[], lines_to_draw=hadron_betaps)
                            

#events passing the corner cut
pr_beta_redefined_corner_set1 = [hhs['proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tmatched_corner'],
                                 hhs['proton_beta_p_redefined_pcal_st_htcc_ftof_pcal_goodelectron_tmatched_corner'],
                                 hhs['proton_beta_p_redefined_htcc_st_htcc_ftof_pcal_goodelectron_tmatched_corner'],
                                 hhs['proton_beta_p_redefined_ftof_st_htcc_ftof_pcal_goodelectron_tmatched_corner'],
                                 hhs['proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tmatched_corner'],
                                 hhs['proton_beta_p_redefined_vt_st_htcc_ftof_pcal_goodelectron_tmatched_corner']
                                 ]

pr_beta_redefined_corner_names_set1 = ['EB Beta vs P Passed Time CornerCut; p ; beta',
                            'PCAL VzCorrStartTime - Beta vs p Passed Time CornerCut; p ; beta',
                            'HTCC VzCorrStartTime - Beta vs p Passed Time CornerCut; p; beta',
                            'FTOFL2 VzCorrStartTime - Beta vs p Passed Time CornerCut; p; beta',
                            'EB StartTime - Beta vs p Passed Time CornerCut; p; beta',
                            'EB Vt - Beta vs p Passed Time CornerCut; p; beta']
                                 
can.SetCanvasSize(1000,1000)
make_plots(can,2,3,mon_out_file_name,pr_beta_redefined_corner_set1,pr_beta_redefined_corner_names_set1, draw_cut=False, cut_pos=[], lines_to_draw=hadron_betaps)

#events not in the corner (failed corner cut)
pr_beta_redefined_corner_set2 = [hhs['proton_beta_p_eb_htcc_ftof_pcal_goodelectron_tunmatched_corner'],
                                 hhs['proton_beta_p_redfined_pcal_st_htcc_ftof_pcal_goodelectron_tunmatched_corner'],
                                 hhs['proton_beta_p_redfined_htcc_st_htcc_ftof_pcal_goodelectron_tunmatched_corner'],
                                 hhs['proton_beta_p_redfined_ftof_st_htcc_ftof_pcal_goodelectron_tunmatched_corner'],
                                 hhs['proton_beta_p_redefined_event_st_htcc_ftof_pcal_goodelectron_tunmatched_corner'],
                                 hhs['proton_beta_p_redfined_vt_st_htcc_ftof_pcal_goodelectron_tunmatched_corner']
                                 ]

pr_beta_redefined_corner_names_set2 = ['EB Beta vs P Failed Time CornerCut; p ; beta',
                            'PCAL VzCorrStartTime - Beta vs p Failed Time CornerCut; p ; beta',
                            'HTCC VzCorrStartTime - Beta vs p Failed Time CornerCut; p; beta',
                            'FTOFL2 VzCorrStartTime - Beta vs p Failed Time CornerCut; p; beta',
                            'EB StartTime - Beta vs p Failed Time CornerCut; p; beta',
                            'EB Vt - Beta vs p Failed Time CornerCut; p; beta']
                                 
can.SetCanvasSize(1000,1000)
make_plots(can,2,3,mon_out_file_name,pr_beta_redefined_corner_set2,pr_beta_redefined_corner_names_set2, draw_cut=False, cut_pos=[], lines_to_draw=hadron_betaps)




can.Print('{}]'.format(mon_out_file_name))
