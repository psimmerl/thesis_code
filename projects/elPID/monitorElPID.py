from ROOT import TFile, TH1F, TH2F, TF1, TLine
from ROOT import TCanvas, gStyle, gPad
from ROOT import kRed, kGreen, kBlue, kPink, kBlack, kViolet
import sys,os


def make_plots(can, canx, cany, fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[], out_file_name=''):
    can.Clear()
    can.Divide(canx,cany)
    cc=1
    root_is_dumb = []
    for hh in hists:

        gStyle.SetOptStat(0)#"e")
        gStyle.SetStatY(0.85);      
        gStyle.SetStatX(0.8);  
        gStyle.SetStatW(0.1);                
        gStyle.SetTextFont(102)
        
        gPad.SetLogz()
        gPad.SetLeftMargin(0.15)
        gPad.SetBottomMargin(0.2)
        #gPad.SetRightMargin(0.15)
        
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.Draw('colz')
        hh.SetLineColor(kBlack)
        hh.GetXaxis().SetTitleSize(0.06);
        hh.GetYaxis().SetTitleSize(0.06);
        hh.GetXaxis().SetLabelOffset(0.001);
        hh.GetXaxis().SetLabelSize(0.05)
        hh.GetYaxis().SetLabelSize(0.05)
        
        can.Update()

        if draw_cut:
            can.Update()
            pad=can.GetPad(cc)
            ymax=gPad.GetUymax()
            cutline = TLine(cut_pos[0], 0.0, cut_pos[1],ymax)
            cutline.SetLineColor(kRed)
            cutline.Draw('same')
            root_is_dumb.append(cutline)

        root_is_dumb.append(hh)
        cc+=1
    can.SaveAs(out_file_name)
    can.Print(fname)

def make_plots_overlap1D(can, canx, cany, fname, hists, hh_titles,draw_cut=False, cut_pos=[],  out_file_name=''):
    can.Clear()
    can.Divide(1,1)
    cc=1
    root_is_dumb = []
    hist_colors = [kRed, kGreen, kBlue, kPink, kBlack, kViolet]
    can.cd(1)
    for hh in hists:
        
        hh.SetTitle(hh_titles[cc-1])
        hh.SetLineColor(hist_colors[cc-1])
        hh.Draw('same')    
        hh.GetXaxis().SetTitleSize(0.05);
        hh.GetYaxis().SetTitleSize(0.05);
        hh.GetXaxis().SetLabelOffset(0.001);
        hh.GetXaxis().SetLabelSize(0.05)
        hh.GetYaxis().SetLabelSize(0.05)

        gStyle.SetOptStat(0)#"e")
        gStyle.SetStatY(0.85);      
        gStyle.SetStatX(0.8);  
        gStyle.SetStatW(0.1);                

        
        gPad.SetLogz()
        gPad.SetLeftMargin(0.15)
        gPad.SetBottomMargin(0.2)
        #gPad.SetRightMargin(0.15)
        
        can.Update()
        pad=can.GetPad(1)
        ymax=gPad.GetUymax()
        
        hh.GetYaxis().SetRangeUser(0., ymax+100);


        root_is_dumb.append(hh)
        cc+=1



    if draw_cut:
        can.Update()
        pad=can.GetPad(1)
        ymax=gPad.GetUymax()
        for cp in cut_pos[0]: # if all same per sector
            cutline = TLine(cp, 0.0, cp,ymax)
            cutline.SetLineColor(kRed)
            cutline.Draw('same')
            root_is_dumb.append(cutline)
    
    can.SaveAs(out_file_name)        
    can.Print(fname)

def make_plots_overlap(can, canx, cany, fname, hists, hh_titles, cut_before, set_axis_range=[],  out_file_name='' ):
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
        #gPad.SetLogz()

        root_is_dumb.append(hh)
        cc+=1

    can.SaveAs(out_file_name)
    can.Print(fname)


def make_plots_per_sector(can, cx, cy,  fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False, cut_pos=[],  out_file_name=''):
    can.Clear()
    can.Divide(cx,cy) #3,6)
    cc=1
    root_is_dumb = []
    ss=1
    for hh in hists:
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.Draw('colz')
        gPad.SetLogz()

        hh.SetLineColor(kBlack)
        hh.GetXaxis().SetTitleSize(0.06);
        hh.GetYaxis().SetTitleSize(0.06);
        hh.GetXaxis().SetLabelOffset(0.001);
        hh.GetXaxis().SetLabelSize(0.05)
        hh.GetYaxis().SetLabelSize(0.05)

        gStyle.SetOptStat(0)#"e")
        gStyle.SetStatY(0.85);      
        gStyle.SetStatX(0.8);  
        gStyle.SetStatW(0.1);                
        
        gPad.SetLogz()
        gPad.SetLeftMargin(0.15)
        gPad.SetBottomMargin(0.2)
        #gPad.SetRightMargin(0.15)
        can.Update()


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

    can.SaveAs(out_file_name)
    can.Print(fname)


def make_plots_per_sector_with_cut_lines(can, cx,cy,  fname, hists, hh_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False,cut_lines=[],  out_file_name=''):
    can.Clear()
    can.Divide(cx,cy)#3,6)
    cc=1
    root_is_dumb = []
    ss=0
    for hh in hists:
        print('drawing to canvas {}'.format(cc))
        can.cd(cc)
        hh.SetTitle(hh_titles[cc-1])
        hh.Draw('colz')
        gPad.SetLogz()

        hh.SetLineColor(kBlack)    
        hh.GetXaxis().SetTitleSize(0.06);
        hh.GetYaxis().SetTitleSize(0.06);
        hh.GetXaxis().SetLabelOffset(0.005);
        hh.GetXaxis().SetLabelSize(0.05)
        hh.GetYaxis().SetLabelSize(0.05)

        gStyle.SetOptStat(0)#"e")
        gStyle.SetStatY(0.85);      
        gStyle.SetStatX(0.8);  
        gStyle.SetStatW(0.1);                
        
        gPad.SetLogz()
        gPad.SetLeftMargin(0.15)
        gPad.SetBottomMargin(0.2)
        #gPad.SetRightMargin(0.15)
        can.Update()


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

    can.SaveAs(out_file_name)
    can.Print(fname)


#gStyle.SetOptStat(00000)
gStyle.SetPalette(55)

beamE = 10.604
hhs = {}
ff = TFile(sys.argv[1])
for kk in ff.GetListOfKeys():
    obj = kk.ReadObj()
    hhs[obj.GetName()] = obj

mon_out_file_name='monitor_eb_electron_pid_cuts.pdf'
can = TCanvas('can','can',4000,4000)
can.Print('{}['.format(mon_out_file_name))


#plot NPHE
hnphe_neg=hhs['nphe_cut0']
hnphe_neg_allbutnphe=hhs['nphe_pass_all_but_nphe']
hnphe_neg_pass=hhs['nphe_passall']
hnphe=[hnphe_neg,hnphe_neg_allbutnphe,hnphe_neg_pass]
hnphe_titles=["Nphe All Neg.; Nphe; counts",
              "Nphe All Cuts but Nphe; Nphe; counts",
              "Nphe Pass All Cuts; Nphe; counts"
              ]

make_plots(can,3,1,mon_out_file_name,hnphe,hnphe_titles,draw_cut=True, cut_pos=[2,2])

hnphe_sec = []
hnphe_sec_titles=[]
nphe_sec_cuts=[]
for ss in range(1,7):
    hnphe_sec.append(hhs['nphe_s'+str(ss)+'_cut0'])
    hnphe_sec_titles.append('Number of photo-electrons in HTCC Sector '+str(ss)+'; Nphe; counts')
    nphe_sec_cuts.append([2,2])

make_plots_per_sector(can,2,3,mon_out_file_name,hnphe_sec,hnphe_sec_titles,draw_cut=True, cut_pos=nphe_sec_cuts, out_file_name='nphe_per_sector_all_neg.pdf') 

#plot VZ
hvz_neg=[]
hvz_titles=[]
hvz_neg_set2=[]
hvz_titles_set2=[]
vz_cut=[]
for ss in range(1,7):
    vz_cut.append([-12,9])
        

for ss in range(1,7):
    hvz_titles.append('Vz All Neg. S'+ str(ss) +'; vz (cm); counts')
    hvz_neg.append(hhs['vz_cut0_s'+str(ss)])
    hvz_titles.append('Vz All Cuts but Vz S'+str(ss)+'; Vz (cm); counts')
    hvz_neg.append(hhs['vz_passall_passall_but_vz_s'+str(ss)])
    hvz_titles.append('Vz Pass All Cuts S'+str(ss)+';Vz (cm); counts')
    hvz_neg.append(hhs['vz_passall_s'+str(ss)])

for ss in range(1,7):
    hvz_titles_set2.append('Reconstructed Vz - All Cuts Applied but z-vertex ) Sector ' + str(ss) +'; z-vertex (cm); counts')
    hvz_neg_set2.append(hhs['vz_passall_passall_but_vz_s'+str(ss)])
    

make_plots_per_sector(can,3,6,mon_out_file_name,hvz_neg,hvz_titles,draw_cut=True, cut_pos=vz_cut )


make_plots_per_sector(can,2,3,mon_out_file_name,hvz_neg_set2,hvz_titles_set2,draw_cut=True, cut_pos=vz_cut, out_file_name='vz_per_sector_allbutvz.pdf')


make_plots_overlap1D(can, 900, 700, mon_out_file_name, hvz_neg_set2, hvz_titles_set2, draw_cut=True, cut_pos=vz_cut,out_file_name='vz_per_sector_allbutvz_overlap.pdf')



# plot EIEO
heieo=[]
heieo_titles=[]
eieo_cut=0.06
heieo=[hhs['eieo_small_cut0'], hhs['eieo_small_passall_but_eieo'], hhs['eieo_small_passall'] ]
heieo_titles=['EiEo All Neg; E_{inner} (GeV); E_{outer} (GeV)',
              'EiEo All Cuts but EiEo; E_{inner} (GeV); E_{outer} (GeV)',
              'EiEo Pass All Cuts; E_{inner} (GeV); E_{outer} (GeV)'
             ]

heieo_sect=[]
heieo_titles_sect=[]
eieo_cuts=[]
for ss in range(1,7):
    heieo_sect.append(hhs['eieo_small_cut0_s'+str(ss)])
    heieo_titles_sect.append('Inner vs Outer Energy Deposited Sector '+str(ss)+';E_{inner} (GeV); E_{outer} (GeV)')
    eieo_cuts.append([0.06,0.06])

make_plots(can,3,1,mon_out_file_name,heieo,heieo_titles,draw_cut=True, cut_pos=[eieo_cut,eieo_cut], min_y=0.0, max_y=2.0)

make_plots_per_sector(can,2,3,mon_out_file_name,heieo_sect,heieo_titles_sect,draw_cut=True, cut_pos=eieo_cuts, out_file_name='eieo_per_sector_allbuteieo.pdf')


hpcalusf = []
hpcalusf_titles = []
hpcalusf_set2 = []
hpcalusf_titles_set2 = []

hpcalvsf = []
hpcalvsf_titles = []
hpcalvsf_set2 = []
hpcalvsf_titles_set2 = []

hpcalwsf = []
hpcalwsf_titles = []
hpcalwsf_set2 = []
hpcalwsf_titles_set2 = []

hpcalu_cut=[ [33,408],
             [26,408],
             [34,408],
             [22,408],
             [27,408],
             [25,408]
         ]
            
#########################33
## PCAL U coordinate 

for ss in range(1,7):
    hpcalusf.append(hhs['pcalusf_cut0_s'+str(ss)])
    hpcalusf_titles.append('PCAL Sampling Fraction Layer U All Neg. S'+str(ss)+'; layer U; Sampling Fraction')
    hpcalusf.append(hhs['pcalusf_passall_but_pcalfid_s'+str(ss)])
    hpcalusf_titles.append('PCAL Sampling Fraction Layer U Pass All But PCAL Fid. S'+str(ss)+' ; layer U; Sampling Fraction')
    hpcalusf.append(hhs['pcalusf_passall_s'+str(ss)])
    hpcalusf_titles.append('PCAL Sampling Fraction Layer U Pass All S'+str(ss)+'; layer U; Sampling Fraction')

    hpcalusf_set2.append(hhs['pcalusf_passall_but_pcalfid_s'+str(ss)])
    hpcalusf_titles_set2.append('PCAL Sampling Fraction Layer U Pass All But PCAL Fid. S'+str(ss)+' ; layer U; Sampling Fraction')

make_plots_per_sector(can,3,6,mon_out_file_name,hpcalusf,hpcalusf_titles,draw_cut=True, cut_pos=hpcalu_cut, min_y=0, max_y=0.3)#

make_plots_per_sector(can,2,3,mon_out_file_name,hpcalusf_set2,hpcalusf_titles_set2,draw_cut=True, cut_pos=hpcalu_cut, min_y=0, max_y=0.3, out_file_name='pcalusf_per_sector_allbutusf.pdf')


#########################33
## PCAL V coordinate 

hpcalv_cut=[ [16,400],
             [10.5,400],
             [17.0,400],
             [14.25,400],
             [18.0,400],
             [11,400]
         ]
            
for ss in range(1,7):
    hpcalvsf.append(hhs['pcalvsf_cut0_s'+str(ss)])
    hpcalvsf_titles.append('PCAL Sampling Fraction Layer V All Neg. S'+str(ss)+'; layer V; Sampling Fraction')
    hpcalvsf.append(hhs['pcalvsf_passall_but_pcalfid_s'+str(ss)])
    hpcalvsf_titles.append('PCAL Sampling Fraction Layer V Pass All But PCAL Fid. S'+str(ss)+' ; layer V; Sampling Fraction')
    hpcalvsf.append(hhs['pcalvsf_passall_s'+str(ss)])
    hpcalvsf_titles.append('PCAL Sampling Fraction Layer V Pass All S'+str(ss)+'; layer V; Sampling Fraction')

    hpcalvsf_set2.append(hhs['pcalvsf_passall_but_pcalfid_s'+str(ss)])
    hpcalvsf_titles_set2.append('PCAL Sampling Fraction Layer V Pass All But PCAL Fid. S'+str(ss)+' ; layer V; Sampling Fraction')


make_plots_per_sector(can,3,6,mon_out_file_name,hpcalvsf,hpcalvsf_titles, draw_cut=True, cut_pos=hpcalv_cut, min_y=0, max_y=0.3)#out_file_name='vsf_per_sector_allbutvsf.pdf'))

make_plots_per_sector(can,2,3,mon_out_file_name,hpcalvsf_set2,hpcalvsf_titles_set2,draw_cut=True, cut_pos=hpcalv_cut, min_y=0, max_y=0.3, out_file_name='pcalvsf_per_sector_allbutvsf.pdf')

#########################33
## PCAL W coordinate 
hpcalw_cut=[ [11, 400],
             [17.5, 400],
             [16.25, 400],
             [7.25, 400],
             [14.5, 400],
             [9.25, 400]
         ]
            
for ss in range(1,7):
    hpcalwsf.append(hhs['pcalwsf_cut0_s'+str(ss)])
    hpcalwsf_titles.append('PCAL Sampling Fraction Layer W All Neg. S'+str(ss)+'; layer W; Sampling Fraction')
    hpcalwsf.append(hhs['pcalwsf_passall_but_pcalfid_s'+str(ss)])
    hpcalwsf_titles.append('PCAL Sampling Fraction Layer W Pass All But PCAL Fid. S'+str(ss)+' ; layer W; Sampling Fraction')
    hpcalwsf.append(hhs['pcalwsf_passall_s'+str(ss)])
    hpcalwsf_titles.append('PCAL Sampling Fraction Layer W Pass All S'+str(ss)+'; layer W; Sampling Fraction')

    hpcalwsf_set2.append(hhs['pcalwsf_passall_but_pcalfid_s'+str(ss)])
    hpcalwsf_titles_set2.append('PCAL Sampling Fraction Layer W Pass All But PCAL Fid. S'+str(ss)+' ; layer W; Sampling Fraction')

make_plots_per_sector(can,3,6,mon_out_file_name,hpcalwsf,hpcalwsf_titles, draw_cut=True, cut_pos=hpcalw_cut, min_y=0, max_y=0.3)

make_plots_per_sector(can,2,3,mon_out_file_name,hpcalwsf_set2,hpcalwsf_titles_set2,draw_cut=True, cut_pos=hpcalw_cut, min_y=0, max_y=0.3, out_file_name='pcalwsf_per_sector_allbutwsf.pdf')

#########################################
## Drift Chambers Region 1
hdcr1=[]
hdcr1_titles=[]


hdcr1=[hhs['dcr1hit_cut0'], hhs['dcr1hit_pass_all_but_dcr1cut'], hhs['dcr1hit_passall'] ]
hdcr1_titles=['DCR1 xy All Neg; x (cm); y (cm)',
              'DCR1 xy All Cuts but DCR1 Fid.; x (cm); y (cm)',
              'DCR1 Pass All Cuts; x (cm); y (cm)'
          ]

hdcr1_set2=[hhs['dcr1hit_cut0'],  hhs['dcr1hit_passall'] ]
hdcr1_titles_set2=['DCR1 xy All Neg; x (cm); y (cm)',
              'DCR1 Pass All Cuts; x (cm); y (cm)'
          ]

can.SetCanvasSize(1800,600)
make_plots_overlap(can,3,1,mon_out_file_name, hdcr1, hdcr1_titles, cut_before=hdcr1[0] )

can.SetCanvasSize(1000,1000)
make_plots_overlap(can,1,1,mon_out_file_name, hdcr1_set2, hdcr1_titles_set2, cut_before=hdcr1[0], out_file_name='dcr1_overlap.pdf' )

make_plots(can,3,1,mon_out_file_name,hdcr1,hdcr1_titles,draw_cut=False)

#########################################
## Drift Chambers Region 2
hdcr2=[]
hdcr2_titles=[]

hdcr2=[hhs['dcr2hit_cut0'], hhs['dcr2hit_pass_all_but_dcr2cut'], hhs['dcr2hit_passall'] ]
hdcr2_titles=['DCR2 xy All Neg; x (cm); y (cm)',
              'DCR2 xy All Cuts but DCR2 Fid.; x (cm); y (cm)',
              'DCR2 Pass All Cuts; x (cm); y (cm)'
          ]

can.SetCanvasSize(1800,600)
make_plots_overlap(can,3,1,mon_out_file_name, hdcr2, hdcr2_titles, cut_before=hdcr2[0] )
make_plots(can,3,1,mon_out_file_name,hdcr2,hdcr2_titles,draw_cut=False)

#########################################
## Drift Chambers Region 3
hdcr3=[]
hdcr3_titles=[]

hdcr3=[hhs['dcr3hit_cut0'], hhs['dcr3hit_pass_all_but_dcr3cut'], hhs['dcr3hit_passall'] ]
hdcr3_titles=['DCR3 xy All Neg; x (cm); y (cm)',
              'DCR3 xy All Cuts but DCR3 Fid.; x (cm); y (cm)',
              'DCR3 Pass All Cuts; x (cm); y (cm)'
          ]

can.SetCanvasSize(1800,600)
make_plots_overlap(can,3,1,mon_out_file_name, hdcr3, hdcr3_titles, cut_before=hdcr3[0] )
make_plots(can,3,1,mon_out_file_name,hdcr3,hdcr3_titles,draw_cut=False)


hsf = []
hsf_titles=[]
hsf_set2 = []
hsf_titles_set2=[]
sigma_range = 3.0
p0mean_inb = [0.105631, 0.11551, 0.112799, 0.109937, 0.116249, 0.119057]
p1mean_inb = [-0.153951, -0.0253273, -0.125718, 0.165414, 0.0768411, 0.0555026]
p2mean_inb = [0.00860091, 0.00706291, 0.00908884, 0.00499666, 0.00448701, 0.00558927]
p3mean_inb = [-0.000821675, -0.000711488, -0.000930922, -0.000298311, -0.000455716, -0.000657084]
p0sigma_inb = [0.0149613, 0.0115116, 0.00580737, 0.0106817, 0.012667, 0.00553471]
p1sigma_inb = [0.00700773, 0.0116193, 0.0202375, 0.0126958, 0.00892239, 0.0216206]

hsfchi=[]
hsfchi_titles=[]


ecsf_cuts=[]

def returnSFCut(p0,p1,p2,p3,p0s,p1s,sect,sig,limit):
    cut_line = TF1('ecsf_'+limit+'_s'+str(sect),'([0]*(1 + x/TMath::Sqrt(x*x + [1])) + x*[2] + x*x*[3]) + [6]*([4] + [5]/TMath::Sqrt(x))',1.4,10.0)
    cut_line.SetParameter(0,p0)
    cut_line.SetParameter(1,p1)
    cut_line.SetParameter(2,p2)
    cut_line.SetParameter(3,p3)
    cut_line.SetParameter(4,p0s)
    cut_line.SetParameter(5,p1s)
    cut_line.SetParameter(6,sig)
    return cut_line

    

sf_sig=3
for ss in range(0,6):    
    print('making plots for sector {}'.format(ss))
    hsf.append(hhs['calsf_cut0_s'+str(ss)])
    hsf.append(hhs['calsf_passall_but_ecsf_s'+str(ss)])
    hsf.append(hhs['calsf_passall_s'+str(ss)])
    hsf_titles.append('Cal. Sampling Fraction All Neg. S'+str(ss+1)+'; p (GeV); Edep/p ')
    hsf_titles.append('Cal. Sampling Fraction Pass All Cuts But SF S'+str(ss+1)+'; p (GeV); Edep/p ')
    hsf_titles.append('Cal. Sampling Fraction Pass All Cuts S'+str(ss+1)+'; p (GeV); Edep/p ')

    sect=ss
    cut_min =returnSFCut(p0mean_inb[sect],p1mean_inb[sect],p2mean_inb[sect],p3mean_inb[sect],p0sigma_inb[sect],p1sigma_inb[sect],sect,-sf_sig,'min')
    cut_max =returnSFCut(p0mean_inb[sect],p1mean_inb[sect],p2mean_inb[sect],p3mean_inb[sect],p0sigma_inb[sect],p1sigma_inb[sect],sect,sf_sig,'max')
    cut_pair=[]
    cut_pair.append(cut_min)
    cut_pair.append(cut_max)
    ecsf_cuts.append(cut_pair)

    hsfchi.append(hhs['calsfchi_passall_s'+str(ss)])
    hsfchi_titles.append('Cal. Sampling Fraction Pull Dist. Pass All Cuts Sector '+str(ss+1)+'; pull; counts')

for ss in range(0,6):
    hsf_set2.append(hhs['calsf_passall_but_ecsf_s'+str(ss)])
    hsf_titles_set2.append('Cal. Sampling Fraction Pass All Cuts But SF S'+str(ss+1)+'; p (GeV); Edep/p ')

can.SetCanvasSize(6000,6000)
make_plots_per_sector_with_cut_lines(can,3,6, mon_out_file_name, hsf, hsf_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=True,cut_lines=ecsf_cuts)

can.SetCanvasSize(6000,6000)
make_plots_per_sector_with_cut_lines(can,2,3, mon_out_file_name, hsf_set2, hsf_titles_set2, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=True,cut_lines=ecsf_cuts, out_file_name='ecsf_per_sector_allbutecsf.pdf')

can.SetCanvasSize(6000,6000)
make_plots_per_sector_with_cut_lines(can,2,3, mon_out_file_name, hsfchi, hsfchi_titles, x_ranges=[], rebinx=1, min_y=0, max_y=0, draw_cut=False,cut_lines=[], out_file_name='ecsfchi_per_sector_passall.pdf')


can.SetCanvasSize(6000,6000)
make_plots(can,1,1,mon_out_file_name,[hhs['iron_cross_xscint_ypcal_passall']],['Cross FTOF vs PCAL; #Delta t FTOF; #Delta t PCAL'],draw_cut=False, cut_pos=[])                                                                             

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



can.Print('{}]'.format(mon_out_file_name))
