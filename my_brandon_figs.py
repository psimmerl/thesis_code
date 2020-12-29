import ROOT, math, sys, argparse
from MyHist import MyHist

ct = "counts"
xvz = "#Delta Vz_{K^{+}}-Vz_{K^{-}}"
xim2, yim2 = "I.M. PrK^{-}", "I.M. K^{+}K^{-}"
xim = ["I.M. PrK^{+}","I.M. PrK^{-}","I.M. K^{+}K^{-}"]
xc = ["Coplanarity "+s for s in ["#Theta_{pr}","#Theta_{K^{+}}","#Theta_{K^{-}}"]]
ep,MM2=": ep#rightarrow e","Missing Mass^{2}"

cg, cb = "PASS_","FAIL_"

pDir = "figs/"

mmcs = ['pass_all_but_missing_energy',
    'fail_at_least_missing_energy',
    'pass_all_but_missing_proton',
    'fail_at_least_missing_proton',
    'pass_all_but_missing_kaonP',
    'fail_at_least_missing_kaonP',
    'pass_all_but_missing_kaonM',
    'fail_at_least_missing_kaonM']
csp =  ['combssize1',
    'pass_all',
    'pass_all_pass_cpl_pro',
    'pass_all_pass_cpl_kp',
    'pass_all_pass_cpl_km',
    'pass_all_pass_cpl_all',
    'pass_all_pass_dvz_kpkm',
    'pass_all_pass_additional',
    'pass_all_pass_additional_pass_phi_mass']
csf =  ['combssize1',
    'fail_all',
    'pass_all_fail_cpl_pro',
    'pass_all_fail_cpl_kp',
    'pass_all_fail_cpl_km',
    'pass_all_fail_cpl_all',
    'pass_all_fail_dvz_kpkm',
    'pass_all_fail_additional',
    'pass_all_pass_additional_fail_phi_mass']
#'p_ele_', 'p_pro_', 'p_kp_', 'p_km_', 
#'phi_ele_', 'phi_pro_', 'phi_kp_', 'phi_km_', 
#'phi_vs_phi_', 'q2_','w_','xb_','t_', 't_q2_','t_xb_','t_w',
#,'mm_ekpX_'
pqwxb = ['p_ele_theta_', 'phi_ele_theta_', 'p_pro_theta_', 'phi_pro_theta_', 
    'p_kp_theta_', 'phi_kp_theta_','p_km_theta_', 'phi_km_theta_','w_q2_','xb_q2_']
cpl = ['cpl_pro_','cpl_kp_','cpl_km_']
mm = ['mm_epkpkmX_','mm_ekpkmX_','mm_epkpX_','mm_epkmX_','missing_energy_']
im = ['im_pkm_','im_pkp_','im_kpkm_', 'im_kpkm_im_pkm_']
dvz = ['dvz_kpkm_']

ffi = "monitor_phi_skim8_005032_pidtype1_cutvar1_rmvres12-inb_16core_incl.root"

met   = MyHist(mm[4]+csp[0],    ffi, "Missing Energy"+ep+"pK^{+}K^{-}X", "Energy (GeV)", ct, "mm")
mep   = MyHist(mm[4]+mmcs[0],    ffi, "Missing Energy"+ep+"pK^{+}K^{-}X", "Energy (GeV)", ct, "mm")
mef   = MyHist(mm[4]+mmcs[1],    ffi, "Missing Energy"+ep+"pK^{+}K^{-}X", "Energy (GeV)", ct, "mm")

prot  = MyHist(mm[1]+csp[0], ffi, MM2+ep+"K^{+}K^{-}X", "Proton Mass (GeV^{2})", ct, "mm")
prop  = MyHist(mm[1]+mmcs[2], ffi, MM2+ep+"k^{+}K^{-}X", "Proton Mass (GeV^{2})", ct, "mm")
prof  = MyHist(mm[1]+mmcs[3], ffi, MM2+ep+"K^{+}K^{-}X", "Proton Mass (GeV^{2})", ct, "mm")

kpt   = MyHist(mm[3]+csp[0], ffi, MM2+ep+"pK^{-}X", "K^{+} Mass (GeV^{2})", ct, "mm")
kpp   = MyHist(mm[3]+mmcs[4], ffi, MM2+ep+"pK^{-}X", "K^{+} Mass (GeV^{2})", ct, "mm")
kpf   = MyHist(mm[3]+mmcs[5], ffi, MM2+ep+"pK^{-}X", "K^{+} Mass (GeV^{2})", ct, "mm")

kmt   = MyHist(mm[2]+csp[0], ffi, MM2+ep+"pK^{+}X", "K^{-} Mass (GeV^{2})", ct, "mm")
kmp   = MyHist(mm[2]+mmcs[6], ffi, MM2+ep+"pK^{+}X", "K^{-} Mass (GeV^{2})", ct, "mm")
kmf   = MyHist(mm[2]+mmcs[7], ffi, MM2+ep+"pK^{+}X", "K^{-} Mass (GeV^{2})", ct, "mm")
    

vz  = MyHist(dvz[0]+csp[1], ffi, "(FD), "+xvz+"Pass", xvz+" (cm)", ct,"dvz")
h34br= MyHist(im[2]+csp[1], ffi, xim[2]+" Pass-Fail "+xvz, xim[2]+" (GeV)", ct,"imk")
h34bp= MyHist(im[2]+csp[6], ffi, xim[2]+" Pass-Fail "+xvz, xim[2]+" (GeV)", ct,"imk")
h34bf= MyHist(im[2]+csf[6], ffi, xim[2]+" Pass-Fail "+xvz, xim[2]+" (GeV)", ct,"imk")
    
pac, pfc, l9 = ", Pass All Cuts", xim[2]+", Pass-Fail Coplan. #theta_{","<9 deg Cuts"
cp = MyHist(cpl[0]+csp[1],ffi, "Proton (FD), "+xc[0]+pac, xc[0], ct,"coplane")
h35bt= MyHist(im[2]+csp[1], ffi, pfc+"pr}"+l9, xim[2]+" (GeV)",       ct,"imk")
h35bp= MyHist(im[2]+csp[2], ffi, pfc+"pr}"+l9, xim[2]+" (GeV)",       ct,"imk")
h35bf= MyHist(im[2]+csf[2], ffi, pfc+"pr}"+l9, xim[2]+" (GeV)",       ct,"imk")
ckp = MyHist(cpl[1]+csp[1], ffi, "K^{+} (FD), "+xc[1]+pac, xc[1], ct,"coplane")
h35dt= MyHist(im[2]+csp[1], ffi, pfc+"K^{+}}"+l9, xim[2]+" (GeV)",       ct,"imk")
h35dp= MyHist(im[2]+csp[3], ffi, pfc+"K^{+}}"+l9, xim[2]+" (GeV)",       ct,"imk")
h35df= MyHist(im[2]+csf[3], ffi, pfc+"K^{+}}"+l9, xim[2]+" (GeV)",       ct,"imk")
ckm = MyHist(cpl[2]+csp[1], ffi, "K^{-} (FD), "+xc[2]+pac, xc[2], ct,"coplane")
h35ft= MyHist(im[2]+csp[1], ffi, pfc+"K^{-}}"+l9, xim[2]+" (GeV)",       ct,"imk")
h35fp= MyHist(im[2]+csp[4], ffi, pfc+"K^{-}}"+l9, xim[2]+" (GeV)",       ct,"imk")
h35ff= MyHist(im[2]+csf[4], ffi, pfc+"K^{-}}"+l9, xim[2]+" (GeV)",       ct,"imk")
    
h36r= MyHist(im[2]+csp[1], ffi, xim[2]+", Pass-Fail All Coplan. Cuts", xim[2]+" (GeV)", ct,"imk")
h36p= MyHist(im[2]+csp[5], ffi, xim[2]+", Pass-Fail All Coplan. Cuts", xim[2]+" (GeV)", ct,"imk")
h36f= MyHist(im[2]+csf[5], ffi, xim[2]+", Pass-Fail All Coplan. Cuts", xim[2]+" (GeV)", ct,"imk")
    

im2 = MyHist(im[3]+csp[5], ffi, "(FD), "+yim2+" vs "+xim2+", Pass All Cuts", \
    xim2+" (GeV)", yim2+" (GeV)","im2D")


ffo = "monitor_phi_epKpKm_RGA_OUT_pidtype1_cutvar1_rmvres12-inb_16core_incl.root"

h37i = MyHist(im[0]+csp[8], ffi, "(FD) inbNoutb, "+xim[1]+", Pass All,--<I.M.K^{+}K^{-}<--",\
          xim[1]+" (GeV)", ct,"imk")


h37o = MyHist(im[0]+csp[8], ffo, "(FD) outb, "+xim[1]+", Pass All,--<I.M.K^{+}K^{-}<--",\
          xim[1]+" (GeV)", ct,"imk")
h37i.Add(h37o)


h38i = MyHist(im[2]+csp[7], ffi, "Invariant Mass of Charged Kaons from epK^{+}K^{-}", \
    xim[2]+" (GeV)", ct+" inbNoutb","imk")
h38o = MyHist(im[2]+csp[7], ffo, "Invariant Mass of Charged Kaons from epK^{+}K^{-}", \
    xim[2]+" (GeV)", ct+" inbNoutb","imk")
h38i.Add(h38o)


met.setRange(  [-1.5,1.5] );mep.setRange(  [-1.5,1.5] );mef.setRange(  [-1.5,1.5] )
prot.setRange( [0,2] );prop.setRange( [0,2] );prof.setRange( [0,2] )
kpt.setRange( [-.5,1.5] );kpp.setRange( [-.5,1.5] );kpf.setRange( [-.5,1.5] )
kmt.setRange( [-.5,1.5] );kmp.setRange( [-.5,1.5] );kmf.setRange( [-.5,1.5] )

met.setCanvas(True,2,2)
met.Draw( "same",pDir+"32b.png",1,ROOT.kBlack,print_fit=False)
mep.Draw( "same",pDir+"32b.png",1,ROOT.kBlue,print_fit=False)
mef.Draw( "same",pDir+"32b.png",1,ROOT.kRed,print_fit=False)

prot.Draw("same",pDir+"32b.png",2,ROOT.kBlack,print_fit=False)
prop.Draw("same",pDir+"32b.png",2,ROOT.kBlue,print_fit=False)
prof.Draw("same",pDir+"32b.png",2,ROOT.kRed,print_fit=False)

kpt.Draw( "same",pDir+"32b.png",4,ROOT.kBlack,print_fit=False)
kpp.Draw( "same",pDir+"32b.png",4,ROOT.kBlue,print_fit=False)
kpf.Draw( "same",pDir+"32b.png",4,ROOT.kRed,print_fit=False)

kmt.Draw( "same",pDir+"32b.png",3,ROOT.kBlack,print_fit=False)
kmp.Draw( "same",pDir+"32b.png",3,ROOT.kBlue,print_fit=False)
kmf.Draw( "same",pDir+"32b.png",3,ROOT.kRed,print_fit=False)

ppars=prot.gaussFit( 0.938 )
mepars=met.gaussFit( 0.0 )
kppars=kpt.gaussFit( 0.494**2 )
kmpars=kmt.gaussFit( 0.494**2 )

met.setRange(  [-1.5,1.5] );mep.setRange(  [-1.5,1.5] );mef.setRange(  [-1.5,1.5] )
prot.setRange( [0,2] );prop.setRange( [0,2] );prof.setRange( [0,2] )
kpt.setRange( [-.5,1.5] );kpp.setRange( [-.5,1.5] );kpf.setRange( [-.5,1.5] )
kmt.setRange( [-.5,1.5] );kmp.setRange( [-.5,1.5] );kmf.setRange( [-.5,1.5] )

met.setCanvas(True,2,2)
met.Draw( "same",pDir+"32a.png",1, lines=[mepars[1]+4*mepars[2],mepars[1]-4*mepars[2]])
prot.Draw("same",pDir+"32a.png",2, lines=[ppars[1]+4*ppars[2],ppars[1]-4*ppars[2]])
kpt.Draw( "same",pDir+"32a.png",4, lines=[kppars[1]+4*kppars[2],kppars[1]-4*kppars[2]])
kmt.Draw( "same",pDir+"32a.png",3, lines=[kmpars[1]+4*kmpars[2],kmpars[1]-4*kmpars[2]])

im2.setCanvas(True); im2.Draw("COLZ", pDir+"33.png")

vzpars = vz.gaussFit( 0.0 )
vz.setCanvas(True);
vz.Draw("same", pDir+"34a.png",lines=[vzpars[1]+4*vzpars[2],vzpars[1]-4*vzpars[2]])

h34br.setCanvas(True,1,1) 
h34br.Draw( "same",pDir+"34b.png",1,ROOT.kBlack)
h34bp.Draw( "same",pDir+"34b.png",1,ROOT.kBlue)
h34bf.Draw( "same",pDir+"34b.png",1,ROOT.kRed)


cp.setCanvas(True); cp.Draw( "",pDir+"35a.png",lines=[9])
ckp.setCanvas(True); ckp.Draw( "",pDir+"35c.png",lines=[9])
ckm.setCanvas(True); ckm.Draw( "",pDir+"35e.png",lines=[9])

h35bt.setCanvas(True)
h35bt.Draw("same",pDir+"35b.png",1,ROOT.kBlack)
h35bp.Draw("same",pDir+"35b.png",1,ROOT.kBlue)
h35bf.Draw("same",pDir+"35b.png",1,ROOT.kRed)

h35dt.setCanvas(True)
h35dt.Draw("same",pDir+"35d.png",1,ROOT.kBlack)
h35dp.Draw("same",pDir+"35d.png",1,ROOT.kBlue)
h35df.Draw("same",pDir+"35d.png",1,ROOT.kRed)

h35ft.setCanvas(True)
h35ft.Draw("same",pDir+"35f.png",1,ROOT.kBlack)
h35fp.Draw("same",pDir+"35f.png",1,ROOT.kBlue)
h35ff.Draw("same",pDir+"35f.png",1,ROOT.kRed)

h36r.setCanvas(True)
h36r.Draw("same",pDir+"36.png",1,ROOT.kBlack)
h36p.Draw("same",pDir+"36.png",1,ROOT.kBlue)
h36f.Draw("same",pDir+"36.png",1,ROOT.kRed)

h37i.setCanvas(True); h37i.Draw("",pDir+"37.png")
h38i.gaussFit( 1.02 )
h38i.setRange([0.96,1.14])
h38i.setCanvas(True); h38i.Draw("",pDir+"38.png")



#q2w,q2xb,phip,thetap


q2wti  = MyHist('w_q2_'+csp[0], ffi, "Q^{2} vs W, cs1 inb",     "W (GeV)", "Q^{2} (GeV^{2})", "q2w")
q2Xbti = MyHist('xb_q2_'+csp[0],ffi, "Q^{2} vs x_{B}, cs1 inb", "x_{B}",   "Q^{2} (GeV^{2})","q2Xb")
q2wci  = MyHist('w_q2_'+csp[8], ffi, "Q^{2} vs W, all inb",     "W (GeV)", "Q^{2} (GeV^{2})", "q2w")
q2Xbci = MyHist('xb_q2_'+csp[8],ffi, "Q^{2} vs x_{B}, all inb", "x_{B}",   "Q^{2} (GeV^{2})","q2Xb")

q2wto  = MyHist('w_q2_'+csp[0], ffo, "Q^{2} vs W, cs1 outb",     "W (GeV)", "Q^{2} (GeV^{2})", "q2w")
q2Xbto = MyHist('xb_q2_'+csp[0],ffo, "Q^{2} vs x_{B}, cs1 outb", "x_{B}",   "Q^{2} (GeV^{2})","q2Xb")
q2wco  = MyHist('w_q2_'+csp[8], ffo, "Q^{2} vs W, all outb",     "W (GeV)", "Q^{2} (GeV^{2})", "q2w")
q2Xbco = MyHist('xb_q2_'+csp[8],ffo, "Q^{2} vs x_{B}, all outb", "x_{B}",   "Q^{2} (GeV^{2})","q2Xb")


for p, P in zip([0, 2, 4, 6],[s+", "for s in ["Electron","Proton","K^{+}","K^{-}"]]):
    ThPti = MyHist(pqwxb[p]+csp[0], ffi, P+"#Theta vs P, cs1 inb","P (GeV)","#Theta","ThP")
    PhiPti= MyHist(pqwxb[p+1]+csp[0], ffi, P+"#Theta vs #Phi, cs1 inb","#Phi","#Theta","ThPhi")
    ThPci = MyHist(pqwxb[p]+csp[8], ffi, P+"#Theta vs P, all inb","P (GeV)","#Theta","ThP")
    PhiPci= MyHist(pqwxb[p+1]+csp[8], ffi, P+"#Theta vs #Phi, all inb","#Phi","#Theta","ThPhi")
    ThPti.setCanvas(True,2,2); pp = [ThPti,PhiPti,ThPci,PhiPci]
    for t in range(4): pp[t].Draw("COLZ",pDir+"phase_inb_"+P[0:-2]+".png",t+1)

    ThPto = MyHist(pqwxb[p]+csp[0], ffo, P+"#Theta vs P, cs1 outb","P (GeV)","#Theta","ThP")
    PhiPto= MyHist(pqwxb[p+1]+csp[0], ffo, P+"#Theta vs #Phi, cs1 outb","#Phi","#Theta","ThPhi")
    ThPco = MyHist(pqwxb[p]+csp[8], ffo, P+"#Theta vs P, all outb","P (GeV)","#Theta","ThP")
    PhiPco= MyHist(pqwxb[p+1]+csp[8], ffo, P+"#Theta vs #Phi, all outb","#Phi","#Theta","ThPhi")
    ThPto.setCanvas(True,2,2); pp = [ThPto,PhiPto,ThPco,PhiPco]
    for t in range(4): pp[t].Draw("COLZ",pDir+"phase_outb_"+P[0:-2]+".png",t+1)
    


q2 = [q2wti, q2Xbti, q2wci, q2Xbci, q2wto, q2Xbto, q2wco, q2Xbco]

for i in range(8):
    if i%4==0: q2[i].setCanvas(True, 2, 2)
    bt = "inb" if i < 4 else "outb"
    q2[i].Draw("COLZ",pDir+"q2_"+bt+".png", (i%4) + 1)


