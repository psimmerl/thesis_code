#!/usr/bin/python

import sys, ROOT
from ROOT import (kTRUE, gROOT)

gROOT.SetBatch(kTRUE);

mmcs = ['pass_all_but_missing_energy',
        'fail_at_least_missing_energy',
        'pass_all_but_missing_proton',
        'fail_at_least_missing_proton',
        'pass_all_but_missing_kaonP',
        'fail_at_least_missing_kaonP',
        'pass_all_but_missing_kaonM',
        'fail_at_least_missing_kaonM']
cs =  ['combssize1',
  'pass_all','fail_all']\
  +['pass_all'+cc for cc in [
  '_pass_cpl_pro',                 '_fail_cpl_pro',
  '_pass_cpl_kp',                  '_fail_cpl_kp',
  '_pass_cpl_km',                  '_fail_cpl_km',
  '_pass_cpl_all',                 '_fail_cpl_all',
  '_pass_dvz_kpkm',                '_fail_dvz_kpkm',
  '_pass_additional',              '_fail_additional',
  '_pass_additional_pass_phi_mass','_pass_additional_fail_phi_mass']]

#'p_ele_', 'p_pro_', 'p_kp_', 'p_km_', 
#'phi_ele_', 'phi_pro_', 'phi_kp_', 'phi_km_', 
#'phi_vs_phi_', 'q2_','w_','xb_','t_', 't_q2_','t_xb_','t_w',
#,'mm_ekpX_'
phasewxb = ['p_ele_theta_', 'phi_ele_theta_', 
            'p_pro_theta_', 'phi_pro_theta_', 
            'p_kp_theta_', 'phi_kp_theta_',
            'p_km_theta_', 'phi_km_theta_',
            'w_q2_','xb_q2_']
cpl = ['cpl_pro_','cpl_kp_','cpl_km_']
mm = ['mm_epkpkmX_','mm_ekpkmX_','mm_epkpX_','mm_epkmX_','missing_energy_']
im = ['im_pkm_','im_pkp_','im_kpkm_', 'im_kpkm_im_pkm_']
dvz = ['dvz_kpkm_']

hists = [ct + c for ct in phasewxb+cpl+mm+im+dvz for c in cs+mmcs]


def getObjects(directory, objs={}):
  for kk in directory.GetListOfKeys():
    obj = kk.ReadObj()
    try:
      getObjects(obj, objs)
    except:
      path = directory.GetPath().split(':')[-1].strip('/')
      objs[path+'/'+kk.GetName()] = obj
  return objs


ff = ROOT.TFile(sys.argv[1])
c1 = ROOT.TCanvas("c1","c1",1100,800)
c1.SetBatch(kTRUE)
pdfname = sys.argv[1].split('/')[-1]+'.pdf'
c1.Print(pdfname+'[')
objs = getObjects(ff)
for kk in sorted(objs):
  
  #if 'phi_trento_vs_phi_ele_pass_all' not in objs[kk].GetTitle(): continue
  if all([hist != objs[kk].GetTitle() for hist in hists]):
    print("Skipping")
    continue
  
  opt = ''
  if 'TH2' in objs[kk].__class__.__name__:
    opt = 'colz'
  if 'TGraph' in objs[kk].__class__.__name__:
    opt = 'AP'
  objs[kk].Draw(opt)
  c1.Print(pdfname)
c1.Print(pdfname+']')
  
