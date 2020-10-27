#!/opt/coatjava/6.5.3/bin/run-groovy
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import org.jlab.groot.data.TDirectory
import org.jlab.jroot.ROOTFile
import org.jlab.jroot.TNtuple
import java.io.File

import pid.electron.ElectronFromEvent
import pid.hadron.HadronID
import examples.epkpkm
import examples.epkpX
import event.EventConverter
import event.Event


LorentzVector.metaClass.static.withPID = {pid,x,y,z->
    double m = PDGDatabase.getParticleMass(pid)
    double e = Math.sqrt(x*x+y*y+z*z+m*m)
    return new LorentzVector(x,y,z,e)
}

LorentzVector.metaClass.plus = {
    LorentzVector a = new LorentzVector(delegate)
    a.add(it)
    return a
}

LorentzVector.metaClass.minus = {
    LorentzVector a = new LorentzVector(delegate)
    a.sub(it)
    return a
}

LorentzVector.metaClass.sector = null

Vector3.metaClass.plus = {
    Vector3 a = new Vector3(delegate)
    a.add(it)
    return a
}

Vector3.metaClass.minus = {
    Vector3 a = new Vector3(delegate)
    a.sub(it)
    return a
}


def hists = new ConcurrentHashMap()
def beam = LorentzVector.withPID(11,0,0,10.604)
def target = LorentzVector.withPID(2212,0,0,0)
def sig = 5.0


// load cut function for mm2 cuts
def loadCuts(file, left_sig_range, right_sig_range){
    cuts=[ ]
    if( file.exists() ){
	file.eachLine {  
	    line -> println "line : $line"; 
	    def mean = Double.parseDouble(line.split(" ")[0])
	    def sig = Double.parseDouble(line.split(" ")[1])
	    println("--- Cut Parameters: " + mean + " " + sig)	    
	    def lims = [mean - left_sig_range*sig, mean + right_sig_range*sig]
	    cuts = lims	
	}
    }
    return cuts
}

// cut base path
def base_path_cut = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/analysis_code_fork0/projects/exclusive_phi/epkpkm_top/"
def epkpkmxe_cut = loadCuts(new File(base_path_cut+"excl_phi_me_limits_pidtype1.txt"),4,4)
def ekpkmX_cut = loadCuts(new File(base_path_cut+"excl_phi_mm2_pro_limits_pidtype1.txt"),4,4)
def epkmX_cut = loadCuts(new File(base_path_cut+"excl_phi_mm2_kp_limits_pidtype1.txt"),4,4)
def epkpX_cut = loadCuts(new File(base_path_cut+"excl_phi_mm2_km_limits_pidtype1.txt"),4,4)

println(' Event Cut parameters ')
println(' --> me ' + epkpkmxe_cut)
println(' --> mm2 pr ' + ekpkmX_cut)
println(' --> mm2 kp ' + epkmX_cut)
println(' --> mm2 km ' + epkpX_cut)


excl_cuts = [
    epkpkmxe  : epkpkmxe_cut, //[0.08 - 0.0398*sig, 0.08 + 0.0398*sig],    
    epkpX : epkpX_cut,        // [0.248 - 0.059*sig, 0.248 + 0.059*sig],    
    epkmX : epkmX_cut,        //[0.248 - 0.055*sig, 0.248 + 0.055*sig],
    ekpkmX : ekpkmX_cut,      //[0.9431 - 0.0719*sig, 0.9431 + 0.0719*sig],
    el_p_min : 0.5,
    pro_p_max : 3.5,
    kp_p_max : 3.5,
    km_p_max : 3.5,
    pt_max : 0.12,
    q2_min : 1,
    w_min : 2
    
]

println(' old cut values me ' + excl_cuts.epkpkmxe)
println(' old cut values mm2 pro ' + excl_cuts.ekpkmX)
println(' old cut values mm2 kp ' + excl_cuts.epkmX)
println(' old cut values mm2 km ' + excl_cuts.epkpX)

tighter_kin_bounds = [
    theta_ele : [0, 45],
    theta_ele_ftof : [5, 30],
    theta_pro : [0, 70],
    theta_pro_ftof : [0, 60],
    theta_sum : [0, 120],
    p_ele     : [0, 10.5],
    p_ele_ftof : [0.0, 10.5],
    p_pro     : [0.0, 5.5],
    p_pro_ftof : [0.0, 6.0],
    phi_ele : [-180, 180],
    phi_pro : [-180, 180],

    p_kp_ftof : [0, 10],
    theta_kp_ftof : [0, 60],
    
    p_kp : [0, 10],
    theta_kp : [0, 60],

    p_kaon : [0, 5.5],
    theta_kaon : [0, 45],
    phi_kaon : [-180, 180],


    p_km_ftof : [0, 10],
    theta_km_ftof : [0, 60],
    p_km : [0, 10],
    theta_km : [0, 60],

    p_phi : [ 0, 4],
    e_phi : [ 0, 9],

    w         : [0.5, 5],
    w_val     : [0.6, 4.7],
    phi       : [-30, 330],
    angle_ep  : [160, 180],
    angle_ep_diff : [160, 200],
    q2        : [1.2, 4.5],
    vz        : [-20, 15],
    dvz        : [-20, 20],
    missing_mass_small: [-0.3, 0.3],
    prot_chi2 : [-10, 10],
    prot_chi2_small : [-15, 15],
    beta : [0, 1.2],

    q2 : [0, 12],
    xb : [0, 1.15],
    phi_trento : [-180, 180],
    minus_t : [0, 6],
    im_phi : [0, 2.5],
    im_pkm : [0, 5],
    cone_angle: [0,8],
    missing_energy : [-1, 1],
    combs: [1,10],
    coplan_ang : [0, 10],
    pt : [0, 5],
    
    mm2: [-1.5, 1.5],
    mm2_big : [0, 5 ],    
    mass:[0.8, 1.5],
    prokm: [1, 5.5],
    ekp: [1, 5],
    kpkm_big: [0.75, 2.2],
    mm_ep: [0, 5],
    epkp_big : [0.0, 0.7],
    mm_epkm_big:[0.0, 0.7],
    pt_big: [0, 0.25],
    missing_energy_big: [-0.2, 0.5],
    prkm_big : [0.5, 5],

    epkpkm_big : [-0.015, 0.015],

    res_p : [-1,1],
    res_theta : [-1,1],
    res_phi : [-1,1],
    res_phys : [-1,1],

    cos_cm : [-1, 1]   
]

lim = tighter_kin_bounds

def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}




histos = new ConcurrentHashMap()
histoBuilders = [
    //p_ele      : { title -> limited_h1(title, 200, lim.p_ele) },
    //theta_ele  : { title -> limited_h1(title, 200, lim.theta_ele) },
    phi_ele    : { title -> limited_h1(title, 200, lim.phi_ele) },

    p_pro      : { title -> limited_h1(title, 200, lim.p_pro) },
    theta_pro  : { title -> limited_h1(title, 200, lim.theta_pro) },
    phi_pro    : { title -> limited_h1(title, 200, lim.phi_pro) },

    p_kaon      : { title -> limited_h1(title, 200, lim.p_kaon) },
    theta_kaon  : { title -> limited_h1(title, 200, lim.theta_kaon) },
    phi_kaon    : { title -> limited_h1(title, 200, lim.phi_kaon) },

    w        : { title -> limited_h1(title, 200, lim.w) },
    w_val    : { title -> limited_h1(title, 200, lim.w_val) },
    p_ele    : { title -> limited_h1(title, 75, lim.p_ele) },
    p_pro    : { title -> limited_h1(title, 200, lim.p_pro) },
    vz       : { title -> limited_h1(title, 200, lim.vz) },
    angle_ep : { title -> limited_h1(title, 200, lim.angle_ep) },
    theta_p  : { title -> limited_h1(title, 200, lim.theta_pro_ftof) },
    theta_ele: { title -> limited_h1(title, 75, lim.theta_ele) },
    pro_chi2 : { title -> limited_h1(title, 200, lim.prot_chi2_small) },
    missing_mass_small: { title -> limited_h1(title, 200, lim.missing_mass_small) }, // added by brandon c . to validate
    combs    : { title -> limited_h1(title, 200, lim.combs) },
    coplan_angle :{ title -> limited_h1(title, 200, lim.coplan_ang) },

    im_ekpkmX : { title -> limited_h1(title, 200, lim.mm2 ) },
    im_epkX : { title -> limited_h1(title, 200, lim.mm2 ) },
    im_ekpX : { title -> limited_h1(title, 200, lim.mm2_big ) },
    missing_e : { title -> limited_h1(title, 200, lim.mm2 ) },
    missing_p :{ title -> limited_h1(title, 200, lim.mm2 ) },
    im_pkm : { title -> limited_h1(title, 200, lim.mm2_big ) },
    im_kpkm : { title -> limited_h1(title, 150, lim.mass ) },
    phitrento_small: { title -> limited_h1(title, 150, lim.phi_trento ) },
    phitrento: { title -> limited_h1(title, 10, lim.phi_trento ) },

    res_p : { title -> limited_h1(title, 500, lim.res_p ) },
    res_theta : { title -> limited_h1(title, 500, lim.res_theta ) },
    res_phi : { title -> limited_h1(title, 500, lim.res_phi ) },
    
    res_physics : { title -> limited_h1(title, 500, lim.res_phys ) },
    q2 : { title -> limited_h1(title, 500, lim.q2 ) },
    t  : { title -> limited_h1(title, 500, lim.minus_t ) },
    w  : { title -> limited_h1(title, 500, lim.w ) },
    xb : { title -> limited_h1(title, 500, lim.xb ) },
    pt : { title -> limited_h1(title, 500, lim.pt ) },
 
    cm_cos : { title -> limited_h1(title, 10, lim.cos_cm ) },
    
]

histoBuilders2 = [
    p_ele_theta    : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.theta_ele) },
    p_pro_theta    : { title -> limited_h2(title, 200, 200, lim.p_pro, lim.theta_pro_ftof) },
    p_kp_theta     : { title -> limited_h2(title, 200, 200, lim.p_kp, lim.theta_kp) },
    p_km_theta     : { title -> limited_h2(title, 200, 200, lim.p_km, lim.theta_km) },
    p_vz           : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.vz)},
    p_dvz          : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.dvz)},
		      
    phi_theta_ele : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_ele) },
    phi_theta_proton : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_pro_ftof) },
    phi_theta_kp : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_kp) },
    phi_theta_km : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_km) },
        
    phi_vz           : { title -> limited_h2(title, 200, 200, lim.phi, lim.vz) },
    phi_dvz          : { title -> limited_h2(title, 200, 200, lim.phi, lim.dvz) },
    phi_dp           : { title -> limited_h2(title, 100, 100, lim.phi, lim.dp_ele) },
    phi_theta        : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_ele) },
    phi_w            : { title -> limited_h2(title, 200, 200, lim.phi, lim.w) },
    phi_ele_p        : { title -> limited_h2(title, 200, 200, lim.phi, lim.p_ele) },
    phi_pro_p        : { title -> limited_h2(title, 200, 200, lim.phi, lim.p_pro) },
    phi_kp_p        : { title -> limited_h2(title, 200, 200, lim.phi, lim.p_kp) },
    phi_km_p        : { title -> limited_h2(title, 200, 200, lim.phi, lim.p_km) },

    beta_p : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.beta) },
    
    theta_dvz       : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dvz) },
    theta_vz       : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.vz) },
    theta_el_theta_pro : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.theta_ele) },
    theta_ele_vz     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.vz) },
    theta_pro_vz     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.vz) },
    theta_w_ele      : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.w) },
    theta_ele_dp     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dp_ele) },
    theta_pro_dp     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dp_pro) },
    
    w_q2             : { title -> limited_h2(title, 200, 200, lim.w, lim.q2) },
    xb_q2            : { title -> limited_h2(title, 200, 200, lim.xb, lim.q2) },
    t_q2             : { title -> limited_h2(title, 200, 200, lim.minus_t, lim.q2) },
    phitrento_q2     : { title -> limited_h2(title, 200, 200, lim.phi_trento, lim.q2) },
    phitrento_t      : { title -> limited_h2(title, 200, 200, lim.phi_trento, lim.minus_t) },
    phitrento_xb     : { title -> limited_h2(title, 200, 200, lim.phi_trento, lim.xb) },
    phitrento_imkpkm : { title -> limited_h2(title, 200, 200, lim.phi_trento, lim.mass) },

    t_xb             : { title -> limited_h2(title, 200, 200, lim.minus_t, lim.xb) },
    t_w              : { title -> limited_h2(title, 200, 200, lim.minus_t, lim.w) },
    mntm_t           : { title -> limited_h2(title, 200, 200, lim.p_kp, lim.minus_t) },
    mntm_q2          : { title -> limited_h2(title, 200, 200, lim.p_kp, lim.q2) },
    
    p_cone_angle     : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.cone_angle) },    
    cone_angle2D     : { title -> limited_h2(title, 200, 200, lim.cone_angle, lim.cone_angle) },    
    p_w_ele          : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.w) },
    p_phi_theta_sum  : { title -> limited_h2(title, 200, 200, lim.p_phi, lim.theta_sum)},
    e_phi_theta_sum  : { title -> limited_h2(title, 250, 250, lim.e_phi, lim.theta_sum)},
    dalitz           : { title -> limited_h2(title, 200, 200, lim.im_phi, lim.pkm) },
    beta_missing_energy   : { title ->limited_h2(title, 200, 200, lim.beta, lim.missing_energy) },
    
    p_ele_theta_ftof      : { title -> limited_h2(title, 200, 200, lim.p_ele_ftof, lim.theta_ele_ftof) },
    p_pro_theta_ftof      : { title -> limited_h2(title, 200, 200, lim.p_pro_ftof, lim.theta_pro_ftof) },
    

    w_theta_sum : { title -> limited_h2(title, 200, 200, lim.w, lim.theta_sum) },
    p_theta_sum : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.theta_sum) },
    vz_theta_sum : { title -> limited_h2(title, 200, 200, lim.vz, lim.theta_sum) },

    w_delta_theta_proton : {title -> limited_h2(title, 200, 200, lim.w, lim.dtheta_pro) },

    w_missing_mass_small: { title -> limited_h2(title, 200, 200, lim.w, lim.missing_mass_small) }, // added by brandon c. to validate
    vz_missing_mass_small: { title -> limited_h2(title, 200,200, lim.vz, lim.missing_mass_small) }, // added by brandon c. to validate

    theta_sum_missing_mass_small : { title -> limited_h2(title, 200, 200, lim.missing_mass_small, lim.theta_sum) },

    p_chi2 : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.prot_chi2)},
    theta_chi2 : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.prot_chi2)},
    beta_chi2 : { title -> limited_h2(title, 200, 200, lim.beta, lim.prot_chi2)},

    pro_chi2_pro_theta : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.prot_chi2)},
    pro_chi2_pro_theta_ftof : { title -> limited_h2(title, 200, 200, lim.theta_pro_ftof, lim.prot_chi2)},
    pro_chi2_theta_sum : { title -> limited_h2(title, 200, 200, lim.prot_chi2, lim.theta_sum)},

    // Y_axis_X_axis naming convention - stupid but too late to fix
    im_prokm_mm_ekp : { title -> limited_h2(title, 200, 200,  lim.ekp, lim.prokm)},
    im_kpkm_mm_ep : { title -> limited_h2(title, 200, 200, lim.mm_ep, lim.kpkm_big )},
    mm_epkp_mm_epkm : {title -> limited_h2(title, 200, 200, lim.mm_epkm_big, lim.epkp_big)},
    pt_missing_energy : {title -> limited_h2(title, 200, 200, lim.missing_energy_big, lim.pt_big)},
    im_kpkm_mm_epkpkmX: {title -> limited_h2(title, 200, 200, lim.epkpkm_big, lim.kpkm_big )},
    im_kpkm_im_prkm : { title -> limited_h2(title, 200, 200, lim.prkm_big, lim.kpkm_big)},

    cm_cos_minut_t : {title->limited_h2(title, 200, 200, lim.minus_t, lim.cos_cm)}
    
            
]

def getPhiShifted(lv, sec) {
    def fi0 = Math.toDegrees(Math.atan2(lv.py(),lv.px()))
    if( fi0 < 0 && sec > 1 ) fi0 +=360
    return fi0
}

def shiftPhi(temp_phi) {
    def phi = Math.toDegrees(temp_phi)
    return (phi > 150) ? phi - 180 : phi + 180
}



def getPKin(beam, target, electron, proton) {

    def eX = beam + target - electron 
    def w = eX.mass()
    def epX = beam + target - electron - proton
    def missing_mass2 = epX.mass2()

    def q = beam - electron 
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))


    def zaxis = new Vector3(0, 0, 1)
    def diff_phi = Math.abs(electron.phi() - proton.phi())

    def enorm = electron.vect().cross(zaxis)
    def pnorm = proton.vect().cross(zaxis)
    def phi = enorm.theta(pnorm)
    //println(' el ' + Math.toDegrees(electron.phi() ))
    //println(' pr ' + Math.toDegrees(proton.phi()))

    def theta_gamma = Math.toDegrees(epX.theta())
    def theta_egamma = epX.vect().theta(electron.vect())

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle: phi, angle_diff: diff_phi,
            missing_mass: missing_mass2, missing_energy: epX.e(), 
	    theta_gamma:theta_gamma, theta_egamma:theta_egamma, epx : epX, ex : eX ]
}


def getPhysicsKin(beam, target, el, pro, kp, km){

    def eX = beam + target - el
    def epX = beam + target - el - pro
    def ekpX = beam + target - el - kp    
    def ekmX = beam + target - el - km
    def epkpkmX = beam + target - el - pro - kp -km
    def epkpX = beam + target - el - pro - kp
    def epkmX = beam + target - el - pro - km
    def ekpkmX = beam + target - el - kp - km

    def cplpro = Math.toDegrees(ekpkmX.vect().theta(pro.vect()))
    def cplkm = Math.toDegrees(epkpX.vect().theta(km.vect()))
    def cplkp = Math.toDegrees(epkmX.vect().theta(kp.vect()))

    def pt = Math.sqrt(epkpkmX.px()**2 + epkpkmX.py()**2)

    def lmda = pro+km
    def vphi = kp+km 
    def q = beam-el  
    def q2 = -q.mass2()  
    def xb = q2/(2*target.mass()*q.e())
    def t = 2*target.mass()*(pro.e()-target.mass()) 
    def nu = q2/(2*target.mass()*xb)
    def y = nu/beam.e() 

    def v3l = beam.vect().cross(el.vect())
    def v3h = pro.vect().cross((beam-el).vect()) 
    def trento = Math.toDegrees(Math.acos(v3l.dot(v3h)/(v3l.mag()*v3h.mag())))
    if(v3l.dot(pro.vect())<0) trento=-trento

    ///#v3l = beam.Vect().Cross(el.Vect())
    def v3hphi = vphi.vect().cross((beam-el).vect()) 
    def trentophi = Math.toDegrees(Math.acos(v3l.dot(v3h)/(v3l.mag()*v3h.mag())))
    if(v3l.dot(vphi.vect())<0) trentophi=-trentophi
    
    
    def temp_boost_phicm = vphi.boostVector()
    def zero_vect = new Vector3(0,0,0)
    def boost_phicm = zero_vect-temp_boost_phicm
    def boost_kp = LorentzVector.withPID(321, kp.px(), kp.py(), kp.pz())
    def boost_km = LorentzVector.withPID(-321, km.px(), km.py(), km.pz())
    boost_kp.boost(boost_phicm)
    boost_km.boost(boost_phicm)
    def cos_theta_cm_kp = Math.cos(boost_kp.vect().theta())
        
    return [eX : eX, epX: epX, ekpX : ekpX, ekmX: ekmX,  epkpkmX: epkpkmX, epkpX: epkpX, epkmX: epkmX,
	    ekpkmX : ekpkmX, cplpro: cplpro, cplkm: cplkm, cplkp: cplkp, pt: pt, lmda: lmda,
	    vphi: vphi, w: eX.mass(), q2: q2, xb: xb, t: t, nu: nu, y: y,
	    trento: trento, trentophi: trentophi, cos_theta_cm_kp: cos_theta_cm_kp]

}

def fillProtonID(event, histos, electron, pro, ipro, comment){

    histos.computeIfAbsent('pro_chi2_'+comment, histoBuilders.pro_chi2).fill(event.chi2pid[ipro]) 
    histos.computeIfAbsent('p_pro_chi2_'+comment, histoBuilders2.p_chi2).fill(pro.p(), event.chi2pid[ipro]) 
    histos.computeIfAbsent('theta_pro_chi2_'+comment, histoBuilders2.theta_chi2).fill(Math.toDegrees(pro.theta()), event.chi2pid[ipro]) 
    histos.computeIfAbsent('beta_pro_chi2_'+comment, histoBuilders2.beta_chi2).fill(event.beta[ipro], event.chi2pid[ipro]) 
    histos.computeIfAbsent('p_pro_beta_'+comment, histoBuilders2.beta_p).fill( pro.p(), event.beta[ipro])

}

def fillKaonPlusID(event, histos, electron, kp, ikp, comment){

    histos.computeIfAbsent('kp_chi2_'+comment, histoBuilders.pro_chi2).fill(event.chi2pid[ikp]) 
    histos.computeIfAbsent('p_kp_chi2_'+comment, histoBuilders2.p_chi2).fill(kp.p(), event.chi2pid[ikp]) 
    histos.computeIfAbsent('theta_kp_chi2_'+comment, histoBuilders2.theta_chi2).fill(Math.toDegrees(kp.theta()), event.chi2pid[ikp]) 
    histos.computeIfAbsent('beta_kp_chi2_'+comment, histoBuilders2.beta_chi2).fill(event.beta[ikp], event.chi2pid[ikp]) 
    histos.computeIfAbsent('p_kp_beta_'+comment, histoBuilders2.beta_p).fill( kp.p(), event.beta[ikp])
}

def fillKaonMinusID(event, histos, electron, km, ikm, comment){

    histos.computeIfAbsent('km_chi2_'+comment, histoBuilders.pro_chi2).fill(event.chi2pid[ikm]) 
    histos.computeIfAbsent('p_km_chi2_'+comment, histoBuilders2.p_chi2).fill(km.p(), event.chi2pid[ikm]) 
    histos.computeIfAbsent('theta_km_chi2_'+comment, histoBuilders2.theta_chi2).fill(Math.toDegrees(km.theta()), event.chi2pid[ikm]) 
    histos.computeIfAbsent('beta_km_chi2_'+comment, histoBuilders2.beta_chi2).fill(event.beta[ikm], event.chi2pid[ikm]) 
    histos.computeIfAbsent('p_km_beta_'+comment, histoBuilders2.beta_p).fill( km.p(), event.beta[ikm])
}


def fillElectronAccp(event, histos, ele, comment){

    histos.computeIfAbsent('p_ele_'+comment, histoBuilders.p_ele).fill(ele.p())
    histos.computeIfAbsent('theta_ele_'+comment, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_ele_'+comment, histoBuilders.phi_ele).fill(Math.toDegrees(ele.phi()))
    
}


def fillProtonAccp(event, histos, pro, comment){

    histos.computeIfAbsent('p_pro_'+comment, histoBuilders.p_pro).fill(pro.p())
    histos.computeIfAbsent('theta_pro_'+comment, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_pro_'+comment, histoBuilders.phi_pro).fill(Math.toDegrees(pro.phi()))

}

def fillKaonPAccp(event, histos, kp, comment){

    histos.computeIfAbsent('p_kp_'+comment, histoBuilders.p_kaon).fill(kp.p())
    histos.computeIfAbsent('theta_kp_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_kp_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(kp.phi()))


}


def fillKaonMAccp(event, histos, km, comment){

    histos.computeIfAbsent('p_km_'+comment, histoBuilders.p_kaon).fill(km.p())
    histos.computeIfAbsent('theta_km_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(km.theta()))
    histos.computeIfAbsent('phi_km_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(km.phi()))


}


//feed in event, hist, 
// electron proton kaonP kaonM are maps
// ele pro kp km are lorentzvectors
def fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, comment ){

    histos.computeIfAbsent('p_ele_'+comment, histoBuilders.p_ele).fill(ele.p())
    histos.computeIfAbsent('theta_ele_'+comment, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_ele_'+comment, histoBuilders.phi_ele).fill(Math.toDegrees(ele.phi()))

    histos.computeIfAbsent('p_pro_'+comment, histoBuilders.p_pro).fill(pro.p())
    histos.computeIfAbsent('theta_pro_'+comment, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_pro_'+comment, histoBuilders.phi_pro).fill(Math.toDegrees(pro.phi()))

    histos.computeIfAbsent('p_kp_'+comment, histoBuilders.p_kaon).fill(kp.p())
    histos.computeIfAbsent('theta_kp_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_kp_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(kp.phi()))

    histos.computeIfAbsent('p_km_'+comment, histoBuilders.p_kaon).fill(km.p())
    histos.computeIfAbsent('theta_km_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(km.theta()))
    histos.computeIfAbsent('phi_km_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(km.phi()))


    histos.computeIfAbsent('p_ele_theta_'+comment, histoBuilders2.p_ele_theta).fill(ele.p(), Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('p_pro_theta_'+comment, histoBuilders2.p_pro_theta).fill(pro.p(), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('p_kp_theta_'+comment, histoBuilders2.p_kp_theta).fill(kp.p(), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('p_km_theta_'+comment, histoBuilders2.p_km_theta).fill(km.p(), Math.toDegrees(km.theta()))
    
    histos.computeIfAbsent('phi_ele_theta_'+comment, histoBuilders2.phi_theta_ele).fill(shiftPhi(ele.phi()) , Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_pro_theta_'+comment, histoBuilders2.phi_theta_proton).fill(shiftPhi(pro.phi()), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_kp_theta_'+comment, histoBuilders2.phi_theta_kp).fill(shiftPhi(kp.phi()), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_km_theta_'+comment, histoBuilders2.phi_theta_km).fill(shiftPhi(km.phi()), Math.toDegrees(km.theta()))

    histos.computeIfAbsent('p_ele_vz_'+comment, histoBuilders2.p_vz).fill(ele.p(), event.vz[electron.index])
    histos.computeIfAbsent('p_pro_vz_'+comment, histoBuilders2.p_vz).fill(pro.p(), event.vz[proton.index])
    histos.computeIfAbsent('p_kp_vz_'+comment, histoBuilders2.p_vz).fill(kp.p(), event.vz[kaonP.index])
    histos.computeIfAbsent('p_km_vz_'+comment, histoBuilders2.p_vz).fill(km.p(), event.vz[kaonM.index])

    histos.computeIfAbsent('phi_ele_vz_'+comment, histoBuilders2.phi_vz).fill(shiftPhi(ele.phi()) , event.vz[electron.index])
    histos.computeIfAbsent('phi_pro_vz_'+comment, histoBuilders2.phi_vz).fill(shiftPhi(pro.phi()) , event.vz[proton.index])
    histos.computeIfAbsent('phi_kp_vz_'+comment, histoBuilders2.phi_vz).fill(shiftPhi(kp.phi()) , event.vz[kaonP.index])
    histos.computeIfAbsent('phi_km_vz_'+comment, histoBuilders2.phi_vz).fill(shiftPhi(km.phi()) , event.vz[kaonM.index])

    histos.computeIfAbsent('theta_ele_vz_'+comment, histoBuilders2.theta_vz).fill(Math.toDegrees(ele.theta()), event.vz[electron.index])
    histos.computeIfAbsent('theta_pro_vz_'+comment, histoBuilders2.theta_vz).fill(Math.toDegrees(pro.theta()), event.vz[proton.index])
    histos.computeIfAbsent('theta_kp_vz_'+comment, histoBuilders2.theta_vz).fill(Math.toDegrees(kp.theta()), event.vz[kaonP.index])
    histos.computeIfAbsent('theta_km_vz_'+comment, histoBuilders2.theta_vz).fill(Math.toDegrees(km.theta()), event.vz[kaonM.index])
    
    //vertex differences
    histos.computeIfAbsent('p_ele_dvz_elepro_'+comment, histoBuilders2.p_dvz).fill(ele.p(), (event.vz[electron.index] - event.vz[proton.index]))
    histos.computeIfAbsent('p_pro_dvz_elepro_'+comment, histoBuilders2.p_dvz).fill(pro.p(), (event.vz[electron.index] - event.vz[proton.index]))
    histos.computeIfAbsent('p_kp_dvz_kpkm_'+comment, histoBuilders2.p_dvz).fill(kp.p(), (event.vz[kaonP.index] - event.vz[kaonM.index]))
    histos.computeIfAbsent('p_km_dvz_kpkm_'+comment, histoBuilders2.p_dvz).fill(km.p(), (event.vz[kaonP.index] - event.vz[kaonM.index]))

    histos.computeIfAbsent('theta_ele_dvz_elepro_'+comment, histoBuilders2.theta_dvz).fill(Math.toDegrees(ele.theta()), (event.vz[electron.index] - event.vz[proton.index]))
    histos.computeIfAbsent('theta_pro_dvz_elepro_'+comment, histoBuilders2.theta_dvz).fill(Math.toDegrees(pro.theta()), (event.vz[electron.index] - event.vz[proton.index]))
    histos.computeIfAbsent('theta_kp_dvz_kpkm_'+comment, histoBuilders2.theta_dvz).fill(Math.toDegrees(kp.theta()), (event.vz[kaonP.index] - event.vz[kaonM.index]))
    histos.computeIfAbsent('theta_km_dvz_kpkm_'+comment, histoBuilders2.theta_dvz).fill(Math.toDegrees(km.theta()), (event.vz[kaonP.index] - event.vz[kaonM.index]))

    histos.computeIfAbsent('phi_ele_dvz_elepro_'+comment, histoBuilders2.phi_dvz).fill(shiftPhi(ele.phi()) , (event.vz[electron.index] - event.vz[proton.index]))
    histos.computeIfAbsent('phi_pro_dvz_elepro_'+comment, histoBuilders2.phi_dvz).fill(shiftPhi(pro.phi()) , (event.vz[electron.index] - event.vz[proton.index]))
    histos.computeIfAbsent('phi_kp_dvz_kpkm_'+comment, histoBuilders2.phi_dvz).fill(shiftPhi(kp.phi()) , (event.vz[kaonP.index] - event.vz[kaonM.index]))
    histos.computeIfAbsent('phi_km_dvz_kpkm_'+comment, histoBuilders2.phi_dvz).fill(shiftPhi(km.phi()) , (event.vz[kaonP.index] - event.vz[kaonM.index]))
    
    //beta
    histos.computeIfAbsent('p_pro_beta_'+comment, histoBuilders2.beta_p).fill( pro.p(), event.beta[proton.index])
    histos.computeIfAbsent('p_kp_beta_'+comment, histoBuilders2.beta_p).fill( kp.p(), event.beta[kaonP.index])
    histos.computeIfAbsent('p_km_beta_'+comment, histoBuilders2.beta_p).fill( km.p(), event.beta[kaonM.index])
    

    // kinematics
    histos.computeIfAbsent('p_kp_theta_sum_'+comment, histoBuilders2.p_theta_sum).fill( kp.p(), Math.toDegrees(kp.theta()) + Math.toDegrees(km.theta()))
    histos.computeIfAbsent('p_km_theta_sum_'+comment, histoBuilders2.p_theta_sum).fill( km.p(), Math.toDegrees(kp.theta()) + Math.toDegrees(km.theta()))
    histos.computeIfAbsent('w_theta_sum_'+comment, histoBuilders2.w_theta_sum).fill( phykin.w, Math.toDegrees(kp.theta()) + Math.toDegrees(km.theta()))
    histos.computeIfAbsent('mass_phi_theta_sum_'+comment, histoBuilders2.p_phi_theta_sum).fill( phykin.vphi.mass(), Math.toDegrees(kp.theta()) + Math.toDegrees(km.theta()))
    histos.computeIfAbsent('e_phi_theta_sum_'+comment, histoBuilders2.e_phi_theta_sum).fill( phykin.vphi.e(), Math.toDegrees(kp.theta()) + Math.toDegrees(km.theta()))
    histos.computeIfAbsent('e_pkm_theta_sum_'+comment, histoBuilders2.e_phi_theta_sum).fill( (pro+km).e(), Math.toDegrees(kp.theta()) + Math.toDegrees(km.theta()))
    
    histos.computeIfAbsent('cpl_pro_'+comment,histoBuilders.coplan_angle).fill( phykin.cplpro )
    histos.computeIfAbsent('cpl_pro_'+comment,histoBuilders.coplan_angle).fill( phykin.cplpro )
    histos.computeIfAbsent('cpl_km_'+comment,histoBuilders.coplan_angle).fill( phykin.cplpro )
    
    //physics observables
    histos.computeIfAbsent('w_q2_'+comment,histoBuilders2.w_q2).fill( phykin.w, phykin.q2 )
    histos.computeIfAbsent('xb_q2_'+comment,histoBuilders2.xb_q2).fill( phykin.xb, phykin.q2 )
    histos.computeIfAbsent('t_q2_'+comment,histoBuilders2.t_q2).fill( phykin.t, phykin.q2 )
    histos.computeIfAbsent('phitrento_q2_'+comment,histoBuilders2.phitrento_q2).fill( phykin.trentophi, phykin.q2 )
    histos.computeIfAbsent('phitrento_t_'+comment,histoBuilders2.phitrento_t).fill( phykin.trentophi, phykin.t)
    histos.computeIfAbsent('phitrento_xb_'+comment,histoBuilders2.phitrento_xb).fill( phykin.trentophi, phykin.xb)

    histos.computeIfAbsent('t_xb_'+comment,histoBuilders2.t_xb).fill( phykin.t, phykin.xb )
    histos.computeIfAbsent('t_w_'+comment,histoBuilders2.t_w).fill( phykin.t, phykin.w )

    histos.computeIfAbsent('p_kp_t_'+comment,histoBuilders2.mntm_t).fill( kp.p(),  phykin.t)
    histos.computeIfAbsent('p_pro_t_'+comment,histoBuilders2.mntm_t).fill( pro.p(),  phykin.t)
    histos.computeIfAbsent('p_ekpkmX_t_'+comment,histoBuilders2.mntm_t).fill( phykin.ekpkmX.p(),  phykin.t)
    histos.computeIfAbsent('p_kp_q2_'+comment,histoBuilders2.mntm_q2).fill( kp.p(),  phykin.q2)

    histos.computeIfAbsent('cm_cos_theta_kp_'+comment,histoBuilders.cm_cos).fill(phykin.cos_theta_cm_kp)
    histos.computeIfAbsent('cm_cos_theta_kp_vs_t_'+comment,histoBuilders2.cm_cos_minut_t).fill(phykin.t, phykin.cos_theta_cm_kp)
    
    def helicity_state = event.helicity < 0 ? 'neg' : 'pos'
    histos.computeIfAbsent('phitrento_small_'+helicity_state+'_'+comment,histoBuilders.phitrento_small).fill(phykin.trentophi)
    histos.computeIfAbsent('phitrento_'+helicity_state+'_'+comment,histoBuilders.phitrento).fill(phykin.trentophi)    
    
    //mass info
    // 2D
    histos.computeIfAbsent('im_prokm_mm_ekp_'+comment,histoBuilders2.im_prokm_mm_ekp).fill(phykin.ekpX.mass(), (pro+km).mass() )
    histos.computeIfAbsent('im_prokp_mm_ekm_'+comment,histoBuilders2.im_prokm_mm_ekp).fill( phykin.ekmX.mass(), (pro+kp).mass() )
    histos.computeIfAbsent('im_kpkm_mm_ep_'+comment,histoBuilders2.im_kpkm_mm_ep).fill( phykin.epX.mass(), (kp+km).mass())
    histos.computeIfAbsent('mm_epkp_mm_epkm_'+comment,histoBuilders2.mm_epkp_mm_epkm).fill( phykin.epkmX.mass(), phykin.epkpX.mass())
    histos.computeIfAbsent('pT_missing_energy_'+comment,histoBuilders2.pt_missing_energy).fill( phykin.epkpkmX.e(), phykin.pt )
    histos.computeIfAbsent('im_kpkm_mm_epkpkmX_'+comment,histoBuilders2.im_kpkm_mm_epkpkmX).fill( phykin.epkpkmX.mass() , (kp+km).mass())
    histos.computeIfAbsent('im_kpkm_im_pkm_'+comment,histoBuilders2.im_kpkm_im_prkm).fill( (pro+km).mass(), (kp+km).mass() )
    
    //mass info
    // 1d
    histos.computeIfAbsent('mm_epkpkmX_'+comment,histoBuilders.missing_e).fill(phykin.epkpkmX.e())
    histos.computeIfAbsent('mm_ekpkmX_'+comment,histoBuilders.im_ekpkmX).fill(phykin.ekpkmX.mass2())
    histos.computeIfAbsent('mm_epkpX_'+comment,histoBuilders.im_epkX).fill(phykin.epkpX.mass2())
    histos.computeIfAbsent('mm_epkmX_'+comment,histoBuilders.im_epkX).fill(phykin.epkmX.mass2())
    histos.computeIfAbsent('mm_ekpX_'+comment,histoBuilders.im_ekpX).fill(phykin.ekpX.mass2())
    histos.computeIfAbsent('missing_energy_'+comment,histoBuilders.missing_e).fill(phykin.epkpkmX.e())
    histos.computeIfAbsent('missing_p_'+comment,histoBuilders.missing_p).fill(phykin.pt)
    histos.computeIfAbsent('im_pkm_'+comment,histoBuilders.im_pkm).fill( (pro+km).mass() )
    histos.computeIfAbsent('im_pkp_'+comment,histoBuilders.im_pkm).fill( (pro+kp).mass() )
    histos.computeIfAbsent('im_kpkm_'+comment,histoBuilders.im_kpkm).fill( (kp+km).mass() )
    histos.computeIfAbsent('phitrento_imkpkm_'+comment,histoBuilders2.phitrento_imkpkm).fill( phykin.trentophi, (kp+km).mass() )

    //physics
    // 1d
    histos.computeIfAbsent('q2_'+comment, histoBuilders.q2).fill(phykin.q2)
    histos.computeIfAbsent('t_'+comment, histoBuilders.t).fill(phykin.t)
    histos.computeIfAbsent('w_'+comment, histoBuilders.w).fill(phykin.w)
    histos.computeIfAbsent('xb_'+comment, histoBuilders.xb).fill(phykin.xb)
    histos.computeIfAbsent('pt_'+comment, histoBuilders.pt).fill(phykin.pt)
    histos.computeIfAbsent('phitrento_'+comment, histoBuilders.phitrento_small).fill(phykin.trentophi)
   
        
}

//fill when event for phi is present in simulation 
def fillSimulationResolutions(histos, event, ele, pro, kp, km, gen_ele, gen_pro, gen_kp, gen_km, comment ){
    
    histos.computeIfAbsent('p_ele_res_'+comment, histoBuilders.res_p).fill(ele.p() - gen_ele.p())
    histos.computeIfAbsent('theta_ele_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(ele.theta() - gen_ele.theta()))
    histos.computeIfAbsent('phi_ele_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(ele.phi() - gen_ele.phi()))

    histos.computeIfAbsent('p_pro_res_'+comment, histoBuilders.res_p).fill(pro.p() - gen_pro.p())
    histos.computeIfAbsent('theta_pro_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(pro.theta() - gen_pro.theta()))
    histos.computeIfAbsent('phi_pro_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(pro.phi() - gen_pro.phi()))
    
    histos.computeIfAbsent('p_kp_res_'+comment, histoBuilders.res_p).fill(kp.p() - gen_kp.p())
    histos.computeIfAbsent('theta_kp_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(kp.theta() - gen_kp.theta()))
    histos.computeIfAbsent('phi_kp_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(kp.phi() - gen_kp.phi()))
    
    histos.computeIfAbsent('p_kp_res_'+comment, histoBuilders.res_p).fill(km.p() - gen_km.p())
    histos.computeIfAbsent('theta_km_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(km.theta() - gen_km.theta()))
    histos.computeIfAbsent('phi_km_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(km.phi() - gen_km.phi()))
    
    
}

def fillSimulationHistograms(histos, beam, target, ele, pro, kp, km, comment){

    def phykin = getPhysicsKin(beam, target, ele, pro, kp, km)

    histos.computeIfAbsent('p_ele_'+comment, histoBuilders.p_ele).fill(ele.p())
    histos.computeIfAbsent('theta_ele_'+comment, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_ele_'+comment, histoBuilders.phi_ele).fill(Math.toDegrees(ele.phi()))

    histos.computeIfAbsent('p_pro_'+comment, histoBuilders.p_pro).fill(pro.p())
    histos.computeIfAbsent('theta_pro_'+comment, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_pro_'+comment, histoBuilders.phi_pro).fill(Math.toDegrees(pro.phi()))

    histos.computeIfAbsent('p_kp_'+comment, histoBuilders.p_kaon).fill(kp.p())
    histos.computeIfAbsent('theta_kp_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_kp_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(kp.phi()))

    histos.computeIfAbsent('p_km_'+comment, histoBuilders.p_kaon).fill(km.p())
    histos.computeIfAbsent('theta_km_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(km.theta()))
    histos.computeIfAbsent('phi_km_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(km.phi()))

    histos.computeIfAbsent('p_pro_theta_'+comment, histoBuilders2.p_pro_theta).fill(pro.p(), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('p_kp_theta_'+comment, histoBuilders2.p_kp_theta).fill(kp.p(), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('p_km_theta_'+comment, histoBuilders2.p_km_theta).fill(km.p(), Math.toDegrees(km.theta()))
    
    histos.computeIfAbsent('phi_ele_theta_'+comment, histoBuilders2.phi_theta_ele).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_pro_theta_'+comment, histoBuilders2.phi_theta_proton).fill(Math.toDegrees(pro.phi()), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_kp_theta_'+comment, histoBuilders2.phi_theta_kp).fill(Math.toDegrees(kp.phi()), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_km_theta_'+comment, histoBuilders2.phi_theta_km).fill(Math.toDegrees(km.phi()), Math.toDegrees(km.theta()))
    
    histos.computeIfAbsent('phi_ele_p_'+comment, histoBuilders2.phi_ele_p).fill(Math.toDegrees(ele.phi()), ele.p())
    histos.computeIfAbsent('phi_pro_p_'+comment, histoBuilders2.phi_pro_p).fill(Math.toDegrees(pro.phi()), pro.p())
    histos.computeIfAbsent('phi_kp_p_'+comment, histoBuilders2.phi_kp_p).fill(Math.toDegrees(kp.phi()), kp.p())
    histos.computeIfAbsent('phi_km_p_'+comment, histoBuilders2.phi_km_p).fill(Math.toDegrees(km.phi()), km.p())

    histos.computeIfAbsent('q2_'+comment, histoBuilders.q2).fill(phykin.q2)
    histos.computeIfAbsent('t_'+comment, histoBuilders.t).fill(phykin.t)
    histos.computeIfAbsent('w_'+comment, histoBuilders.w).fill(phykin.w)
    histos.computeIfAbsent('xb_'+comment, histoBuilders.xb).fill(phykin.xb)
    histos.computeIfAbsent('pt_'+comment, histoBuilders.pt).fill(phykin.pt)
    histos.computeIfAbsent('phitrento_'+comment, histoBuilders.phitrento_small).fill(phykin.trentophi)
    
    histos.computeIfAbsent('w_q2_'+comment,histoBuilders2.w_q2).fill( phykin.w, phykin.q2 )
    histos.computeIfAbsent('xb_q2_'+comment,histoBuilders2.xb_q2).fill( phykin.xb, phykin.q2 )
    histos.computeIfAbsent('t_q2_'+comment,histoBuilders2.t_q2).fill( phykin.t, phykin.q2 )
    histos.computeIfAbsent('phitrento_q2_'+comment,histoBuilders2.phitrento_q2).fill( phykin.trentophi, phykin.q2 )
    histos.computeIfAbsent('phitrento_t_'+comment,histoBuilders2.phitrento_t).fill( phykin.trentophi, phykin.t)
    histos.computeIfAbsent('phitrento_xb_'+comment,histoBuilders2.phitrento_xb).fill( phykin.trentophi, phykin.xb)


    histos.computeIfAbsent('t_xb_'+comment,histoBuilders2.t_xb).fill( phykin.t, phykin.xb )
    histos.computeIfAbsent('t_w_'+comment,histoBuilders2.t_w).fill( phykin.t, phykin.w )
 
}

def createFileOut( f_out ){
    return new File(f_out)
}

def writeOutEvent( t_out, temp_l ){
    //t_out.fill(temp_l[0], temp_l[1], temp_l[2])   

    t_out.fill(*temp_l)//.collect{partbank.getFloat(it, eleind)}
    //t_out.fill(temp_l)   
}

def write_to_file = true
def f_out = null
def t1=null
if( write_to_file ){
    f_out = new ROOTFile("phi_mass_results_pidtype1-inb.root")
    t1 = f_out.makeNtuple('clas12','data','phimass:helicity:trento:q2:xb:t')
}

// 1 for SIMULATION, 0 for DATA
def DATATYPE = 0
def proc_method = -1
data_prefix=null
if( DATATYPE == 1 ){
    data_prefix='mc'
    proc_method = 3
}
else{
    data_prefix='data'
    proc_method = 2
}

def topology = 'epkpkm'
def field_setting = 'inbending'
// cut lvl meanings: 0 loose, 1 med, 2 tight
el_cut_strictness_lvl=["ecal_cut_lvl":1,
		       "nphe_cut_lvl":1,
		       "vz_cut_lvl":1,
		       "min_u_cut_lvl":1,
		       "min_v_cut_lvl":1,
		       "min_w_cut_lvl":1,
		       "max_u_cut_lvl":1,
		       "max_v_cut_lvl":1,
		       "max_w_cut_lvl":1,
		       "dcr1_cut_lvl":1,
		       "dcr2_cut_lvl":1,
		       "dcr3_cut_lvl":1,
		       "anti_pion_cut_lvl":1
]

def electron_from_event = new ElectronFromEvent()
electron_from_event.setElectronCutStrictness(el_cut_strictness_lvl)
electron_from_event.setElectronCutParameters(field_setting)
def myElectronCutStrategies = [
    electron_from_event.passElectronStatus,
    electron_from_event.passElectronChargeCut,
    electron_from_event.passElectronEBPIDCut,
    //on top of event  builder
    //electron_from_event.passElectronMinMomentum,
    //electron_from_event.passElectronVertexCut,
    //electron_from_event.passElectronPCALFiducialCut,
    //electron_from_event.passElectronEIEOCut,
    //electron_from_event.passElectronSamplingFractionCut
]

def proton_cand = new HadronID(hadron_id:2212)
def kaonP_cand = new HadronID(hadron_id:321)
def kaonM_cand = new HadronID(hadron_id:-321)

def myProtonCuts = [
    proton_cand.passHadronEBPIDCut,    
    //proton_cand.passTrajChi2Cut
    //proton_cand.passDCFiducialCutChi2Region1,
    //proton_cand.passDCFiducialCutChi2Region2,
    //proton_cand.passDCFiducialCutChi2Region3
]

def myKaonPCuts = [
    kaonP_cand.passHadronEBPIDCut,
    //kaonP_cand.passTrajChi2Cut
    //kaonP_cand.passDCFiducialCutChi2Region1,
    //kaonP_cand.passDCFiducialCutChi2Region2,
    //kaonP_cand.passDCFiducialCutChi2Region3
]

def myKaonMCuts = [    
    kaonM_cand.passHadronEBPIDCut,
    //kaonM_cand.passTrajChi2Cut
    //kaonM_cand.passDCFiducialCutChi2Region1,
    //kaonM_cand.passDCFiducialCutChi2Region2,
    //kaonM_cand.passDCFiducialCutChi2Region3
]


GParsPool.withPool 16, {
    args.eachParallel{fname->
	def reader = new HipoDataSource()
	reader.open(fname)

	cc=0
	while(reader.hasEvent()){
	    if (cc % 5000 == 0) {
                println("Processing " + cc)
            }

	    def event = EventConverter.convert(reader.getNextEvent())
	    
	    def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()
	    def e_present= false
	    def rec_electron_present = false
	    def rec_proton_present = false
	    def rec_kaonP_present = false
	    def rec_kaonM_present = false

	    for( ee in e_cand){
		def e_index = ee.key
		def e_results = ee.value
		def cut_index = 0
		def temp_ele = LorentzVector.withPID(11, event.px[e_index], event.py[e_index], event.pz[e_index])
		if( e_index == 0 && !e_results.contains(false) && Math.toDegrees(temp_ele.theta()) < 35 ){
		    //histos.computeIfAbsent('p_ele_pass_cut_'+cut_index, histoBuilders.p_ele).fill(temp_ele.p())
		    fillElectronAccp(event, histos, temp_ele, 'sim_rec_all')
		    rec_electron_present=true		    
		}
	    }


	    // when checking acceptance fill histograms with ALL detected particles independent of detecting the other ones
	    if( rec_electron_present ){
		(0..<event.npart).findAll{event.pid[it]==2212 && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{ipro->
		    def temp_pro = LorentzVector.withPID(2212, event.px[ipro], event.py[ipro], event.pz[ipro])
		    //define strictly to angular coverage of FD
		    if( Math.toDegrees(temp_pro.theta()) < 35 ){
			fillProtonAccp(event, histos, temp_pro, 'sim_rec_all')
			rec_proton_present = true
		    }
		}

		(0..<event.npart).findAll{event.pid[it]==321 && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{ikp->
		    def temp_kp = LorentzVector.withPID(321, event.px[ikp], event.py[ikp], event.pz[ikp])
		    //define strictly to angular coverage of FD
		    if( Math.toDegrees(temp_kp.theta()) < 35 ){
			fillKaonPAccp(event, histos, temp_kp, 'sim_rec_all')
			rec_kaonP_present=true
		    }		    
		}

		(0..<event.npart).findAll{event.pid[it]==-321 && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{ikm->
		    def temp_km = LorentzVector.withPID(-321, event.px[ikm], event.py[ikm], event.pz[ikm])
		    //define strictly to angular coverage of FD
		    if( Math.toDegrees(temp_km.theta()) < 35 ){
			fillKaonMAccp(event, histos, temp_km, 'sim_rec_all')
			rec_kaonM_present = true
		    }
		}
	    }

	    // boolean flag used to get resolution with simulations
	    def reconstructed_event_present = false
	    def rec_ele, rec_pro, rec_kp, rec_km = null

	    //get electron proton pair first
	    def ep = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
		def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    
		return [index:indx, lv:ele, name:'el', status:event.status[indx] ]
		//return ele
	    }.collectMany{ mele -> 
		(0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		    def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 	
		    def mpro = [index:ipro, lv:pro, name:'pro', status:event.status[ipro]]
		    def ctof = (event.pid[ipro] >= 2000 && event.pid[ipro] < 4000) ? "FTOF" : "CTOF"
		    def pr_fid_stat = myProtonCuts.collect{ pr_test -> pr_test(event,ipro) }
		    //fillProtonID(event, histos, mele.lv, pro, ipro, 'eb_'+ctof)
		    return [mele,mpro]
		}
	    }
	    
	    if( ep.size() == 1 ){	    
		//get Kp and Km now
		def combs = ep.collectMany{mele,mpro->            	
		    //(0..<event.npart).findAll{event.pid[it]==321 && !myKaonPCuts.collect{ kp_test -> kp_test(event,it) }.contains(false) }.findResults{ikp->
		    (0..<event.npart).findAll{ !myKaonPCuts.collect{ kp_test -> kp_test(event,it) }.contains(false) }.findResults{ikp->
			def kp = LorentzVector.withPID(321,event.px[ikp], event.py[ikp], event.pz[ikp]) 		
 			def mkp = [index:ikp, lv:kp,name:'kp',status:event.status[ikp]]
			def ctof = (event.pid[ikp] >= 2000 && event.pid[ikp] < 4000) ? "FTOF" : "CTOF"
			//fillKaonPlusID(event, histos, mele.lv, kp, ikp, 'eb_'+ctof)
  			def kp_fid_stat = myKaonPCuts.collect{ kp_test -> kp_test(event,ikp) }
			return [mele,mpro,mkp]
		    }
		}.collectMany{mele,mpro,mkp->  
		    //(0..<event.npart).findAll{event.pid[it]==-321}.findResults{ikm->
		    (0..<event.npart).findAll{!myKaonMCuts.collect{km_test->km_test(event,it)}.contains(false)}.findResults{ikm->
			def km = LorentzVector.withPID(-321,event.px[ikm], event.py[ikm], event.pz[ikm]) 						
			def mkm = [index:ikm, lv:km, name:'km',status:event.status[ikm]]
			def ctof = (event.pid[ikm] >= 2000 && event.pid[ikm] < 4000) ? "FTOF" : "CTOF"
			def km_fid_stat = myKaonMCuts.collect{ km_test -> km_test(event,ikm) }
			//fillKaonMinusID(event, histos, mele.lv, km, ikm, 'eb_'+ctof)
			return [mele,mpro,mkp,mkm]
		    }
		}.findResults{ele,pro,kp,km->	    
		    //returns a list of map elements
		    return [ele,pro,kp,km]
		}
		
		//check how many combinations occur
 		histos.computeIfAbsent('combsize',histoBuilders.combs).fill(combs.size())
		
		// look into physics observables for the final state
		def all_in_fd = false
		if( combs.size() > 0 ){
		    def check_for_duplicates=false
		    def kp_ftof_p = 0
		    for(comb in combs){
			def (el, pro, kp, km ) = comb
			//def ctof = (event.status[it] >= 2000 && event.status[it] < 4000) ? "FTOF" : "CTOF"
			def el_ctof = (el.status < 0) ? "FTOF" : "CTOF"
			def pr_ctof = (pro.status >= 2000 && pro.status < 4000) ? "FTOF" : "CTOF"
			def kp_ctof = (kp.status >= 2000 && kp.status < 4000) ? "FTOF" : "CTOF"
			def km_ctof = (km.status >= 2000 && km.status < 4000) ? "FTOF" : "CTOF"
			
			if( kp_ctof == "FTOF" ){
			    check_for_duplicates=true
			    kp_ftof_p = kp.lv.p()
			}
			if( combs.size() == 2 && kp_ctof == "CTOF" && check_for_duplicates ){
			    //println('combs size ' + combs.size() + ' kp ftof p ' + kp_ftof_p + ' ' + kp.lv.p())
			}
			    

			def phykin =  getPhysicsKin(beam, target, el.lv, pro.lv, kp.lv, km.lv)
			
			if(el_ctof == "FTOF" && pr_ctof == "FTOF" && 
			   kp_ctof == "FTOF" && km_ctof == "FTOF" ){
			    all_in_fd=true
			    
			}		
		    }
		    if( all_in_fd ) { histos.computeIfAbsent('combsize_all_in_fd',histoBuilders.combs).fill(combs.size()) }
		    if( combs.size() == 1 && all_in_fd ){
			def (electron, proton, kaonP, kaonM ) = combs[0]
			def ele = electron.lv
			def pro = proton.lv
			def kp = kaonP.lv
			def km = kaonM.lv
			def phykin =  getPhysicsKin(beam, target, ele, pro, kp, km)
			
			def comment='combssize1'
			//kinematics for all events with only one possible combination of proton, kp, km. i.e. phi candidate
			fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, comment)
			fillProtonID(event, histos, ele, pro, proton.index, comment)
			fillKaonPlusID(event, histos, ele, kp, kaonP.index, comment)
			fillKaonMinusID(event, histos, ele, km, kaonM.index, comment)
			
			//now we want to apply cuts to select the phi events 
			//apply cuts on epkpkmX.e(), and missing proton, missing kaonP and missing kaonM
			def pass_epkpkmXe = phykin.epkpkmX.e() > excl_cuts.epkpkmxe[0] && phykin.epkpkmX.e() < excl_cuts.epkpkmxe[1]
			def pass_epkpX = phykin.epkpX.mass2() > excl_cuts.epkpX[0] && phykin.epkpX.mass2() < excl_cuts.epkpX[1]
			def pass_epkmX = phykin.epkmX.mass2() > excl_cuts.epkmX[0] && phykin.epkmX.mass2() < excl_cuts.epkmX[1]
			def pass_ekpkmX = phykin.ekpkmX.mass2() > excl_cuts.ekpkmX[0] && phykin.ekpkmX.mass2() < excl_cuts.ekpkmX[1]
			def pass_el_p_min = ele.p() > excl_cuts.el_p_min
			def pass_pro_p_min = pro.p() < excl_cuts.pro_p_max
			def pass_kp_p_min = kp.p() < excl_cuts.kp_p_max
			def pass_km_p_min = km.p() < excl_cuts.km_p_max
			def pass_pt_max = phykin.pt < excl_cuts.pt_max
			def pass_q2_min = phykin.q2 > excl_cuts.q2_min
			def pass_w_min = phykin.w > excl_cuts.w_min
			//restrict to defined FD coverage
			def pass_pro_theta_max = Math.toDegrees(pro.theta()) < 35
			def pass_kp_theta_max = Math.toDegrees(kp.theta()) < 35
			def pass_km_theta_max = Math.toDegrees(km.theta()) < 35
			
			//make plots for passing all but one of the missing mass/energy cuts
			// see ekpkmX first - missnig proton
			if( pass_epkpkmXe && pass_epkpX && pass_epkmX ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_but_missing_proton')
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'fail_at_least_missing_proton')
			}
			// see  missnig kp
			if( pass_epkpkmXe && pass_epkpX && pass_ekpkmX  ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_but_missing_kaonP')			    
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'fail_at_least_missing_kaonP')
			}
			// see  missnig km
			if( pass_epkpkmXe && pass_epkmX && pass_ekpkmX  ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_but_missing_kaonM')
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'fail_at_least_missing_kaonM')
			}
			// see missing energy
			if(  pass_epkpX && pass_epkmX && pass_ekpkmX  ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_but_missing_energy')
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'fail_at_least_missing_energy')
			}


			//restrict by cutting on all variables above
			if( pass_epkpkmXe && pass_epkpX && pass_epkmX && pass_ekpkmX ){ 
			    reconstructed_event_present=true
			    
			    rec_ele = ele
			    rec_pro = pro
			    rec_kp  = kp
			    rec_km  = km //maybe fill as map?
			    
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all')
			    fillProtonID(event, histos, ele, pro, proton.index, 'pass_all')
			    fillKaonPlusID(event, histos, ele, kp, kaonP.index, 'pass_all')
			    fillKaonMinusID(event, histos, ele, km, kaonM.index, 'pass_all')
			    
			    def pass_not_resonant1 = (pro+km).mass() < 1.5 || (pro+km).mass() > 1.58
			    def pass_not_resonant2 = (pro+km).mass() < 1.78 || (pro+km).mass() > 1.9
			    // check what happens when we remove resonances vs keep them
			    //removed writing out here and included the additional cuts which help remove background contr. from peak

			    //if( t1 != null && event.helicity != 0 ){// && pass_not_resonant2) {
			    //writeOutEvent(t1, [ (kp+km).mass(), event.helicity, phykin.trentophi, phykin.q2, phykin.xb, phykin.t ] )
			    //}
			    
			    // check the impact of a cut on the maximum momentum of proton and kaons on IM of KPKM
			    if( pass_pro_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_low_pro_p')
			    }
			    if( !pass_pro_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_high_pro_p')
			    }
			    if( pass_kp_p_min && pass_km_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_low_kaonPM_p')
			    }
			    if( !pass_kp_p_min || !pass_km_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_high_kaonPM_p')
			    }
			    
			    // limit theta to only defined FD coverage less than 35 degrees
			    if( pass_pro_theta_max && pass_kp_theta_max && pass_km_theta_max ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_pass_thetaFD')
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_fail_thetaFD')
			    }
			    
	 		    if( (kp+km).mass() < 1.050 && (kp+km).mass() > 1.0100 ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_cut_phi_mass')
				fillProtonID(event, histos, ele, pro, proton.index, 'pass_all_cut_phi_mass')
				fillKaonPlusID(event, histos, ele, kp, kaonP.index, 'pass_all_cut_phi_mass')
				fillKaonMinusID(event, histos, ele, km, kaonM.index, 'pass_all_cut_phi_mass')
			    }

			    			   
			    // now add additional cut on the pid chi2 from EB for hadrons to see if it removes background near the signal			     
 			    for( ii in 1..6){
				def pass_pro_chi2 = Math.abs(event.chi2pid[ proton.index ]) < ii 
				def pass_kp_chi2 = Math.abs(event.chi2pid[ kaonP.index ]) < ii
				def pass_km_chi2 = Math.abs(event.chi2pid[ kaonM.index ]) < ii
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_raw_${ii}sigcut")
				if( pass_pro_chi2 && pass_kp_chi2 && pass_km_chi2 ){
				    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_${ii}sigcut")
				}
				else{
				    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_fail_${ii}sigcut")
				}
			    }

			    //additional cuts that could be used - delta vz of kp-km
			    def pass_dvz_kpkm = (event.vz[kaonP.index] - event.vz[kaonM.index]) < 5 && (event.vz[kaonP.index] - event.vz[kaonM.index]) > -7
			    def pass_dvz_ep = Math.abs(event.vz[electron.index] - event.vz[proton.index]) < 6
			    // pass just dvz kpkm check is phi mass signal improves
			    if( pass_dvz_kpkm ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_dvz_kpkm")
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_fail_dvz_kpkm")
			    }
			    // pass just dvz electron proton
			    if( pass_dvz_ep ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_dvz_ep")
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_fail_dvz_ep")
			    }
			    //passs both the delta vz cuts
			    if( pass_dvz_kpkm && pass_dvz_ep ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_dvz_both")
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_fail_dvz_both")
			    }


			    // now include specific additional cuts to select final events.
			    def pass_pro_chi2 = Math.abs(event.chi2pid[ proton.index ]) < 6 
			    def pass_kp_chi2 = Math.abs(event.chi2pid[ kaonP.index ]) < 6
			    def pass_km_chi2 = Math.abs(event.chi2pid[ kaonM.index ]) < 6			    
			    if( pass_dvz_kpkm &&  pass_pro_theta_max && pass_kp_theta_max && pass_km_theta_max && 
			       pass_kp_p_min && pass_km_p_min && pass_pro_chi2 && pass_kp_chi2 && pass_km_chi2 ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_additional")

				if( t1 != null && event.helicity != 0 ){// && pass_not_resonant2) {
				    writeOutEvent(t1, [ (kp+km).mass(), event.helicity, phykin.trentophi, phykin.q2, phykin.xb, phykin.t ] )
				}
				// now to have a look at the result of additional cuts + cut on phi mass peak 
				if( (kp+km).mass() < 1.050 && (kp+km).mass() > 1.0100 ){
				    
				    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_additional_pass_phi_mass")
				}

			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_fail_additional")
			    }

			}
												
		    }

		}
	    }

 	    // This is a good phi event, let's see if this is simulation.
	    // If so, find the generated partciles and make more histograms.
	    // Fill out simulation histograms for all generated events.
	    // Secondly, fill out histograms for when a good phi event is reconstructed.
	    // we can later divide the reconstructed number of events by the generated for acceptance calculations
	    // and other items

	    (0..<event.mc_npart).find { j ->
		event.mc_pid[j] == 11
	    }?.each { gidx ->
		def gen_ele = LorentzVector.withPID(11, event.mc_px[gidx], event.mc_py[gidx], event.mc_pz[gidx])
		//def gen_sector = (int) getGeneratedSector(Math.toDegrees(gen_ele.phi()))				

		(1..<event.mc_npart).find { k ->
		    event.mc_pid[k] == 2212
		}?.each { pidx ->
		    def gen_pro = LorentzVector.withPID(2212, event.mc_px[pidx], event.mc_py[pidx], event.mc_pz[pidx])

 		    (2..<event.mc_npart).find { k ->
			event.mc_pid[k] == 321
		    }?.each { kpidx ->
 			def gen_kp = LorentzVector.withPID(321, event.mc_px[kpidx], event.mc_py[kpidx], event.mc_pz[kpidx])

			(3..<event.mc_npart).find { k ->
			    event.mc_pid[k] == -321
			}?.each { kmidx ->
 			    def gen_km = LorentzVector.withPID(-321, event.mc_px[kmidx], event.mc_py[kmidx], event.mc_pz[kmidx])

			    if( reconstructed_event_present ){
				fillSimulationResolutions(histos, event, rec_ele, rec_pro, rec_kp, rec_km, gen_ele, gen_pro, gen_kp, gen_km, 'sim_resolutions' )
				fillSimulationHistograms(histos, beam, target, rec_ele, rec_pro, rec_kp, rec_km, 'sim_rec')
				//use this histogram to compare to the rec and generated final state ones filled by the method above
				fillSimulationHistograms(histos, beam, target, gen_ele, gen_pro, gen_kp, gen_km, 'sim_gen_per_rec')			  
			    }

			    fillSimulationHistograms(histos, beam, target, gen_ele, gen_pro, gen_kp, gen_km, 'sim_gen')
			    //fillEventHistogprams(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, '')

			    if( rec_electron_present ){
				fillElectronAccp(event, histos, gen_ele, 'sim_gen_ele_present')
			    }
			    if( rec_proton_present ){
				fillProtonAccp(event, histos, gen_pro, 'sim_gen_pro_present')
			    }
			    if( rec_kaonP_present ){
				fillKaonPAccp(event, histos, gen_kp, 'sim_gen_kp_present')
			    }
			    if( rec_kaonM_present ){
				fillKaonMAccp(event, histos, gen_km, 'sim_gen_km_present')
			    }
			    
			    
			}
		    }
		}
	    }   
	    cc+=1
	}
	
	reader.close()
    }
}

println('Finished Event Loop')

def outname = args[0].split('/')[-1].tokenize( '.' )[0]
println('appending ' + outname + ' to outputfile name ' )

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("monitor.hipo")

def root_out = new ROOTFile("monitor_phi_${outname}_pidtype1-accp-inb.root")
histos.each{root_out.writeDataSet(it.value)}
root_out.close()

//println(' Writing out txt file if needed')
if( f_out != null ){
    t1.write()
    f_out.close()
}

