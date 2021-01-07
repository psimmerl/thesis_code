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



def field_type = args[0]//"outb"
args.remove(0)

def gen_asy = -0.12//-0.07282


def Random r = new Random(); 

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
//def base_path_cut = "/u/home/psimmerl/thesis_code/projects/exclusive_phi/epkpkm_top/"
def epkpkmxe_cut = loadCuts(new File(base_path_cut+"excl_phi_me_limits_pidtype1_"+field_type+".txt"),4,4)
def ekpkmX_cut = loadCuts(new File(base_path_cut+"excl_phi_mm2_pro_limits_pidtype1_"+field_type+".txt"),4,4)
def epkmX_cut = loadCuts(new File(base_path_cut+"excl_phi_mm2_kp_limits_pidtype1_"+field_type+".txt"),4,4)
def epkpX_cut = loadCuts(new File(base_path_cut+"excl_phi_mm2_km_limits_pidtype1_"+field_type+".txt"),4,4)


def epkpkmXe_cut_range = [1 : loadCuts(new File(base_path_cut+"excl_phi_me_limits_pidtype1_"+field_type+".txt"), 1, 1),
			  2 : loadCuts(new File(base_path_cut+"excl_phi_me_limits_pidtype1_"+field_type+".txt"), 2, 2),
			  3 : loadCuts(new File(base_path_cut+"excl_phi_me_limits_pidtype1_"+field_type+".txt"), 3, 3),
			  4 : loadCuts(new File(base_path_cut+"excl_phi_me_limits_pidtype1_"+field_type+".txt"), 4, 4)
			  ]

def ekpkmX_cut_range = [1 : loadCuts(new File(base_path_cut+"excl_phi_mm2_pro_limits_pidtype1_"+field_type+".txt"), 1, 1),
			2 : loadCuts(new File(base_path_cut+"excl_phi_mm2_pro_limits_pidtype1_"+field_type+".txt"), 2, 2),
			3 : loadCuts(new File(base_path_cut+"excl_phi_mm2_pro_limits_pidtype1_"+field_type+".txt"), 3, 3),
			4 : loadCuts(new File(base_path_cut+"excl_phi_mm2_pro_limits_pidtype1_"+field_type+".txt"), 4, 4)
                        ]

def epkmX_cut_range = [1 : loadCuts(new File(base_path_cut+"excl_phi_mm2_kp_limits_pidtype1_"+field_type+".txt"), 1, 1),
		       2 : loadCuts(new File(base_path_cut+"excl_phi_mm2_kp_limits_pidtype1_"+field_type+".txt"), 2, 2),
		       3 : loadCuts(new File(base_path_cut+"excl_phi_mm2_kp_limits_pidtype1_"+field_type+".txt"), 3, 3),
		       4 : loadCuts(new File(base_path_cut+"excl_phi_mm2_kp_limits_pidtype1_"+field_type+".txt"), 4, 4)
                       ]

def epkpX_cut_range = [1 : loadCuts(new File(base_path_cut+"excl_phi_mm2_km_limits_pidtype1_"+field_type+".txt"), 1, 1),
		       2 : loadCuts(new File(base_path_cut+"excl_phi_mm2_km_limits_pidtype1_"+field_type+".txt"), 2, 2),
		       3 : loadCuts(new File(base_path_cut+"excl_phi_mm2_km_limits_pidtype1_"+field_type+".txt"), 3, 3),
		       4 : loadCuts(new File(base_path_cut+"excl_phi_mm2_km_limits_pidtype1_"+field_type+".txt"), 4, 4)
                       ]

println(epkpkmXe_cut_range)
println(ekpkmX_cut_range)

println(' Event Cut parameters ')
println(' --> me ' + epkpkmxe_cut)
println(' --> mm2 pr ' + ekpkmX_cut)
println(' --> mm2 kp ' + epkmX_cut)
println(' --> mm2 km ' + epkpX_cut)


excl_cuts_inb = [
    epkpkmxe  : [0.08 - 0.0398*sig, 0.08 + 0.0398*sig], //epkpkmxe_cut,   
    epkpX     : [0.248 - 0.059*sig, 0.248 + 0.059*sig], //epkpX_cut,   
    epkmX     : [0.248 - 0.055*sig, 0.248 + 0.055*sig], //epkmX_cut,
    ekpkmX    : [0.9431 - 0.0719*sig, 0.9431 + 0.0719*sig], //ekpkmX_cut,
    el_p_min  : 0.5,
    pro_p_max : 3.5,
    kp_p_max  : 3.5,
    km_p_max  : 3.5,
    pt_max    : 0.12,
    q2_min    : 1,
    w_min     : 2,
    kpkm_mass : [1.0107, 1.0287],
    delta_vz_kpkm : [-7, 5],
    pro_chi2 : [0,6],
    kp_chi2 : [0,6],
    km_chi2 : [0,6],
    cpl_pro_max : [0, 9],
    cpl_kp_max : [0, 9],
    cpl_km_max : [0, 9]
]

excl_cuts_outb = [
    epkpkmxe  : epkpkmxe_cut, 
    epkpX     : epkpX_cut,    
    epkmX     : epkmX_cut,    
    ekpkmX    : ekpkmX_cut,   
    el_p_min  : 0.5,
    pro_p_max : 3.5,
    kp_p_max  : 3.5,
    km_p_max  : 3.5,
    pt_max    : 0.12,
    q2_min    : 1,
    w_min     : 2,
    kpkm_mass : [1.0110, 1.0281],
    delta_vz_kpkm : [-6, 9],
    pro_chi2 : [0,6],
    kp_chi2 : [0,6],
    km_chi2 : [0,6],
    cpl_pro_max : [0, 9],
    cpl_kp_max : [0, 9],
    cpl_km_max : [0, 9]
]



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
    phi_gen   : [-180, 180],
    angle_ep  : [160, 180],
    angle_ep_diff : [160, 200],
    q2        : [1.2, 4.5],
    vz        : [-20, 15],
    dvz        : [-20, 20],
    missing_mass_small: [-0.3, 0.3],
    prot_chi2 : [-10, 10],
    prot_chi2_small : [-15, 15],
    beta : [0, 1.2],
    
    mm_range1 : [-0.05, 0.05],
    mm_range2 : [0.0, 2.0],
    mm_range3 : [-0.5, 1.5],

    q2 : [0, 12],
    xb : [0, 1.15],
    phi_trento : [-180, 180],
    minus_t : [0, 6],
    im_phi : [0, 2.5],
    im_pkm : [1.3, 3.0],
    cone_angle: [0,8],
    missing_energy : [-1, 1],
    combs: [1,10],
    coplan_ang : [0, 10],
    delta_theta : [-10, 10],
    pt : [0, 5],
    
    mm2_small:[-0.15, 0.15],
    mm2: [-1.5, 1.5],
    mm2_big : [0, 5 ],    
    mass:[0.8, 1.5],
    mass_small:[0.95, 1.12],
    prokm: [1, 5.5],
    ekp: [1, 5],
    kpkm_big: [0.75, 2.2],
    kpkm_small : [0.8, 1.5],
    mm_ep: [0, 5],
    epkp_big : [0.0, 0.7],
    mm_epkm_big:[0.0, 0.7],
    pt_big: [0, 0.25],
    missing_energy_big: [-0.2, 0.5],
    prkm_big : [0.5, 5],

    epkpkm_big : [-0.015, 0.015],

    res_p_tiny : [-0.1, 0.1],
    res_p_small : [-0.2,0.2],
    res_p_med : [-0.5,0.5],
    res_p_big : [-1, 1],
    res_p_extra : [-3,3],

    res_theta : [-0.2,0.2],
    res_phi : [-4,4],
    res_phys : [-0.2,0.2],

    cos_cm : [-1, 1]   
]

lim = tighter_kin_bounds

def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}

my_phi_trento = new H1F('my_phi_trento','my_phi_trento', 10, -180.0, 180.0)


histos = new ConcurrentHashMap()
histoBuilders = [
    p_ele      : { title -> limited_h1(title, 200, lim.p_ele) },
    theta_ele  : { title -> limited_h1(title, 200, lim.theta_ele) },
    phi_ele    : { title -> limited_h1(title, 200, lim.phi_ele) },

    p_pro      : { title -> limited_h1(title, 200, lim.p_pro) },
    theta_pro  : { title -> limited_h1(title, 200, lim.theta_pro) },
    phi_pro    : { title -> limited_h1(title, 200, lim.phi_pro) },

    p_kaon      : { title -> limited_h1(title, 200, lim.p_kaon) },
    theta_kaon  : { title -> limited_h1(title, 200, lim.theta_kaon) },
    phi_kaon    : { title -> limited_h1(title, 200, lim.phi_kaon) },

    w        : { title -> limited_h1(title, 200, lim.w) },
    w_val    : { title -> limited_h1(title, 200, lim.w_val) },
    p_ele    : { title -> limited_h1(title, 200, lim.p_ele) },
    p_pro    : { title -> limited_h1(title, 200, lim.p_pro) },
    vz       : { title -> limited_h1(title, 200, lim.vz) },
    angle_ep : { title -> limited_h1(title, 200, lim.angle_ep) },
    theta_p  : { title -> limited_h1(title, 200, lim.theta_pro_ftof) },
    theta_ele: { title -> limited_h1(title, 200, lim.theta_ele) },
    pro_chi2 : { title -> limited_h1(title, 200, lim.prot_chi2_small) },
    missing_mass_small: { title -> limited_h1(title, 200, lim.missing_mass_small) }, // added by brandon c . to validate
    combs    : { title -> limited_h1(title, 200, lim.combs) },
    coplan_angle :{ title -> limited_h1(title, 200, lim.coplan_ang) },
    delta_theta :{ title -> limited_h1(title, 200, lim.delta_theta) },    

    im_epkpkmX_range1 : { title -> limited_h1(title, 200, lim.mm_range1 ) },
    im_ekpkmX_range2 : { title -> limited_h1(title, 200, lim.mm_range2 ) },
    im_epkpX_range3 : { title -> limited_h1(title, 200, lim.mm_range3 ) },
    im_epkmX_range3 : { title -> limited_h1(title, 200, lim.mm_range3 ) },
        
    im_epkpkmX : { title -> limited_h1(title, 200, lim.mm2_small ) },
    im_ekpkmX  : { title -> limited_h1(title, 200, lim.mm2 ) },
    im_epkX : { title -> limited_h1(title, 200, lim.mm2 ) },
    im_ekpX : { title -> limited_h1(title, 200, lim.mm2_big ) },
    missing_e : { title -> limited_h1(title, 200, lim.mm2 ) },
    missing_p :{ title -> limited_h1(title, 200, lim.mm2 ) },
    im_pkm : { title -> limited_h1(title, 200, lim.im_pkm ) },
    im_kpkm : { title -> limited_h1(title, 150, lim.mass ) },
    im_kpkm_small : { title -> limited_h1(title, 300, lim.kpkm_small ) },
    phitrento_small: { title -> limited_h1(title, 150, lim.phi_trento ) },
    phitrento: { title -> limited_h1(title, 10, lim.phi_trento ) },
    phitrentoConfig2: { title -> limited_h1(title, 12, lim.phi_trento ) },

    res_p : { title -> limited_h1(title, 500, lim.res_p_small ) },
    res_theta : { title -> limited_h1(title, 500, lim.res_p_extra) },
    res_phi : { title -> limited_h1(title, 500, lim.res_phi ) },
    
    res_physics : { title -> limited_h1(title, 500, lim.res_phys ) },
    q2 : { title -> limited_h1(title, 500, lim.q2 ) },
    t  : { title -> limited_h1(title, 500, lim.minus_t ) },
    w  : { title -> limited_h1(title, 500, lim.w ) },
    xb : { title -> limited_h1(title, 500, lim.xb ) },
    pt : { title -> limited_h1(title, 500, lim.pt ) },
 
    cm_cos : { title -> limited_h1(title, 10, lim.cos_cm ) },

    res_q2 : { title-> limited_h1(title, 200, lim.res_p_small )},
    res_xb : { title-> limited_h1(title, 200, lim.res_p_tiny )},
    res_t : { title-> limited_h1(title, 200,  lim.res_p_small )},
    res_trento : { title-> limited_h1(title, 200, lim.res_p_extra )},

    dvz   : { title -> limited_h1(title, 200, lim.dvz)},
    
]

histoBuilders2 = [
    p_ele_theta    : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.theta_ele) },
    p_pro_theta    : { title -> limited_h2(title, 200, 200, lim.p_pro, lim.theta_pro_ftof) },
    p_kp_theta     : { title -> limited_h2(title, 200, 200, lim.p_kp, lim.theta_kp) },
    p_km_theta     : { title -> limited_h2(title, 200, 200, lim.p_km, lim.theta_km) },
    p_vz           : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.vz)},
    p_dvz          : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.dvz)},
		      
    phi_theta_ele : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_ele) },
    phi_theta_ele_gen : { title -> limited_h2(title, 100, 100, lim.phi_gen, lim.theta_ele) },
    phi_theta_proton : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_pro_ftof) },
    phi_theta_kp : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_kp) },
    phi_theta_km : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_km) },
    phi_theta_proton_gen : { title -> limited_h2(title, 100, 100, lim.phi_gen, lim.theta_pro_ftof) },
    phi_theta_kp_gen : { title -> limited_h2(title, 100, 100, lim.phi_gen, lim.theta_kp) },
    phi_theta_km_gen : { title -> limited_h2(title, 100, 100, lim.phi_gen, lim.theta_km) },
    phi_vs_phi   : { title -> limited_h2(title, 200, 200, lim.phi_ele, lim.phi_ele) },
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
    phitrento_imkpkm : { title -> limited_h2(title, 10, 150, lim.phi_trento, lim.mass) }, 

    t_xb             : { title -> limited_h2(title, 200, 200, lim.minus_t, lim.xb) },
    t_w              : { title -> limited_h2(title, 200, 200, lim.minus_t, lim.w) },
    mntm_t           : { title -> limited_h2(title, 200, 200, lim.p_kp, lim.minus_t) },
    theta_t           : { title -> limited_h2(title, 200, 200, lim.theta_kp, lim.minus_t) },
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

    // see impact of momentum cut on the kpkm mass
    p_imkpkm : { title -> limited_h2(title, 200, 200, lim.p_pro, lim.mass) },

    cm_cos_minut_t : {title->limited_h2(title, 200, 200, lim.minus_t, lim.cos_cm)},


    //mass resolution items
    res_pro_p : { title -> limited_h2(title, 200, 200, lim.p_pro, lim.res_p_big )},
    res_kp_p : { title -> limited_h2(title, 200, 200, lim.p_kp, lim.res_p_big )},
    res_km_p : { title -> limited_h2(title, 200, 200, lim.p_km, lim.res_p_big )},
 

    //resolution itesm for 2d
    res_p_delta_p     : {title -> limited_h2(title, 200, 200, lim.p_ele, lim.res_p_small) },
    res_theta_delta_p : {title -> limited_h2(title, 200, 200, lim.theta_ele, lim.res_p_small) },
    res_theta_delta_theta : {title -> limited_h2(title, 200, 200, lim.theta_ele, lim.res_p_extra) },
    
    res_xb_delta_xb : { title -> limited_h2(title, 200, 200, lim.xb, lim.res_p_tiny) },
    res_q2_delta_q2 : { title -> limited_h2(title, 200, 200, lim.q2, lim.res_p_small) },
    res_t_delta_t : { title -> limited_h2(title, 200, 200, lim.minus_t, lim.res_p_small) },
    

            
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

    def cplpro = (ekpkmX.vect().theta(pro.vect()))
    def cplkm = (epkpX.vect().theta(km.vect()))
    def cplkp = (epkmX.vect().theta(kp.vect()))
    
    def delta_theta_pro = Math.toDegrees(ekpkmX.theta() - pro.theta())
    def delta_theta_kp = Math.toDegrees(epkmX.theta() - kp.theta())
    def delta_theta_km = Math.toDegrees(epkpX.theta() - km.theta())
    

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
	    ekpkmX : ekpkmX, cplpro: cplpro, cplkm: cplkm, cplkp: cplkp, delta_theta_pro: delta_theta_pro,
	    delta_theta_kp: delta_theta_kp, delta_theta_km: delta_theta_km, pt: pt, lmda: lmda,
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


//feed in event, hist, 
// electron proton kaonP kaonM are maps
// ele pro kp km are lorentzvectors
def fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, comment ){

    histos.computeIfAbsent('p_ele_'+comment, histoBuilders.p_ele).fill(ele.p())
    histos.computeIfAbsent('theta_ele_'+comment, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_ele_'+comment, histoBuilders.phi_ele).fill(Math.toDegrees(ele.phi()))
    histos.computeIfAbsent('vz_ele_'+comment, histoBuilders.vz).fill( event.vz[electron.index])

    histos.computeIfAbsent('p_pro_'+comment, histoBuilders.p_pro).fill(pro.p())
    histos.computeIfAbsent('theta_pro_'+comment, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_pro_'+comment, histoBuilders.phi_pro).fill(Math.toDegrees(pro.phi()))
    histos.computeIfAbsent('vz_pro_'+comment, histoBuilders.vz).fill( event.vz[proton.index])

    histos.computeIfAbsent('p_kp_'+comment, histoBuilders.p_kaon).fill(kp.p())
    histos.computeIfAbsent('theta_kp_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_kp_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(kp.phi()))
    histos.computeIfAbsent('vz_kp_'+comment, histoBuilders.vz).fill( event.vz[kaonP.index])

    histos.computeIfAbsent('p_km_'+comment, histoBuilders.p_kaon).fill(km.p())
    histos.computeIfAbsent('theta_km_'+comment, histoBuilders.theta_kaon).fill(Math.toDegrees(km.theta()))
    histos.computeIfAbsent('phi_km_'+comment, histoBuilders.phi_kaon).fill(Math.toDegrees(km.phi()))
    histos.computeIfAbsent('vz_km_'+comment, histoBuilders.vz).fill( event.vz[kaonM.index])

    histos.computeIfAbsent('p_ele_theta_'+comment, histoBuilders2.p_ele_theta).fill(ele.p(), Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('p_pro_theta_'+comment, histoBuilders2.p_pro_theta).fill(pro.p(), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('p_kp_theta_'+comment, histoBuilders2.p_kp_theta).fill(kp.p(), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('p_km_theta_'+comment, histoBuilders2.p_km_theta).fill(km.p(), Math.toDegrees(km.theta()))
    
    histos.computeIfAbsent('phi_ele_theta_'+comment, histoBuilders2.phi_theta_ele).fill(shiftPhi(ele.phi()) , Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_pro_theta_'+comment, histoBuilders2.phi_theta_proton).fill(shiftPhi(pro.phi()), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_kp_theta_'+comment, histoBuilders2.phi_theta_kp).fill(shiftPhi(kp.phi()), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_km_theta_'+comment, histoBuilders2.phi_theta_km).fill(shiftPhi(km.phi()), Math.toDegrees(km.theta()))
    histos.computeIfAbsent('phi_vs_phi_'+comment, histoBuilders2.phi_vs_phi).fill( Math.toDegrees(pro.phi()), Math.toDegrees(ele.phi()))
    histos.computeIfAbsent('phi_trento_vs_phi_ele_'+comment, histoBuilders2.phi_vs_phi).fill( phykin.trentophi, Math.toDegrees(ele.phi()))
    histos.computeIfAbsent('phi_trento_vs_phi_pro_'+comment, histoBuilders2.phi_vs_phi).fill( phykin.trentophi, Math.toDegrees(pro.phi()))

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
    histos.computeIfAbsent('dvz_kpkm_'+comment, histoBuilders.dvz).fill((event.vz[kaonP.index] - event.vz[kaonM.index]))
    


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
    histos.computeIfAbsent('cpl_kp_'+comment,histoBuilders.coplan_angle).fill( phykin.cplkp )
    histos.computeIfAbsent('cpl_km_'+comment,histoBuilders.coplan_angle).fill( phykin.cplkm )
    
    // look at difference in cal and measured angle
    histos.computeIfAbsent('delta_theta_pro_'+comment,histoBuilders.delta_theta).fill( phykin.delta_theta_pro )
    histos.computeIfAbsent('delta_theta_kp_'+comment,histoBuilders.delta_theta).fill( phykin.delta_theta_kp )
    histos.computeIfAbsent('delta_theta_km_'+comment,histoBuilders.delta_theta).fill( phykin.delta_theta_km )   
    
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
    histos.computeIfAbsent('p_ele_t_'+comment,histoBuilders2.mntm_t).fill( ele.p(),  phykin.t)
    histos.computeIfAbsent('theta_pro_t_'+comment,histoBuilders2.theta_t).fill( Math.toDegrees(pro.theta()),  phykin.t)

    histos.computeIfAbsent('p_ekpkmX_t_'+comment,histoBuilders2.mntm_t).fill( phykin.ekpkmX.p(),  phykin.t)
    histos.computeIfAbsent('p_kp_q2_'+comment,histoBuilders2.mntm_q2).fill( kp.p(),  phykin.q2)

    histos.computeIfAbsent('cm_cos_theta_kp_'+comment,histoBuilders.cm_cos).fill(phykin.cos_theta_cm_kp)
    histos.computeIfAbsent('cm_cos_theta_kp_vs_t_'+comment,histoBuilders2.cm_cos_minut_t).fill(phykin.t, phykin.cos_theta_cm_kp)
    
    def helicity_state = helicity < 0 ? 'neg' : 'pos'
    histos.computeIfAbsent('phitrento_small_'+helicity_state+'_'+comment,histoBuilders.phitrento_small).fill(phykin.trentophi)
    histos.computeIfAbsent('phitrento_'+helicity_state+'_'+comment,histoBuilders.phitrento).fill(phykin.trentophi)    
    histos.computeIfAbsent('phitrento_config2_'+helicity_state+'_'+comment,histoBuilders.phitrentoConfig2).fill(phykin.trentophi)    
    
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
    histos.computeIfAbsent('mm_epkpkmX_'+comment,histoBuilders.im_epkpkmX).fill(phykin.epkpkmX.mass2())
    histos.computeIfAbsent('mm_ekpkmX_'+comment,histoBuilders.im_ekpkmX_range2).fill(phykin.ekpkmX.mass2())
    histos.computeIfAbsent('mm_epkpX_'+comment,histoBuilders.im_epkpX_range3).fill(phykin.epkpX.mass2())
    histos.computeIfAbsent('mm_epkmX_'+comment,histoBuilders.im_epkmX_range3).fill(phykin.epkmX.mass2())
    histos.computeIfAbsent('mm_ekpX_'+comment,histoBuilders.im_ekpX).fill(phykin.ekpX.mass2())
    histos.computeIfAbsent('missing_energy_'+comment,histoBuilders.missing_e).fill(phykin.epkpkmX.e())
    histos.computeIfAbsent('missing_p_'+comment,histoBuilders.missing_p).fill(phykin.pt)
    histos.computeIfAbsent('im_pkm_'+comment,histoBuilders.im_pkm).fill( (pro+km).mass() )
    histos.computeIfAbsent('im_pkp_'+comment,histoBuilders.im_pkm).fill( (pro+kp).mass() )
    histos.computeIfAbsent('im_kpkm_'+comment,histoBuilders.im_kpkm).fill( (kp+km).mass() )
    histos.computeIfAbsent('phitrento_imkpkm_'+comment,histoBuilders2.phitrento_imkpkm).fill( phykin.trentophi, (kp+km).mass())
    
    //mass info specialized plots
    histos.computeIfAbsent('mm_epkpkmX_mod_'+comment,histoBuilders.im_epkpkmX_range1).fill(phykin.epkpkmX.mass2())
    histos.computeIfAbsent('mm_ekpkmX_mod_'+comment,histoBuilders.im_ekpkmX_range2).fill(phykin.ekpkmX.mass2())
    histos.computeIfAbsent('mm_epkpX_mod_'+comment,histoBuilders.im_epkpX_range3).fill(phykin.epkpX.mass2())
    histos.computeIfAbsent('mm_epkmX_mod_'+comment,histoBuilders.im_epkmX_range3).fill(phykin.epkmX.mass2())
    

    //physics
    // 1d
    histos.computeIfAbsent('q2_'+comment, histoBuilders.q2).fill(phykin.q2)
    histos.computeIfAbsent('t_'+comment, histoBuilders.t).fill(phykin.t)
    histos.computeIfAbsent('w_'+comment, histoBuilders.w).fill(phykin.w)
    histos.computeIfAbsent('xb_'+comment, histoBuilders.xb).fill(phykin.xb)
    histos.computeIfAbsent('pt_'+comment, histoBuilders.pt).fill(phykin.pt)
    histos.computeIfAbsent('phitrento_'+comment, histoBuilders.phitrento_small).fill(phykin.trentophi)

    // check missing mass distributions in various phi trento bins
    def phi_trento_bin = my_phi_trento.getXaxis().getBin(phykin.trentophi) + 1
    histos.computeIfAbsent('mm_epkpkmX_trentoB'+phi_trento_bin+'_'+comment,histoBuilders.im_epkpkmX).fill(phykin.epkpkmX.mass2())
    histos.computeIfAbsent('mm_ekpkmX_trentoB'+phi_trento_bin+'_'+comment,histoBuilders.im_ekpkmX_range2).fill(phykin.ekpkmX.mass2())
    histos.computeIfAbsent('mm_epkpX_trentoB'+phi_trento_bin+'_'+comment,histoBuilders.im_epkpX_range3).fill(phykin.epkpX.mass2())
    histos.computeIfAbsent('mm_epkmX_trentoB'+phi_trento_bin+'_'+comment,histoBuilders.im_epkmX_range3).fill(phykin.epkmX.mass2())
    histos.computeIfAbsent('mm_ekpX_trentoB'+phi_trento_bin+'_'+comment,histoBuilders.im_ekpX).fill(phykin.ekpX.mass2())
    histos.computeIfAbsent('missing_energy_trentoB'+phi_trento_bin+'_'+comment,histoBuilders.missing_e).fill(phykin.epkpkmX.e())
    
    histos.computeIfAbsent('p_ele_theta_trentoB'+phi_trento_bin+'_'+comment, histoBuilders2.p_ele_theta).fill(ele.p(), Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('p_pro_theta_trentoB'+phi_trento_bin+'_'+comment, histoBuilders2.p_pro_theta).fill(pro.p(), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('p_kp_theta_trentoB'+phi_trento_bin+'_'+comment, histoBuilders2.p_kp_theta).fill(kp.p(), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('p_km_theta_trentoB'+phi_trento_bin+'_'+comment, histoBuilders2.p_km_theta).fill(km.p(), Math.toDegrees(km.theta()))
            
}

def fillEventSpecial(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, comment ){
    def helicity_state = helicity < 0 ? 'neg' : 'pos'
    //def phi_trento_bin = my_phi_trento.getXaxis().getBin(phykin.trentophi) + 1
    //println(' phi trento bin ' + phi_trento_bin + ' angle ' +  phykin.trentophi)
    histos.computeIfAbsent('phitrento_'+helicity_state+'_'+comment,histoBuilders.phitrento).fill(phykin.trentophi)    
    histos.computeIfAbsent('phitrento_config2_'+helicity_state+'_'+comment,histoBuilders.phitrentoConfig2).fill(phykin.trentophi)    
    //histos.computeIfAbsent('im_kpkm_trento_bin'+phi_trento_bin+'_'+helicity_state+'_'+comment, histoBuilders.im_kpkm).fill( (kp+km).mass() )
        

}

//fill data resoltuion through the missing particle 

def fillDataResolutions( event, histos, ele, pro, kp, km, phykin, comment ){

    histos.computeIfAbsent('res_ekpkmX_mntm_vs_p_'+comment, histoBuilders2.res_pro_p).fill(pro.p(), -( phykin.ekpkmX.p() - pro.p()))
    histos.computeIfAbsent('res_epkpX_mntm_vs_p_'+comment, histoBuilders2.res_km_p).fill( km.p(), -(phykin.epkpX.p() - km.p() ))
    histos.computeIfAbsent('res_epkmX_mntm_vs_p_'+comment, histoBuilders2.res_kp_p).fill( kp.p(), -(phykin.epkmX.p() - kp.p() ))

}

//fill when event for phi is present in simulation 
def fillSimulationResolutions(histos, event, ele, pro, kp, km, gen_ele, gen_pro, gen_kp, gen_km, phykin, sim_phykin, comment ){
    
   // electron res
    histos.computeIfAbsent('p_ele_res_'+comment, histoBuilders.res_p).fill(ele.p() - gen_ele.p())
    histos.computeIfAbsent('p_ele_delta_p_'+comment, histoBuilders2.res_p_delta_p).fill(ele.p(), ele.p() - gen_ele.p())
    histos.computeIfAbsent('theta_ele_delta_p_'+comment, histoBuilders2.res_theta_delta_p).fill(Math.toDegrees(ele.theta()), ele.p() - gen_ele.p())
    histos.computeIfAbsent('theta_ele_delta_theta_'+comment, histoBuilders2.res_theta_delta_theta).fill(Math.toDegrees(ele.theta()), Math.toDegrees(ele.theta() - gen_ele.theta()))

    histos.computeIfAbsent('theta_ele_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(ele.theta() - gen_ele.theta()))
    histos.computeIfAbsent('phi_ele_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(ele.phi() - gen_ele.phi()))

    //proton resolution 2d
    histos.computeIfAbsent('p_pro_delta_p_'+comment, histoBuilders2.res_p_delta_p).fill(pro.p(), pro.p() - gen_pro.p())
    histos.computeIfAbsent('theta_pro_delta_p_'+comment, histoBuilders2.res_theta_delta_p).fill(Math.toDegrees(pro.theta()), pro.p() - gen_pro.p())
    histos.computeIfAbsent('theta_pro_delta_theta_'+comment, histoBuilders2.res_theta_delta_theta).fill(Math.toDegrees(pro.theta()), Math.toDegrees(pro.theta() - gen_pro.theta()))
    
    // 1d
    histos.computeIfAbsent('p_pro_res_'+comment, histoBuilders.res_p).fill(pro.p() - gen_pro.p())
    histos.computeIfAbsent('theta_pro_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(pro.theta() - gen_pro.theta()))
    histos.computeIfAbsent('phi_pro_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(pro.phi() - gen_pro.phi()))

    // kaonP 2d resolutions
    histos.computeIfAbsent('p_kp_delta_p_'+comment, histoBuilders2.res_p_delta_p).fill(kp.p(), kp.p() - gen_kp.p())
    histos.computeIfAbsent('theta_kp_delta_p_'+comment, histoBuilders2.res_theta_delta_p).fill(Math.toDegrees(kp.theta()), kp.p() - gen_kp.p())
    histos.computeIfAbsent('theta_kp_delta_theta_'+comment, histoBuilders2.res_theta_delta_theta).fill(Math.toDegrees(kp.theta()), Math.toDegrees(kp.theta() - gen_kp.theta()))
    
    // 1d
    histos.computeIfAbsent('p_kp_res_'+comment, histoBuilders.res_p).fill(kp.p() - gen_kp.p())
    histos.computeIfAbsent('theta_kp_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(kp.theta() - gen_kp.theta()))
    histos.computeIfAbsent('phi_kp_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(kp.phi() - gen_kp.phi()))
    
    // kaonM 2d resolutions
    histos.computeIfAbsent('p_km_delta_p_'+comment, histoBuilders2.res_p_delta_p).fill(km.p(), km.p() - gen_km.p())
    histos.computeIfAbsent('theta_km_delta_p_'+comment, histoBuilders2.res_theta_delta_p).fill(Math.toDegrees(km.theta()), km.p() - gen_km.p())
    histos.computeIfAbsent('theta_km_delta_theta_'+comment, histoBuilders2.res_theta_delta_theta).fill(Math.toDegrees(km.theta()), Math.toDegrees(km.theta() - gen_km.theta()))
    
    histos.computeIfAbsent('p_km_res_'+comment, histoBuilders.res_p).fill(km.p() - gen_km.p())
    histos.computeIfAbsent('theta_km_res_'+comment, histoBuilders.res_theta).fill(Math.toDegrees(km.theta() - gen_km.theta()))
    histos.computeIfAbsent('phi_km_res_'+comment, histoBuilders.res_phi).fill(Math.toDegrees(km.phi() - gen_km.phi()))    

    histos.computeIfAbsent('delta_q2_'+comment, histoBuilders.res_q2).fill( phykin.q2 - sim_phykin.q2)
    histos.computeIfAbsent('delta_xb_'+comment, histoBuilders.res_xb).fill( phykin.xb - sim_phykin.xb)
    histos.computeIfAbsent('delta_t_'+comment, histoBuilders.res_t).fill( phykin.t - sim_phykin.t)
    histos.computeIfAbsent('delta_phi_trento_'+comment, histoBuilders.res_trento).fill( phykin.trentophi - sim_phykin.trentophi)
    
    histos.computeIfAbsent('q2_delta_q2_'+comment, histoBuilders2.res_q2_delta_q2).fill(phykin.q2, phykin.q2 - sim_phykin.q2)
    histos.computeIfAbsent('xb_delta_xb_'+comment, histoBuilders2.res_xb_delta_xb).fill(phykin.xb, phykin.xb - sim_phykin.xb)
    histos.computeIfAbsent('t_delta_t_'+comment, histoBuilders2.res_t_delta_t).fill(phykin.t, phykin.t - sim_phykin.t)
    
    
}

def fillSimulationHistograms(histos, beam, target, ele, pro, kp, km, helicity, comment){

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

    histos.computeIfAbsent('p_ele_theta_'+comment, histoBuilders2.p_ele_theta).fill(ele.p(), Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('p_pro_theta_'+comment, histoBuilders2.p_pro_theta).fill(pro.p(), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('p_kp_theta_'+comment, histoBuilders2.p_kp_theta).fill(kp.p(), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('p_km_theta_'+comment, histoBuilders2.p_km_theta).fill(km.p(), Math.toDegrees(km.theta()))
    
    histos.computeIfAbsent('phi_ele_theta_'+comment, histoBuilders2.phi_theta_ele_gen).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta()))
    histos.computeIfAbsent('phi_pro_theta_'+comment, histoBuilders2.phi_theta_proton_gen).fill(Math.toDegrees(pro.phi()), Math.toDegrees(pro.theta()))
    histos.computeIfAbsent('phi_kp_theta_'+comment, histoBuilders2.phi_theta_kp_gen).fill(Math.toDegrees(kp.phi()), Math.toDegrees(kp.theta()))
    histos.computeIfAbsent('phi_km_theta_'+comment, histoBuilders2.phi_theta_km_gen).fill(Math.toDegrees(km.phi()), Math.toDegrees(km.theta()))
    histos.computeIfAbsent('phi_vs_phi_'+comment, histoBuilders2.phi_vs_phi).fill(Math.toDegrees(ele.phi()),  Math.toDegrees(pro.phi()) )
    
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

    def helicity_state = helicity < 0 ? 'neg' : 'pos'
    histos.computeIfAbsent('phitrento_'+helicity_state+'_'+comment,histoBuilders.phitrento).fill(phykin.trentophi)
    histos.computeIfAbsent('phitrento_config2_'+helicity_state+'_'+comment,histoBuilders.phitrentoConfig2).fill(phykin.trentophi)    

}


def smearParticleVector( V4, r, pid){//, pscale, thetascale, phiscale){
    // original from fx are 2, 2.5, 3.5 
    def inM = V4.mass();
    //true generated values
    def sP  = V4.p();
    def sTh = V4.theta();
    def sPh = V4.phi();


    //calculate resolutions
    def sThD = Math.toDegrees(sTh)
    def momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
    def momS2 = 0.02*sThD;
    def momR  = 0.01 * Math.sqrt( Math.pow(momS1*sP,2) + Math.pow(momS2,2) );
    momR *= 2.0//pscale;

    def theS1 = 0.004*sThD + 0.1;
    def theS2 = 0;
    def theR  = Math.sqrt(Math.pow(theS1*Math.sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + Math.pow(theS2,2) );
    theR *= 2.5//thetascale;

    def phiS1 = 0.85-0.015*sThD;
    def phiS2 = 0.17-0.003*sThD;
    def phiR  = Math.sqrt(Math.pow(phiS1*Math.sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + Math.pow(phiS2,2) );
    phiR *= 3.5//phiscale;

    //overwrite EB
    sPh += Math.toRadians( phiR ) * r.nextGaussian() 
    sTh += Math.toRadians( theR ) * r.nextGaussian()
    sP  += momR  * r.nextGaussian() *  V4.p() ;
    // EB_rec_mom = GEN_mom + resolution_momentum x gaussian x GEN_mom
    // EB_rec_ang = GEN_ang + resolution_angle x gaussian 

    def px_smear = sP * Math.sin(sTh) * Math.cos(sPh)
    def py_smear = sP * Math.sin(sTh) * Math.sin(sPh)
    def pz_smear = sP * Math.cos(sTh)

    def v4 = LorentzVector.withPID(pid, px_smear, py_smear, pz_smear)
    return v4
    
}


def generateFakeHelicity( gen_asy, gen_phykin, gen_phykin_smear, r ){
    def gen_hel = -1
    def gen_hel_smear = -1

    def gen_bsa = gen_asy*Math.sin( Math.toRadians(gen_phykin.trentophi) )
    def gen_bsa_smear = gen_asy*Math.sin( Math.toRadians(gen_phykin_smear.trentophi) )
    //sample for uniform to get the helicity, use same value for with/ without smear
    def par_to_pass = r.nextDouble()

    if( par_to_pass < (1+gen_bsa)/2.0 ){
	gen_hel=1
    }
    if( par_to_pass < (1+gen_bsa_smear)/2.0 ){
	gen_hel_smear=1
    }


    return [gen_hel, gen_hel_smear]
}

def createFileOut( f_out ){
    return new File(f_out)
}

def writeOutEvent( t_out, temp_l ){
    //t_out.fill(temp_l[0], temp_l[1], temp_l[2])   

    t_out.fill(*temp_l)//.collect{partbank.getFloat(it, eleind)}
    //t_out.fill(temp_l)   
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


//// fix field_type based on run number
//def s_start_run = ((args[0].split('/')[-1].tokenize( '.' )[0]).split('_')[-1])
//def start_run = (s_start_run.replaceFirst("^0*", "")).toInteger()


//println(' Start run number is '+start_run)
if( field_type == "inb" ){
    excl_cuts=excl_cuts_inb
}
else if( field_type == "outb"){
    excl_cuts=excl_cuts_outb   
}
println(' old cut values me ' + excl_cuts.epkpkmxe)
println(' old cut values mm2 pro ' + excl_cuts.ekpkmX)
println(' old cut values mm2 kp ' + excl_cuts.epkmX)
println(' old cut values mm2 km ' + excl_cuts.epkpX)

//note that pidtype 3 is ALL RGA analysis note cuts are used
// note pidtype4a is all rga cuts - found correction needed to cal sf

def write_to_file = true
def f_out = null
def t1=null
if( write_to_file ){
    f_out = new ROOTFile("phi_mass_results_pidtype1_cutvar1_rmvres12-${field_type}_16core_incl.root")
    t1 = f_out.makeNtuple('clas12','data','phimass:helicity:trento:q2:xb:t:w:run')
}


def topology = 'epkpkm'
def field_setting = 'inbending'
if ( field_type  == 'outb' ){
    field_setting = 'outbending'
}
// cut lvl meanings: 0 loose, 1 med, 2 tight
el_cut_strictness_lvl=["ecal_cut_lvl":1,
		       "nphe_cut_lvl":1,
		       "vz_cut_lvl":1,
		       "min_u_cut_lvl":1,
		       "min_v_cut_lvl":0,
		       "min_w_cut_lvl":0,
		       "max_u_cut_lvl":1,
		       "max_v_cut_lvl":0,
		       "max_w_cut_lvl":0,
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
    electron_from_event.passElectronVertexCut, // sector independent but field dependent
    electron_from_event.passElectronPCALFiducialCut, 
    electron_from_event.passElectronAntiPionCut, // save as rga note
    electron_from_event.passElectronEIEOCut, // technically pcal energy cut at 0.07
    electron_from_event.passElectronSamplingFractionCut,
    electron_from_event.passElectronDCFiducialCuts
]

def proton_cand = new HadronID(hadron_id:2212)
def kaonP_cand = new HadronID(hadron_id:321)
def kaonM_cand = new HadronID(hadron_id:-321)

def myProtonCuts = [
    proton_cand.passHadronEBPIDCut,    
    proton_cand.passTrajChi2Cut,
    proton_cand.passDCFiducialCutChi2Region1,
    proton_cand.passDCFiducialCutChi2Region2,
    proton_cand.passDCFiducialCutChi2Region3
]

def myKaonPCuts = [
    kaonP_cand.passHadronEBPIDCut,
    kaonP_cand.passTrajChi2Cut,
    kaonP_cand.passDCFiducialCutChi2Region1,
    kaonP_cand.passDCFiducialCutChi2Region2,
    kaonP_cand.passDCFiducialCutChi2Region3
]

def myKaonMCuts = [    
    kaonM_cand.passHadronEBPIDCut,
    kaonM_cand.passTrajChi2Cut,
    kaonM_cand.passDCFiducialCutChi2Region1,
    kaonM_cand.passDCFiducialCutChi2Region2,
    kaonM_cand.passDCFiducialCutChi2Region3
]

def me_all_counter = 0
def me_pass_counter = 0
def me_fail_counter = 0

GParsPool.withPool 8, {
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

	    // boolean flag used to get resolution with simulations
	    def reconstructed_event_present = false
	    def rec_ele, rec_pro, rec_kp, rec_km, phykin = null
	    def gen_ele_smear, gen_pro_smear, gen_kp_smear, gen_km_smear
	    def helicity, gen_helicity, gen_helicity_smear = null
	    def gen_phykin, gen_phykin_with_smear = null
	    
	    //create the generated lorentz vector particles here
	    // this block is specific to smearing the simulated LV to better match data
 	    if( event.mc_npart > 0 ){
		gen_ele_smear = LorentzVector.withPID(11, event.mc_px[0], event.mc_py[0], event.mc_pz[0])
		gen_pro_smear = LorentzVector.withPID(2212, event.mc_px[1], event.mc_py[1], event.mc_pz[1])
		gen_kp_smear = LorentzVector.withPID(321, event.mc_px[2], event.mc_py[2], event.mc_pz[2])
		gen_km_smear = LorentzVector.withPID(-321, event.mc_px[3], event.mc_py[3], event.mc_pz[3])

		gen_phykin =  getPhysicsKin(beam, target, gen_ele_smear, gen_pro_smear, gen_kp_smear, gen_km_smear)		
		//smear the generated particles for when the reconstructed ones pass all
 		gen_ele_smear = smearParticleVector(gen_ele_smear, r, 11)//, 2.0, 2.5, 3.5)
		gen_pro_smear = smearParticleVector(gen_pro_smear, r, 2212)//, 2.0, 2.5, 3.5)
		gen_kp_smear = smearParticleVector(gen_kp_smear, r, 321)//, 2.0, 2.5, 3.5)
		gen_km_smear = smearParticleVector(gen_km_smear, r, -321)//, 2.0, 2.5, 3.5)
		gen_phykin_with_smear =  getPhysicsKin(beam, target, gen_ele_smear, gen_pro_smear, gen_kp_smear, gen_km_smear)		

		(gen_helicity, gen_helicity_smear) = generateFakeHelicity( gen_asy, gen_phykin, gen_phykin_with_smear, r )
		helicity = gen_helicity
	    }		
	    else{
		helicity = event.helicity
	    }
	     
		
	    
	    //#################################################################################################################################
	    // data analysis code

	    //get electron proton pair first
	    def ep = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
		def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    

		return [index:indx, lv:ele, name:'el', status:event.status[indx] ]
		//return ele
	    }.collectMany{ mele -> 
		//(0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		(0..<event.npart).findAll{event.pid[it]==2212 && !myProtonCuts.collect{ pr_test -> pr_test(event,it) }.contains(false)}.findResults{ipro->
		    def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 	
		    def mpro = [index:ipro, lv:pro, name:'pro', status:event.status[ipro]]
		    def ctof = (event.pid[ipro] >= 2000 && event.pid[ipro] < 4000) ? "FTOF" : "CTOF"
		    def pr_fid_stat = myProtonCuts.collect{ pr_test -> pr_test(event,ipro) }
		    fillProtonID(event, histos, mele.lv, pro, ipro, 'eb_'+ctof)
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
			fillKaonPlusID(event, histos, mele.lv, kp, ikp, 'eb_'+ctof)
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
			fillKaonMinusID(event, histos, mele.lv, km, ikm, 'eb_'+ctof)
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
		def pass_delta_vz_el_hadron = false
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
			
			// cuts on the vertex position - from rga analysis note august 2020
			// easier to do delta vz here - > difference in electron and proton, kp and km			
			def pass_delta_vz_el_pr = Math.abs( event.vz[el.index] - event.vz[pro.index] ) < 20
			def pass_delta_vz_el_kp = Math.abs( event.vz[el.index] - event.vz[kp.index] ) < 20
			def pass_delta_vz_el_km = Math.abs( event.vz[el.index] - event.vz[km.index] ) < 20
		 	pass_delta_vz_el_hadron = (pass_delta_vz_el_pr && pass_delta_vz_el_kp && pass_delta_vz_el_km)
			
			if( kp_ctof == "FTOF" ){
			    check_for_duplicates=true
			    kp_ftof_p = kp.lv.p()
			}
			if( combs.size() == 2 && kp_ctof == "CTOF" && check_for_duplicates ){
			    //println('combs size ' + combs.size() + ' kp ftof p ' + kp_ftof_p + ' ' + kp.lv.p())
			}
			    

			//def phykin =  getPhysicsKin(beam, target, el.lv, pro.lv, kp.lv, km.lv)
			
			if(el_ctof == "FTOF" && pr_ctof == "FTOF" && 
			   kp_ctof == "FTOF" && km_ctof == "FTOF" ){
			    all_in_fd=true			    
			}		
		    }
		    if( all_in_fd ) { histos.computeIfAbsent('combsize_all_in_fd',histoBuilders.combs).fill(combs.size()) }

		    if( combs.size() == 1 && all_in_fd && pass_delta_vz_el_hadron ){
			def (electron, proton, kaonP, kaonM ) = combs[0]
			def ele = electron.lv
			def pro = proton.lv
			def kp = kaonP.lv
			def km = kaonM.lv
			phykin =  getPhysicsKin(beam, target, ele, pro, kp, km)
			
			def comment='combssize1'
			//kinematics for all events with only one possible combination of proton, kp, km. i.e. phi candidate
			fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity,  phykin, comment)
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

			//coplanarity cuts
			def pass_cpl_pro = phykin.cplpro < excl_cuts.cpl_pro_max[1]
			def pass_cpl_kp = phykin.cplkp < excl_cuts.cpl_kp_max[1]
			def pass_cpl_km = phykin.cplkm < excl_cuts.cpl_km_max[1]

			//delta theta cuts ( calculated - measured )
			def pass_delta_theta_pro = Math.abs(phykin.delta_theta_pro) < 6
			def pass_delta_theta_kp = Math.abs(phykin.delta_theta_kp) < 6
			def pass_delta_theta_km = Math.abs(phykin.delta_theta_km) < 6

			//make plots for passing all but one of the missing mass/energy cuts
			// see ekpkmX first - missing proton
			if( pass_epkpkmXe && pass_epkpX && pass_epkmX ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_but_missing_proton')
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'fail_at_least_missing_proton')
			}
			// see  missnig kp
			if( pass_epkpkmXe && pass_epkpX && pass_ekpkmX  ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_but_missing_kaonP')			    
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'fail_at_least_missing_kaonP')
			}
			// see  missnig km
			if( pass_epkpkmXe && pass_epkmX && pass_ekpkmX  ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_but_missing_kaonM')
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'fail_at_least_missing_kaonM')
			}
			// see missing energy
			if(  pass_epkpX && pass_epkmX && pass_ekpkmX  ){
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_but_missing_energy')
			}
			else{
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'fail_at_least_missing_energy')
			}

			// pass or fails only missing energy cut --> special criteria for joo
			me_all_counter = me_all_counter+1
			if( pass_epkpkmXe ){
			    me_pass_counter = me_pass_counter + 1
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_only_missing_energy')
			}
			else{
			    me_fail_counter = me_fail_counter + 1
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'fail_only_missing_energy')
			}
			
			
			// investigate background on kpkm mass spectrum
			// cut on the exclusivity mass values from 1 to 4 sigma while keeping others at 4 sigma (nominal cuts used above)
			// should have 16 plots from this, do this for special histograms 
			/*
			for( ii in 1..4 ){
 			    def pass_epkpkmXe_sig_lvl = epkpkmXe_cut_range[ii][0] < phykin.epkpkmX.e() && epkpkmXe_cut_range[ii][1] > phykin.epkpkmX.e()
 			    def pass_ekpkmX_sig_lvl = ekpkmX_cut_range[ii][0] < phykin.ekpkmX.mass2() && ekpkmX_cut_range[ii][1] > phykin.ekpkmX.mass2()
 			    def pass_epkpX_sig_lvl = epkpX_cut_range[ii][0] < phykin.epkpX.mass2() && epkpX_cut_range[ii][1] > phykin.epkpX.mass2()
 			    def pass_epkmX_sig_lvl = epkmX_cut_range[ii][0] < phykin.epkmX.mass2() && epkmX_cut_range[ii][1] > phykin.epkmX.mass2()

			    //lvls on missing energy
			    if(pass_epkpkmXe_sig_lvl && pass_epkpX && pass_epkmX && pass_ekpkmX ){
				fillEventSpecial(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_epkpkmxe_lvl'+ii)
			    }

			    //lvlvs on missing Km
			    if(pass_epkpkmXe && pass_epkpX_sig_lvl && pass_epkmX && pass_ekpkmX ){
				fillEventSpecial(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_epkpX_lvl'+ii)
			    }

			    // lvls on missing Kp
			    if(pass_epkpkmXe && pass_epkpX && pass_epkmX_sig_lvl && pass_ekpkmX ){
				fillEventSpecial(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_epkmX_lvl'+ii)
			    }

			    // lvls on missing proton
			    if(pass_ekpkmX_sig_lvl && pass_epkpkmXe && pass_epkpX && pass_epkmX ){
				fillEventSpecial(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_ekpkmX_lvl'+ii)
			    }
			    
			    // all have the same sigma lvl cut from 1 to 4 to look at the total effect
			    if(pass_epkpkmXe_sig_lvl && pass_epkpX_sig_lvl && pass_epkmX_sig_lvl && pass_ekpkmX_sig_lvl ){
				fillEventSpecial(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, 'pass_all_lvl'+ii)
			    }
			}
			*/
			    

			//restrict by cutting on all variables above
			if( pass_epkpkmXe && pass_epkpX && pass_epkmX && pass_ekpkmX ){ 
			    rec_ele = ele
			    rec_pro = pro
			    rec_kp  = kp
			    rec_km  = km //maybe fill as map?
			    
			    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all')
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
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity,  phykin, 'pass_all_low_pro_p')
			    }
			    if( !pass_pro_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_high_pro_p')
			    }
			    if( pass_kp_p_min && pass_km_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_low_kaonPM_p')
			    }
			    if( !pass_kp_p_min || !pass_km_p_min ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_high_kaonPM_p')
			    }


			    if( pass_cpl_pro ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km,helicity, phykin, 'pass_all_pass_cpl_pro')
			    }
			    if( !pass_cpl_pro ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_cpl_pro')
			    }
			    if( pass_cpl_kp ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_cpl_kp')
			    }
			    if( !pass_cpl_kp ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_cpl_kp')
			    }
			    if( pass_cpl_km ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_cpl_km')
			    }
			    if( !pass_cpl_km ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_cpl_km')
			    }
			    if( pass_cpl_pro && pass_cpl_kp && pass_cpl_km ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_cpl_all')
			    }
			    if( !(pass_cpl_pro && pass_cpl_kp && pass_cpl_km) ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_cpl_all')
			    }

			    // analyze the impact of the delta theta cuts on the final results 
			    // the dist. goes out to about 8, so let's try cutting far out from the peak to see if it helps remove background we missed
			    if( pass_delta_theta_pro ){
			    	fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_delta_theta_pro')
			    }
			    if( !pass_delta_theta_pro){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_delta_theta_pro')
			    }
			    if( pass_delta_theta_kp ){
			    	fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_delta_theta_kp')
			    }
			    if( !pass_delta_theta_kp){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_delta_theta_kp')
			    }
			    if( pass_delta_theta_km ){
			    	fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_delta_theta_km')
			    }
			    if( !pass_delta_theta_km){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_delta_theta_km')
			    }
			    if( pass_delta_theta_pro && pass_delta_theta_kp && pass_delta_theta_km ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_delta_theta_all')
			    }
			    if( !(pass_delta_theta_pro &&pass_delta_theta_kp && pass_delta_theta_km) ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_delta_theta_all')
			    }
				
			    
			    // limit theta to only defined FD coverage less than 35 degrees
			    if( pass_pro_theta_max && pass_kp_theta_max && pass_km_theta_max ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_pass_thetaFD')
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_fail_thetaFD')
			    }
			    
	 		    if( (kp+km).mass() < excl_cuts.kpkm_mass[1] && (kp+km).mass() > excl_cuts.kpkm_mass[0] ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 'pass_all_cut_phi_mass')
				fillProtonID(event, histos, ele, pro, proton.index, 'pass_all_cut_phi_mass')
				fillKaonPlusID(event, histos, ele, kp, kaonP.index, 'pass_all_cut_phi_mass')
				fillKaonMinusID(event, histos, ele, km, kaonM.index, 'pass_all_cut_phi_mass')
			    }

			    			   
			    // now add additional cut on the pid chi2 from EB for hadrons to see if it removes background near the signal			     
 			    for( ii in 1..6){
				def pass_pro_chi2 = Math.abs(event.chi2pid[ proton.index ]) < ii 
				def pass_kp_chi2 = Math.abs(event.chi2pid[ kaonP.index ]) < ii
				def pass_km_chi2 = Math.abs(event.chi2pid[ kaonM.index ]) < ii
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_raw_${ii}sigcut")
				if( pass_pro_chi2 && pass_kp_chi2 && pass_km_chi2 ){
				    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_pass_${ii}sigcut")
				}
				else{
				    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_fail_${ii}sigcut")
				}
			    }

			    //additional cuts that could be used - delta vz of kp-km

			    
			    def pass_dvz_kpkm = (event.vz[kaonP.index] - event.vz[kaonM.index]) < excl_cuts.delta_vz_kpkm[1] && (event.vz[kaonP.index] - event.vz[kaonM.index]) > excl_cuts.delta_vz_kpkm[0]
			    def pass_dvz_ep = Math.abs(event.vz[electron.index] - event.vz[proton.index]) < 6
			    // pass just dvz kpkm check is phi mass signal improves
			    if( pass_dvz_kpkm ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_pass_dvz_kpkm")
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_fail_dvz_kpkm")
			    }
			    // pass just dvz electron proton
			    if( pass_dvz_ep ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_pass_dvz_ep")
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_fail_dvz_ep")
			    }
			    //passs both the delta vz cuts
			    if( pass_dvz_kpkm && pass_dvz_ep ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_pass_dvz_both")
			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_fail_dvz_both")
			    }


			    // now include specific additional cuts to select final events.
			    def pass_pro_chi2 = Math.abs(event.chi2pid[ proton.index ]) < excl_cuts.pro_chi2[1] 
			    def pass_kp_chi2 = Math.abs(event.chi2pid[ kaonP.index ]) < excl_cuts.kp_chi2[1]
			    def pass_km_chi2 = Math.abs(event.chi2pid[ kaonM.index ]) < excl_cuts.km_chi2[1]			    
			    if( pass_dvz_kpkm && pass_pro_chi2 && pass_kp_chi2 && pass_km_chi2 &&
				pass_cpl_pro && pass_cpl_kp && pass_cpl_km ){
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_pass_additional")

				// assign helicity based on whether we are using simulation or not
				//def helicity = event.helicity
				if( event.mc_npart > 0 ){
				    helicity = gen_helicity_smear
				}				
				if( t1 != null && helicity != 0 && pass_not_resonant1 && pass_not_resonant2 ) {
				    writeOutEvent(t1, [ (kp+km).mass(), helicity, phykin.trentophi, phykin.q2, phykin.xb, phykin.t, phykin.w, event.run_number ] )
				}
				// now to have a look at the result of additional cuts + cut on phi mass peak 
				if( (kp+km).mass() < excl_cuts.kpkm_mass[1] && (kp+km).mass() > excl_cuts.kpkm_mass[0] ){
				    reconstructed_event_present=true
				    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, 
							"pass_all_pass_additional_pass_phi_mass")				     
				    fillDataResolutions( event, histos, ele, pro, kp, km, phykin, "pass_all_pass_additional_pass_phi_mass" )

 				    // this block is specific to smearing the simulated results to better match data
 				    if( event.mc_npart > 0 ){
					//def gen_ele_smear = LorentzVector.withPID(11, event.mc_px[0], event.mc_py[0], event.mc_pz[0])
					///def gen_pro_smear = LorentzVector.withPID(2212, event.mc_px[1], event.mc_py[1], event.mc_pz[1])
					///def gen_kp_smear = LorentzVector.withPID(321, event.mc_px[2], event.mc_py[2], event.mc_pz[2])
					///def gen_km_smear = LorentzVector.withPID(-321, event.mc_px[3], event.mc_py[3], event.mc_pz[3])
					//gen_ele_smear = smearParticleVector(gen_ele_smear, r, 11)//, 2.0, 2.5, 3.5)
					//gen_pro_smear = smearParticleVector(gen_pro_smear, r, 2212)//, 2.0, 2.5, 3.5)
					//gen_kp_smear = smearParticleVector(gen_kp_smear, r, 321)//, 2.0, 2.5, 3.5)
					//gen_km_smear = smearParticleVector(gen_km_smear, r, -321)//, 2.0, 2.5, 3.5)
					gen_phykin_with_smear =  getPhysicsKin(beam, target, gen_ele_smear, gen_pro_smear, gen_kp_smear, gen_km_smear)					
					fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, gen_ele_smear, gen_pro_smear, gen_kp_smear, gen_km_smear, helicity, 
							    gen_phykin_with_smear, "pass_all_pass_additional_pass_phi_mass_with_smear")
				    }

				    // cut on the theta proton angle in intervals of 1 degree to check the effect on the P_ele and -t				   
				    /*for( ii in 1..7 ){ 
					if( Math.toDegrees(pro.theta()) < 35-ii ){
					    fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_additional_pass_phi_mass_thetacut${ii}")
					}
				    }
				    if( pro.p() > 0.8 ){
					fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, "pass_all_pass_additional_pass_phi_mass_lowproPcut")					
				    }*/

				}

			    }
			    else{
				fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, helicity, phykin, "pass_all_fail_additional")
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
			    def sim_phykin = getPhysicsKin(beam, target, gen_ele, gen_pro, gen_kp, gen_km)
			    if( reconstructed_event_present ){

				fillSimulationResolutions(histos, event, rec_ele, rec_pro, rec_kp, rec_km, gen_ele, gen_pro, gen_kp, gen_km, phykin, sim_phykin, 'sim_resolutions' )
				fillSimulationHistograms(histos, beam, target, rec_ele, rec_pro, rec_kp, rec_km, gen_helicity_smear, 'sim_rec')
				// r is the instance of the random number class - didnt global variable it 
				def smear_ele = gen_ele_smear //smearParticleVector(smear_ele,r, 11 )
				def smear_pro = gen_pro_smear //smearParticleVector(smear_pro,r, 2212 )
				def smear_kp = gen_km_smear //smearParticleVector(smear_kp,r, 321 )
				def smear_km = gen_kp_smear //smearParticleVector(smear_km,r, -321 )

				def sim_phykin_with_smear = getPhysicsKin(beam, target, smear_ele, smear_pro, smear_kp, smear_km)
				fillSimulationResolutions(histos, event, rec_ele, rec_pro, rec_kp, rec_km, smear_ele, smear_pro, smear_kp, smear_km, phykin, sim_phykin_with_smear, 'sim_resolutions_with_smear' )
				
			    }
			    fillSimulationHistograms(histos, beam, target, gen_ele, gen_pro, gen_kp, gen_km, gen_helicity,'sim_gen')
			    //fillEventHistograms(event, histos, electron, proton, kaonP, kaonM, ele, pro, kp, km, phykin, '')
 			    
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

def root_out = new ROOTFile("monitor_phi_${outname}_pidtype1_cutvar1_rmvres12-${field_type}_16core_incl.root")
histos.each{root_out.writeDataSet(it.value)}
root_out.close()

//println(' Writing out txt file if needed')
if( f_out != null ){
    t1.write()
    f_out.close()
}


// cut var0
  // file name is monitor_phi_X_pidtype-inb.root

// cut var1 
   /// removed the cut on pr, kp, km polar angle < 35
   /// added loose coplanarity cut on each particle instead



// print out pass rate information
println('00000 ----- PASS RATE INFO FOR EXCLUSIVITY CUTS ON MASSES ------ 00000 ')
def pass_only_missing_energy = histos['q2_pass_only_missing_energy'].getEntries()
def fail_only_missing_energy = histos['q2_fail_only_missing_energy'].getEntries()
def all_missing_energy = histos['q2_combssize1'].getEntries()
println('total number of events with only 1 el, 1 pr, 1 kp, 1 km ' + me_all_counter)
println('--> total passing missing energy cut '  + me_pass_counter)
println('--> total failing missing energy cut '  + me_fail_counter)
println('--> sum of pass fail ' + (me_pass_counter +  me_fail_counter))


println(' ----- using get entries --- ' )
println('hist: total number of events with only 1 el, 1 pr, 1 kp, 1 km ' + all_missing_energy)
println('--> hist: total passing missing energy cut '  + pass_only_missing_energy)
println('--> hist: total failing missing energy cut '  + fail_only_missing_energy)

println(' ----- using get integral --- ' )
println('hist: total number of events with only 1 el, 1 pr, 1 kp, 1 km ' + histos['q2_combssize1'].getIntegral())
println('--> hist: total passing missing energy cut '  + histos['q2_pass_only_missing_energy'].getIntegral() )
println('--> hist: total failing missing energy cut '  + histos['q2_fail_only_missing_energy'].getIntegral() )

