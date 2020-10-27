package mon
import org.jlab.io.hipo.HipoDataSource
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.util.LinkedHashMap
import pid.electron.ElectronFromEvent
import pid.proton.ProtonPID
import pid.hadron.HadronID
import event.Event
import event.EventConverter
import utils.KinTool


class HadronMonV3{

    def beam = 10.604
    def rfPeriod = 4.008
    def rf_large_integer = 1000


    def field_config = ''
    def el_cut_strictness = null

    def electron = new ElectronFromEvent()
    
    def proton = new HadronID(hadron_id:2212)
    def kaonP = new HadronID(hadron_id:321)
    def pip =  new HadronID(hadron_id:211)
    def kaonM = new HadronID(hadron_id:-321)
    def pim =  new HadronID(hadron_id:-211)
    
    def setCutParameters(){
	println(' [HadronMonV3::setCutParameters] -> field configuration is ' + field_config)
	println(' [HadronMonV3::setCutParameters] -> setting cut values for field ' + field_config)
	electron.setElectronCutParameters(field_config)
    }
    
    def setCutStrictness(){
	println(' [HadronMonV3::setCutStrictness] -> cut strictness level is ' + el_cut_strictness )
	electron.setElectronCutStrictness(el_cut_strictness)
    }

    def myElectronCutStrategies = [
	electron.passElectronStatus,
	electron.passElectronChargeCut,
	electron.passElectronEBPIDCut,
	//electron.passElectronNpheCut,
	electron.passElectronVertexCut
	//electron.passElectronPCALFiducialCut,
	//electron.passElectronEIEOCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3
    ]

    def getLocalSector(x,y,z){
	def phi = Math.toDegrees( Math.atan2( y/ (Math.sqrt( x**2 + y**2 + z**2 )), x/Math.sqrt( x**2 + y**2 + z**2 ) ) )
	
	if(phi < 30 && phi >= -30){        return 1}
	else if(phi < 90 && phi >= 30){    return 2}
	else if(phi < 150 && phi >= 90){   return 3}
	else if(phi >= 150 || phi < -150){ return 4}
	else if(phi < -90 && phi >= -150){ return 5}
	else if(phi < -30 && phi >= -90){  return 6}
	
	return 0    
    }

    def getLocalPhi( local_sector, phi_DCr_raw ){
	def phi_DCr = 5000
	if ( local_sector == 1 ) phi_DCr = phi_DCr_raw
	if ( local_sector == 2 ) phi_DCr = phi_DCr_raw - 60
	if ( local_sector == 3 ) phi_DCr = phi_DCr_raw - 120
 	if ( local_sector == 4 && phi_DCr_raw > 0 ) phi_DCr = phi_DCr_raw - 180
	if ( local_sector == 4 && phi_DCr_raw < 0 ) phi_DCr = phi_DCr_raw + 180
	if ( local_sector == 5 ) phi_DCr = phi_DCr_raw + 120
	if ( local_sector == 6 ) phi_DCr = phi_DCr_raw + 60			
	return phi_DCr
    }
	
    def myProtonCuts = [
	proton.passHadronEBPIDCut,
	proton.passDCFiducialCutChi2Region1,
	proton.passDCFiducialCutChi2Region2,
	proton.passDCFiducialCutChi2Region3
    ]

    def myKaonPCuts = [
	kaonP.passHadronEBPIDCut,
	kaonP.passDCFiducialCutChi2Region1,
	kaonP.passDCFiducialCutChi2Region2,
	kaonP.passDCFiducialCutChi2Region3
    ]

    def myPIPCuts = [
	pip.passHadronEBPIDCut,
	pip.passDCFiducialCutChi2Region1,
	pip.passDCFiducialCutChi2Region2,
	pip.passDCFiducialCutChi2Region3
    ]

    def myKaonMCuts = [
	kaonM.passHadronEBPIDCut,
	kaonM.passDCFiducialCutChi2Region1,
	kaonM.passDCFiducialCutChi2Region2,
	kaonM.passDCFiducialCutChi2Region3
    ]

    def myPIMCuts = [
	pim.passHadronEBPIDCut,
	pim.passDCFiducialCutChi2Region1,
	pim.passDCFiducialCutChi2Region2,
	pim.passDCFiducialCutChi2Region3
    ]

    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	vzphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -20, 50 ) },
	vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
	hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
	ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.3) },
	uvwsf : { title -> new H2F("$title","$title", 500, 0.0, 500.0, 500, 0.0, 0.3) },
	tchi2 : { title -> new H1F("$title","$title", 500, 0.0, 500 ) }, 
	tndf : { title -> new H1F("$title","$title", 50, 0.0, 50 ) },
	nphe : { title -> new H1F("$title","$title", 40, 0.0, 40.0 ) },
	eieo : { title -> new H2F("$title","$title", 1000, 0.0, 1.0, 1000, 0.0, 1.0 ) },
	vz   : { title -> new H1F("$title","$title", 400, -40, 40) },
	vxy : { title -> new H2F("$title","$title", 400, -1, 1,400, -1, 1) },
	
	time : { title -> new H1F("$title","$title", 2000, -450, 450) },
	time2d : { title -> new H2F("$title","$title", 61, 0.0, 61.0, 2000, -450, 450) },
	path : { title -> new H1F("$title","$title", 2000, -400, 400) },//
	sectors : { title -> new H1F("$title","$title", 7, 0, 7) },
	helicity : { title -> new H1F("$title","$title", 6, -3 , 3) },
	w    : { title -> new H1F("$title","$title", 500, 0.5, beam+0.3 ) },
	q2w  : { title -> new H2F("$title","$title", 500, 0.5, beam+0.3, 500, 0, 2 * beam ) },  
	passrate : {title -> new H1F("$title","$title", 8, 0.0, 8.0 ) },
	wtheta : { title -> new H2F("$title","$title", 100, 0.0, 2.5, 100, 0.0, 20.0 ) },
 	deltakin : {title -> new H1F("$title","$title", 500, -2, 2  ) },
	deltakin2D : {title -> new H2F("$title","$title", 1000, 0.0, 500, 500, -2, 2  ) },
	deltakin2D_beam : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -2, 2  ) },
	deltat : {title -> new H2F("$title","$title", 1000, 0.0, 8, 1000, -2, 1  ) },
	deltab : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -1, 1  ) },
	betap : {title -> new H2F("$title","$title", 1000, 0.0, beam, 1000, 0.0, 1.4  ) },
	deltab_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -1, 1  ) },
	deltat_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -2, 2  ) },
	calib_time : { title -> new H2F("$title","$title", 1000, -100.0, 100, 1000, -100, 100  ) },
	betacomp : { title -> new H2F("$title","$title", 100, 0.80, 1.2, 100, 0.80, 1.2  ) },
		
	vz_vt : { title -> new H2F("$title","$title", 100, -70, 70, 100, -1.1, 101  ) },
	rf_vt : { title -> new H2F("$title","$title", 100, -70, 70, 100, -1.1, 101  ) },

	vz_p : { title -> new H2F("$title","$title", 600, 0, beam, 600, -50.0, 50.0  ) },
	vz_theta : { title -> new H2F("$title","$title", 600, 0, 120, 600, -50.0, 50.0  ) },
	vz_phi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -50.0, 50.0  ) },	
	vz_beta : { title -> new H2F("$title","$title", 600, -50, 50, 600, 0.0, 1.2  ) },
	
	delta_vz_p : { title -> new H2F("$title","$title", 600, 0, beam, 600, -5.0, 5.0  ) },
	delta_vz_theta : { title -> new H2F("$title","$title", 600, 0, 120, 600, -5.0, 5.0  ) },
	delta_vz_phi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -5.0, 5.0  ) },	
	delta_vz_beta : { title -> new H2F("$title","$title", 600, -5, 5, 600, 0.0, 1.2  ) },
	delta_vz : { title -> new H1F("$title","$title", 1000, -50.0, 50.0  ) },
	
	sum_time : { title -> new H1F("$title","$title", 5000, -500.0, 500.0  ) },
	sum_time2d : { title -> new H2F("$title","$title", 5000, -500.0, 500.0, 10, 0.0, 10.0 ) },
	mon_time : {title -> new H2F("$title","$title", 100,-rfPeriod/2,rfPeriod/2,100,-4,4) },	
	iron_cross : { title -> new H2F("$title","$title", 500, -20.0, 20.0, 500, -20.0, 20.0 ) },

    ]

    def tighter_kin_bounds = [
	theta_ele : [0, 45],
	p_ele     : [0, 10.5],
	prot_chi2 : [0, 100],
	prot_chi2_small : [-10, 10],
	vz        : [-20, 15],
	beta      : [0, 1.4],
	dvz       : [-25, 25],
	local_theta : [0, 60],
	local_phi : [-40, 40]
    ]

    def lim = tighter_kin_bounds
    
    def limited_h1 = { title, nbins, lims ->
	new H1F("$title", "$title", nbins, lims[0], lims[1])
    }
    
    def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
	new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
    }
    
    def histoBuilders1 = [
	pro_chi2 : { title -> limited_h1(title, 200, lim.prot_chi2_small) },
	dvz : {title -> limited_h1(title, 500, lim.dvz) }
    ]

    def histoBuilders2 = [
	p_chi2     : { title -> limited_h2(title, 250, 250, lim.p_ele, lim.prot_chi2_small)},
	theta_chi2 : { title -> limited_h2(title, 250, 250, lim.theta_ele, lim.prot_chi2_small)},
	beta_chi2  : { title -> limited_h2(title, 250, 250, lim.beta, lim.prot_chi2_small)},
	vz_chi2    : { title -> limited_h2(title, 250, 250, lim.vz, lim.prot_chi2_small)},
	pro_chi2_pro_theta : { title -> limited_h2(title, 500, 500, lim.theta_pro, lim.prot_chi2_small)},
	pro_chi2_pro_theta_ftof : { title -> limited_h2(title, 500, 500, lim.theta_pro_ftof, lim.prot_chi2_small)},
	pro_chi2_theta_sum : { title -> limited_h2(title, 500, 500, lim.vz, lim.prot_chi2_small)},
	p_beta      : {title -> limited_h2(title, 1000, 1000, lim.p_ele, lim.beta ) },
	local_hit : { title -> limited_h2(title, 500, 500, lim.local_phi, lim.local_theta) },
	local_hit_rotated : { title -> limited_h2(title, 500, 500, lim.local_theta, lim.local_phi) },
	deltat : {title -> limited_h2(title, 700, 0.0, beam, 700, -2.0 , 2.0  ) }

    ]

    def fillHadronInfo( event, histos, el, had_lv, index, comment){
	//pid chi2 info
	histos.computeIfAbsent('chi2_'+comment, histoBuilders1.pro_chi2).fill(event.chi2pid[index]) 
	histos.computeIfAbsent('p_chi2_'+comment, histoBuilders2.p_chi2).fill(had_lv.mag(), event.chi2pid[index]) 
	histos.computeIfAbsent('theta_chi2_'+comment, histoBuilders2.theta_chi2).fill(Math.toDegrees(had_lv.theta()), event.chi2pid[index]) 
	histos.computeIfAbsent('beta_chi2_'+comment, histoBuilders2.beta_chi2).fill(event.beta[index], event.chi2pid[index]) 
	histos.computeIfAbsent('vz_chi2_'+comment, histoBuilders2.vz_chi2).fill(event.vz[index], event.chi2pid[index]) 
	histos.computeIfAbsent('p_beta_'+comment, histoBuilders2.p_beta).fill( had_lv.mag(), event.beta[index])
	
	//p, theta, vz, phi info
	histos.computeIfAbsent("vz_p_"+comment,histoBuilders.vz_p).fill(had_lv.mag(), had_lv.vz)
	histos.computeIfAbsent("vz_theta_"+comment,histoBuilders.vz_theta).fill(Math.toDegrees(had_lv.theta()), had_lv.vz)
	histos.computeIfAbsent("vz_phi_"+comment,histoBuilders.vz_phi).fill(Math.toDegrees(had_lv.phi()), had_lv.vz)		    
	
	def delta_vz = el.vz - had_lv.vz
	histos.computeIfAbsent("delta_vz_p_"+comment,histoBuilders.delta_vz_p).fill(had_lv.mag(), delta_vz)
	histos.computeIfAbsent("delta_vz_theta_"+comment,histoBuilders.delta_vz_theta).fill(Math.toDegrees(had_lv.theta()), delta_vz)
	histos.computeIfAbsent("delta_vz_phi_"+comment,histoBuilders.delta_vz_phi).fill(Math.toDegrees(had_lv.phi()), delta_vz)    
	histos.computeIfAbsent("delta_vz_"+comment,histoBuilders1.dvz).fill(delta_vz)

	histos.computeIfAbsent("phi_theta_"+comment, histoBuilders.phitheta).fill( Math.toDegrees(had_lv.phi()), Math.toDegrees(had_lv.theta()) )
	
    }
    
    //dc r1 = layer 6, r2 = layer 18, layer 3 = 36
    def fillDriftChamberInfo( event, histos, el, had_lv, index, dc_region, comment ){
	

	def local_sector = -1 
	
	//use dcr2 to determine local variables in phi sector
	//redundant for checking R2, but it is okay
	def region = dc_region
	def layers = null
	def hits = null
	def hit = null
	if( dc_region == 1 ){
	    layers  = 6
	    hits = event.dc1[index]
	    hit = hits.find{it.layer==layers}
	}
	else if(dc_region == 2 ){
	    layers = 18
	    hits = event.dc2[index]
	    hit = hits.find{it.layer==layers}
	}
	else if(dc_region == 3 ){
	    layers = 36
	    hits = event.dc3[index]
	    hit = hits.find{it.layer==layers}
	}
	    

	def hitsR2 = event.dc2[index]
	def hitR2 = hitsR2.find{it.layer == 18 }
	if( hitR2 ) {local_sector=getLocalSector(hitR2.x,hitR2.y,hitR2.z)}

	if( hit ){
 	    histos.computeIfAbsent("dc"+region+"hit_"+comment,histoBuilders.hitxy).fill(hit.x, hit.y) 		    		
	    //determine local coordinate values ( theta and phi)
	    def theta_DCr = Math.toDegrees(  Math.acos(hit.z / Math.sqrt( hit.x**2 + hit.y**2 + hit.z**2 ) ) )
	    def phi_DCr_raw = Math.toDegrees( Math.atan2( hit.y / Math.sqrt( hit.x**2 + hit.y**2 + hit.z**2) , hit.x/Math.sqrt( hit.x**2 + hit.y**2 + hit.z**2) ) )
	    def phi_DCr = getLocalPhi( local_sector, phi_DCr_raw )

	    // by definition require > 0 
	    if( local_sector > 0 ){
		histos.computeIfAbsent("local_hit_dcr"+region+"_s"+local_sector+"_"+comment, histoBuilders2.local_hit).fill ( phi_DCr, theta_DCr )
		histos.computeIfAbsent("local_hit_rotated_dcr"+region+"_s"+local_sector+"_"+comment, histoBuilders2.local_hit_rotated).fill ( theta_DCr, phi_DCr )
	    }
	}

    }
    
    def processEvent( event, hevent ){

	def electron_candidates = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()
	
	def ev_sum_tof_time = 0
	def ev_sum_tof_time_pion_region = 0
	def n_tracks = 0
	def n_tracks_pion_region = 0
 	
 	def hadrons = electron_candidates.findResults{el_candidate -> ! el_candidate.value.contains(false) ? el_candidate.key : null }.findResults{ indx ->
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    	    
	    ele.vz = event.vz[indx]
	    return [indx, ele]
	}.collectMany{ indx, ele -> 
	    //look at charged hadrons (exclude neutrals) for when there is a electron present using EB PID

	    (0..<event.npart).findAll{ event.charge[it] != 0 && ( it != indx ) && (Math.abs( event.vz[indx] - event.vz[it]) < 20) && (event.vz[it]>-7 && event.vz[it]<1) && (event.status[it] >= 2000 && event.status[it] < 4000) }.findResults{index-> // &&  }.findResults{index->		
		
		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	
		def pid = event.pid[index]
		def charge = event.charge[index]
		def charge_sign = charge > 0 ? "pos" : "neg"         
		def ctof = (event.status[index] >= 2000 && event.status[index] < 4000) ? "FTOF" : "CTOF"
		def p = lv.mag()
		def theta = Math.toDegrees(lv.theta())
		def phi = Math.toDegrees(lv.phi())
		def beta = event.beta[index]
		def vz = event.vz[index]
		lv.vz = vz


		def start_time_rf_corrected = event.vt[index] // start time with vertex correction
		def elvz = event.vz[indx]
		def vt = event.vt[index]
 		def run_number = event.run_number
		def delta_vz_abs = Math.abs(vz - elvz)
		def delta_vz = (vz - elvz)

  		def pr_id_stat = myProtonCuts.collect{ pr_test -> pr_test(event,index) }
  		def kp_id_stat = myKaonPCuts.collect{ kp_test -> kp_test(event,index) }
  		def pip_id_stat = myPIPCuts.collect{ pip_test -> pip_test(event,index) }

  		def km_id_stat = myKaonMCuts.collect{ km_test -> km_test(event,index) }
  		def pim_id_stat = myPIMCuts.collect{ pim_test -> pim_test(event,index) }
		
		fillHadronInfo( event, histos, ele, lv, index, charge_sign+"_"+ctof)

		if( event.dc1_status.contains(index) &&  event.dc2_status.contains(index) && ctof=="FTOF" ){		    		    
		    if( charge > 0 ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pos" )
		    if( charge < 0 ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_neg" )
		    if( pr_id_stat[0] ){
			if( pr_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pro_pass_dcr1" )
			if( !pr_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pro_fail_dcr1" )
		    }
		    if( kp_id_stat[0] ){
			if( kp_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_kp_pass_dcr1" )
			if( !kp_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_kp_fail_dcr1" )
		    }
		    if( km_id_stat[0] ){
			if( km_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_km_pass_dcr1" )
			if( !km_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_km_fail_dcr1" )
		    }
		    if( pip_id_stat[0] ){
			if( pip_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pip_pass_dcr1" )
			if( !pip_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pip_fail_dcr1" )

		    }
		    if( pim_id_stat[0] ){
			if( pim_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pim_pass_dcr1" )
			if( !pim_id_stat[1] ) fillDriftChamberInfo( event, histos, ele, lv, index, 1, ctof+"_pass_pim_fail_dcr1" )
		    }
		}

 		if( event.dc2_status.contains(index) && ctof=="FTOF" ){		    
		    if( charge > 0 ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pos" )
		    if( charge < 0 ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_neg" )
		    if( pr_id_stat[0] ){
			if( pr_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pro_pass_dcr2" )
			if( !pr_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pro_fail_dcr2" )
		    }
		    if( kp_id_stat[0] ){
			if( kp_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_kp_pass_dcr2" )
			if( !kp_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_kp_fail_dcr2" )
		    }
		    if( km_id_stat[0] ){
			if( km_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_km_pass_dcr2" )
			if( !km_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_km_fail_dcr2" )
		    }
		    if( pip_id_stat[0] ){
			if( pip_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pip_pass_dcr2" )
			if( !pip_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pip_fail_dcr2" )

		    }
		    if( pim_id_stat[0] ){
			if( pim_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pim_pass_dcr2" )
			if( !pim_id_stat[2] ) fillDriftChamberInfo( event, histos, ele, lv, index, 2, ctof+"_pass_pim_fail_dcr2" )
		    }
		}

		if( event.dc3_status.contains(index) && event.dc2_status.contains(index) && ctof=="FTOF"  ){		    
		    if( charge > 0 ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pos" )
		    if( charge < 0 ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_neg" )
		    if( pr_id_stat[0] ){
			if( pr_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pro_pass_dcr3" )
			if( !pr_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pro_fail_dcr3" )
		    }
		    if( kp_id_stat[0] ){
			if( kp_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_kp_pass_dcr3" )
			if( !kp_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_kp_fail_dcr3" )
		    }
		    if( km_id_stat[0] ){
			if( km_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_km_pass_dcr3" )
			if( !km_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_km_fail_dcr3" )
		    }
		    if( pip_id_stat[0] ){
			if( pip_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pip_pass_dcr3" )
			if( !pip_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pip_fail_dcr3" )

		    }
		    if( pim_id_stat[0] ){
			if( pim_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pim_pass_dcr3" )
			if( !pim_id_stat[3] ) fillDriftChamberInfo( event, histos, ele, lv, index, 3, ctof+"_pass_pim_fail_dcr3" )
		    }
		}
		

		if (event.dc2_status.contains(index) ){		    
   		    def dc_sector = event.dc_sector[index]		    
		    
		    fillHadronInfo( event, histos, ele, lv, index, charge_sign+"_dcr2_s"+dc_sector+"_"+ctof)		    
		    fillHadronInfo( event, histos, ele, lv, index, charge_sign+"_dcr2_"+ctof)		    
		    
		    if(pid == 2212){
			fillHadronInfo( event, histos, ele, lv, index, "pro_s"+dc_sector+"_"+ctof)
			fillHadronInfo( event, histos, ele, lv, index, "pro_"+ctof)
			if(!pr_id_stat.contains(false)){
			    fillHadronInfo( event, histos, ele, lv, index, "pro_s"+dc_sector+"_"+ctof+"_passall")
			    fillHadronInfo( event, histos, ele, lv, index, "pro_"+ctof+"_passall")
			}			    
		    }
		    else if( pid == 321 ){
			fillHadronInfo( event, histos, ele, lv, index, "kp_s"+dc_sector+"_"+ctof)
			fillHadronInfo( event, histos, ele, lv, index, "kp_"+ctof)
			if(!kp_id_stat.contains(false)){
			    fillHadronInfo( event, histos, ele, lv, index, "kp_s"+dc_sector+"_"+ctof+"_passall")
			    fillHadronInfo( event, histos, ele, lv, index, "kp_"+ctof+"_passall")
			}
		    }
		    else if( pid == 211 ){
			fillHadronInfo( event, histos, ele, lv, index, "pip_s"+dc_sector+"_"+ctof)
			fillHadronInfo( event, histos, ele, lv, index, "pip_"+ctof)
			if(!pip_id_stat.contains(false)){
			    fillHadronInfo( event, histos, ele, lv, index, "pip_s"+dc_sector+"_"+ctof+"_passall")
			    fillHadronInfo( event, histos, ele, lv, index, "pip_"+ctof+"_passall")

			}
		    }
		    else if( pid == -321 ){
			fillHadronInfo( event, histos, ele, lv, index, "km_s"+dc_sector+"_"+ctof)
			fillHadronInfo( event, histos, ele, lv, index, "km_"+ctof)
			if(!km_id_stat.contains(false)){
			    fillHadronInfo( event, histos, ele, lv, index, "km_s"+dc_sector+"_"+ctof+"_passall")
			    fillHadronInfo( event, histos, ele, lv, index, "km_"+ctof+"_passall")
			}

		    }
		    else if( pid == -211 ){
			fillHadronInfo( event, histos, ele, lv, index, "pim_s"+dc_sector+"_"+ctof)
			fillHadronInfo( event, histos, ele, lv, index, "pim_"+ctof)
			if(!pim_id_stat.contains(false)){
			    fillHadronInfo( event, histos, ele, lv, index, "pim_s"+dc_sector+"_"+ctof+"_passall")
			    fillHadronInfo( event, histos, ele, lv, index, "pim_"+ctof+"_passall")
			}
		    }	     
		}
	
		if( event.tof_status.contains(index) ){
		    event.tof[index].eachWithIndex{ hh, ii ->
			if( hh.layer == 2 ){ // 1b
			    def ftof_sector = hh.sector
 			    def ftof_t = hh.time
			    def ftof_p = hh.path
			    def ftof_paddle = hh.paddle
			    def ftof_layer = hh.layer
			    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    		   
			    def measured_t  = ftof_t - vt
			    def measured_t_with_starttime=  ftof_t - event.start_time

			    def calc_beta_pr = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(2212)**2)
			    def calc_time_pr = ftof_p/(calc_beta_pr*(29.9792458))
			    def delta_t_pr = measured_t - calc_time_pr
			    def delta_t_pr_with_evst = measured_t_with_starttime - calc_time_pr // convention

			    def calc_beta_kp = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(321)**2)
			    def calc_time_kp = ftof_p/(calc_beta_kp*(29.9792458))
			    def delta_t_kp = measured_t - calc_time_kp
			    def delta_t_kp_with_evst = measured_t_with_starttime - calc_time_kp // convention
			    
			    histos.computeIfAbsent("ftof_deltat_pr",histoBuilders.deltat).fill(p, delta_t_pr)
			    histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pr)

			    histos.computeIfAbsent("ftof_deltat_evst_pr",histoBuilders.deltat).fill(p, delta_t_pr_with_evst)
			    histos.computeIfAbsent("ftof_deltat_evst_pr_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pr_with_evst)
			    
			    histos.computeIfAbsent("ftof_deltat_kp",histoBuilders.deltat).fill(p, delta_t_kp)
			    histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_kp)

			    histos.computeIfAbsent("ftof_deltat_evst_kp",histoBuilders.deltat).fill(p, delta_t_kp_with_evst)
			    histos.computeIfAbsent("ftof_deltat_evst_kp_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_kp_with_evst)

			    if( event.pid[index] == 321 ){
				fillHadronInfo( event, histos, ele, lv, index, charge_sign+"_kp_ftof1b_$ftof_sector"+ctof )
			    }
			}
			
		    }

		}
			
	    }
	}     			   
    }
}
