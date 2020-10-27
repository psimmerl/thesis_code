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
import event.Event
import event.EventConverter
import utils.KinTool


class ElvsPionMon{

    def beam = 10.604
    def target_position = -3.0 //cm
    def speed_of_light = 29.9792458 // cm/ns
    def rfBucketLength = 4.0 //ns
    def lv_beam = LorentzVector.withPID(11,0,0,10.604)
    def target = LorentzVector.withPID(2212,0,0,0)
//el_p_passall_pidchi2
    def el_pid = 11
    def pim_pid = -211
    def el_cut_strictness=[]

    //def rf_offset = [5164:
    //5032: 
    
    //call electron cut constructor
    def field_config = ''
    def electron = new ElectronFromEvent();
    def part_names = ['el','pim']
    
    def setCutParameters(){
	println(' [epkpkm::setCutParameters] -> field configuration is ' + field_config)
	println(' [epkpkm::setCutParameters] -> setting cut values for field ' + field_config)
	electron.setElectronCutParameters(field_config)
    }

    def setCutStrictness(){
	println(' [epkpkm::setCutStrictness] -> cut strictness level is ' + el_cut_strictness )
	electron.setElectronCutStrictness(el_cut_strictness)
    }
    
    def everythingElsePassed(results, i) {
	for ( def j=0; j<results.size(); j++){	 
	    if (i != j && results[j] == false){		
		return false	
	    }	    
	}	
	return true	
    }
    
    def everythingElseFailed(results, i){	
	def reversed = results.collect{ !it }	
	return everythingElsePassed(reversed,i)	
    }

    def getVertexTime( time_at_detector, path_to_detector, calculated_beta ){
	def vertex_time = time_at_detector - path_to_detector/(calculated_beta * speed_of_light) 
	return vertex_time	  
    }
    
    def getStartTime( vertex_time, vertex_pos ){
	def vzCorr = (-target_position + vertex_pos)/speed_of_light
	//def deltatr = - vertex_time - vzCorr 
	//def rfCorr = 0//deltatr % rfBucketLength - rfBucketLength/2.0
	return vertex_time - vzCorr
    }
    
    def max_u_med   = [408, 408, 408, 408, 408, 408]
    def min_u_med   = [33, 26, 34, 22,   27, 25]
    def min_v_med   = [16.0, 10.5, 17.0,  14.25, 18.0, 11.0]
    def max_v_med   = [400, 400, 400, 400, 400, 400]
    def min_w_med   = [11.0, 17.5, 16.25, 7.5,  14.5, 9.25]
    def max_w_med   = [400, 400, 400, 400, 400, 400]

    def s_max_u=max_u_med
    def s_max_v=max_v_med
    def s_max_w=max_w_med
    def s_min_u=min_u_med
    def s_min_v=min_v_med
    def s_min_w=min_w_med
    
    //include PCAL fiducial cut in this study
    def passPCALFiducialCut( event , index ){
	if (event.pcal_status.contains(index)){	  
	    def pcal_sector = event.pcal_sector[index]-1
	    return (event.pcal_u[index].with{it < s_max_u[pcal_sector] && it > s_min_u[pcal_sector] } &&
		    event.pcal_v[index].with{it < s_max_v[pcal_sector] && it > s_min_v[pcal_sector] } &&
		    event.pcal_w[index].with{it < s_max_w[pcal_sector] && it > s_min_w[pcal_sector] })
	}
	return false
    }

    def passMntmRECGENCut( event, index ){
	def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	def mclv = new Vector3(event.mc_px[0],event.mc_py[0],event.mc_pz[0])		
	def res_p = lv.mag() - mclv.p()
	if( res_p < 0.25 && res_p > -0.25 ) return true
	return false
		
    }

    def passThetaRECGENCut( event, index ){
	def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	def mclv = new Vector3(event.mc_px[0],event.mc_py[0],event.mc_pz[0])		
    	def res_theta = Math.toDegrees(lv.theta() - mclv.theta())
	if( res_theta < 1.5 && res_theta > -1.5 ) return true
	return false
    }
    
    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	particle_counter : { title -> new H1F("$title","$title", 9.0, 1.0, 10.0 ) },
	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	hp : {title -> new H1F("$title","$title", 100, 0, 10.00) },
	htheta : {title -> new H1F("$title","$title", 900, 0, 65.0) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	vzphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, -20, 50 ) },
	vztheta : { title -> new H2F("$title","$title", 1000, -100, 100, 900, 0.0, 65.0 ) },
	hitxy : { title -> new H2F("$title","$title", 900, -450.0, 450.0, 900, -450.0, 450.0 ) },
	ecsf : { title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 0.5) },
	ecsf_chi : { title -> new H1F("$title","$title", 100, -5.0, 5.0) },	
	calsf : { title -> new H2F("$title","$title", 500, 0.0, 0.35, 500, 0.0, 0.35) },
	uvwsf : { title -> new H2F("$title","$title", 500, 0.0, 500.0, 500, 0.0, 0.08) },
	uvwhit : { title -> new H1F("$title","$title", 500, 0.0, 100.0) },
	tchi2 : { title -> new H1F("$title","$title", 500, 0.0, 500 ) }, 
	tndf : { title -> new H1F("$title","$title", 75, 0.0, 75 ) },
	nphe : { title -> new H1F("$title","$title", 75, 0.0, 75.0 ) },
	eieo : { title -> new H2F("$title","$title", 1000, 0.0, 2.0, 1000, 0.0, 2.0 ) },
	eieo_small : { title -> new H2F("$title","$title", 1000, 0.0, 0.4, 1000, 0.0, 0.4 ) },
	vz   : { title -> new H1F("$title","$title", 600, -20, 40) },
	vxy : { title -> new H2F("$title","$title", 1000, -1, 1, 1000, -1, 1) },

	time : { title -> new H1F("$title","$title", 2000, -450, 450) },
	time_small : { title -> new H1F("$title","$title",1000, -25, 25) },
	iron_cross : { title -> new H2F("$title","$title", 500, -20.0, 20.0, 500, -20.0, 20.0 ) },

	time2d_small : { title -> new H2F("$title","$title", 100, 0.0, 10.0, 500, -25, 25) },
	time2d_large : { title -> new H2F("$title","$title", 61, 0.0, 61.0, 2000, -450, 450) },
	betatime : { title -> new H2F("$title","$title", 200, 0.0, 1.20, 200, 50, 200) },
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
	hadron_multiplicity : {title -> new H1F("$title","$title", 10, 0.0, 10.0 ) },
	recpart_chi2pid : {title -> new H1F("$title","$title", 200, -5.0, 5.0 ) },
	recpart_chi2pid2D : {title -> new H2F("$title","$title", 200, 0.0,  beam+0.3, 200, -10.0, 10.0 ) },
	betap : { title -> new H2F("$title","$title", 200, 0.0, 8.0, 200, 0.0, 1.20 ) },
	second_moments : {  title -> new H1F("$title","$title", 700, -50.0, 300.0)},
	res_p : { title -> new H2F("$title","$title", 200, 0.0, 10.0, 200, -0.75, 0.75 ) },
	res_theta : { title -> new H2F("$title","$title", 200, 0.0, 40.0, 200, -2.0, 2.0 ) },
	res_phi : { title -> new H2F("$title","$title", 200, -180.0, 180.0, 200, -0.15, 0.15 ) }
	
    ]

    def processEvent( event ){
	
	///println("[ElectronMon::processEvent] - processing event now")
	def run_number = event.run_number 
	//&& passPCALFiducialCut(event,it)
	//def el_cand = (0..<event.npart).findAll{event.charge[it]<0 && Math.abs(event.vz[it]+3.0) < 5 && (Math.abs(event.status[it]) >= 2000 && Math.abs(event.status[it]) < 4000) && passPCALFiducialCut(event,it) && passMntmRECGENCut(event,it) && passThetaRECGENCut(event,it)  }.collect{ ii -> [ii, [(event.pid[ii] == el_pid )] ] }.collectEntries()	
	//def pim_cand = (0..<event.npart).findAll{event.charge[it]<0 && Math.abs(event.vz[it]+3.0) < 5 && (Math.abs(event.status[it]) >= 2000 && Math.abs(event.status[it]) < 4000) && passPCALFiducialCut(event,it) && passMntmRECGENCut(event,it) && passThetaRECGENCut(event,it)}.collect{ ii -> [ii, [(event.pid[ii] == pim_pid )] ] }.collectEntries()	
	//println("[ElectronMon::processEvent] - Found Good Electrons")
	//println(my_el_cuts)
	
	//for( ii in 0..<event.npart){
	 //   println(event.pid[ii])
	//}

 	def rec_mc_assigned_pid_ele = 0
 	def rec_mc_assigned_pid = 0
	def rec_mc_assigned_pid_pim = 0
	
	for( mm in 0..<event.mc_npart){	    
	    
	    // find reconstructed track that best matches the mc generated track
	    def mc_lv = new Vector3(event.mc_px[mm], event.mc_py[mm], event.mc_pz[mm])	    
	    def gen_p = mc_lv.mag()
	    def gen_theta = mc_lv.theta()
	    def gen_phi = mc_lv.phi()
	    def gen_vz = event.mc_vz[mm]
	    for( rr in 0..<event.npart){
		def lv = new Vector3(event.px[rr], event.py[rr], event.pz[rr])
		def rec_pid = event.pid[rr]

		def rec_p = lv.mag()
		def rec_theta = lv.theta()
		def rec_phi = lv.phi()
		def rec_vz = event.vz[rr]

		def delta_p = rec_p - gen_p
		def delta_theta = Math.toDegrees( rec_theta - gen_theta )
		def delta_phi = Math.toDegrees( rec_phi - gen_phi )
		def delta_vz = (rec_vz - gen_vz)

		// plot delta rec -gen
		def particle_name = null
		
		if( rec_pid == 11){ 
		    particle_name = "el"
		}
		else if(rec_pid == -211){
		    particle_name = "pim"
		}
		if( particle_name != null ){
		    histos.computeIfAbsent(particle_name+"_res_p",histoBuilders.res_p).fill(delta_p, gen_p)
		    histos.computeIfAbsent(particle_name+"_res_theta",histoBuilders.res_theta).fill(delta_theta, gen_theta)
		    histos.computeIfAbsent(particle_name+"_res_phi",histoBuilders.res_phi).fill(delta_phi, gen_phi)
		}


		def pass_vz = Math.abs(event.vz[rr]+3.0) < 5
		def pass_status = (Math.abs(event.status[rr]) >= 2000 && Math.abs(event.status[rr]) < 4000)
		def pass_status_trigger = event.status[rr]<0
		
		//identify only those tracks in target and in FD		
		if( Math.abs(delta_p) < 0.25 && Math.abs(delta_theta) < 1.5 && Math.abs(delta_phi) < 1.5 && Math.abs(delta_vz) < 2
		   && pass_vz && pass_status_trigger ){
		    rec_mc_assigned_pid = event.mc_pid[mm]
		    
		    // get calorimeter based information
		    def pcal_edep = 0
		    def eical_edep = 0
		    def eocal_edep = 0
		    
		    if( event.pcal_status.contains(rr) ){
			pcal_edep = event.pcal_energy[rr]
		    }
		    if( event.ecal_inner_status.contains(rr) ){
			eical_edep = event.ecal_inner_energy[rr]		    
		    }
		    if( event.ecal_outer_status.contains(rr) ){
			eocal_edep = event.ecal_outer_energy[rr]			
		    }

		    def cal_edep = pcal_edep + eical_edep + eocal_edep 
		    def p = rec_p
		    // this is reconstructed electron coming from generated electron
		    if( rec_mc_assigned_pid == 11 && event.pid[rr] == 11 ){
			//println(' wow gen electron as rec electron ')
			histos.computeIfAbsent("vz_elASel",histoBuilders.vz).fill(rec_vz)
 			histos.computeIfAbsent("ptheta_elASel",histoBuilders.ptheta).fill(rec_p, Math.toDegrees(rec_theta))
			histos.computeIfAbsent("p_elASel",histoBuilders.hp).fill(rec_p)
			histos.computeIfAbsent("phitheta_elASel",histoBuilders.phitheta).fill(Math.toDegrees(rec_phi), Math.toDegrees(rec_theta))
			histos.computeIfAbsent("calsf_elASel",histoBuilders.ecsf).fill(p, cal_edep/p)
		    }
		
		    // this is a reconstructed electron coming from a generated negative pion
		    if( rec_mc_assigned_pid == -211 && event.pid[rr] == 11 ) {
		   	//println(' wow gen negative pion as rec electron ')
			histos.computeIfAbsent("vz_pimASel",histoBuilders.vz).fill(rec_vz)
 			histos.computeIfAbsent("ptheta_pimASel",histoBuilders.ptheta).fill(rec_p, Math.toDegrees(rec_theta))
			histos.computeIfAbsent("p_pimASel",histoBuilders.hp).fill(rec_p)
			histos.computeIfAbsent("phitheta_pimASel",histoBuilders.phitheta).fill(Math.toDegrees(rec_phi), Math.toDegrees(rec_theta))
			histos.computeIfAbsent("calsf_pimASel",histoBuilders.ecsf).fill(p, cal_edep/p)
		    }

		}		   
	    }

	}

	

	/*
	def part_name_cc = 0
	def identified_particles = [el_cand, pim_cand]
	for( part in identified_particles ){
	    def p_name = part_names[part_name_cc]
	    part_name_cc+=1

	    part.each{ index, value ->
		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		def mclv = LorentzVector.withPID(11,event.mc_px[0],event.mc_py[0],event.mc_pz[0])		

		def p = lv.mag()
		def pid = event.pid[index]
		def vz = event.vz[index]
		def theta = Math.toDegrees(lv.theta())
		def phi = Math.toDegrees(lv.phi())
		def pidchi2 = event.chi2pid[index]
		def vx = event.vx[index]
		def vy = event.vy[index]
		def measured_el = LorentzVector.withPID(11,event.px[index],event.py[index],event.pz[index])
 		def eX = lv_beam+target-measured_el
		def q = lv_beam - measured_el
		def q2 = -q.mass2()
		def xb = q2/(2*target.mass()*q.e())
		def nu = q2/(2*target.mass()*xb)
		def y = nu/lv_beam.e()	    

		def delta_p = 0.100
		def p_min = 0
		def p_max = 0.100
		def p_range = -1
 		
 	        //define values here
		def nphe = null
		def cherenkov_sector = null
		def e_pcal = null
		def pcal_sect = null
		def pcal_u = null
		def pcal_v = null
		def pcal_w = null
		def pcal_t = null
		def pcal_l = null

		def e_eical = null
		def eical_sect = null
		def eical_u = null
		def eical_v = null
		def eical_w = null
		def eical_t = null
		def eical_l = null

		def e_eocal = null
		def eocal_sect = null
		def eocal_u = null
		def eocal_v = null
		def eocal_w = null
		def eocal_t = null
		def eocal_l = null
		
		def pcal_m2u = null
		def pcal_m2v = null
		def pcal_m2w = null
		def pcal_m2 = null

		def eical_m2u = null
		def eical_m2v = null
		def eical_m2w = null
		def eical_m2 = null

		def eocal_m2u = null
		def eocal_m2v = null
		def eocal_m2w = null
		def eocal_m2 = null

		def eidep = 0
		def eodep = 0
		def pcaldep = 0
		def cal_sector = -1
		def ei_sect = -1
		def eo_sect = -1 

		def dc_sector = -1
		def dc_chi2 = null 
		def dc_ndf = null

		if( event.dc1_status.contains(index)) {
		    dc_sector  = event.dc_sector[index]		    
		    dc_chi2 = event.dc_chi2[index]
		    dc_ndf = event.dc_ndf[index]
		}
		
		if (event.cherenkov_status.contains(index)) { 
		    nphe = event.nphe[index]	    
		    cherenkov_sector = event.cherenkov_sector[index]
		}	    

		if( event.pcal_status.contains(index) ){
		    e_pcal = event.pcal_energy[index]
		    pcal_sect = event.pcal_sector[index]
		    pcal_u = event.pcal_u[index]
		    pcal_v = event.pcal_v[index]
		    pcal_w = event.pcal_w[index]
		    pcal_t = event.pcal_time[index]
		    pcal_l = event.pcal_path[index]
		    pcal_m2u = event.pcal_m2u[index]
		    pcal_m2v = event.pcal_m2v[index]
		    pcal_m2w = event.pcal_m2w[index]
		    pcal_m2 = pcal_m2u+pcal_m2v+pcal_m2w
		}

		if( event.ecal_inner_status.contains(index) ){

		    eidep = event.ecal_inner_energy[index]
		    ei_sect = event.ecal_inner_sector[index]
	    	    cal_sector = event.ecal_inner_sector[index] - 1
		    e_eical = event.ecal_inner_energy[index]
		    eical_sect = ei_sect //redundant
		    eical_u = event.ecal_inner_u[index]
		    eical_v = event.ecal_inner_v[index]
		    eical_w = event.ecal_inner_w[index]
		    eical_t = event.ecal_inner_time[index]
		    eical_l = event.ecal_inner_path[index]
 		    eical_m2u = event.ecal_inner_m2u[index]
		    eical_m2v = event.ecal_inner_m2v[index]
		    eical_m2w = event.ecal_inner_m2w[index]		    
		    eical_m2 = eical_m2u + eical_m2v + eical_m2w
		}
		if( event.ecal_outer_status.contains(index) ){
		    eodep = event.ecal_outer_energy[index]
 		    eo_sect = event.ecal_outer_sector[index]
		    cal_sector = event.ecal_outer_sector[index]-1

		    e_eocal = event.ecal_outer_energy[index]
		    eocal_sect = eo_sect //redundant
		    eocal_u = event.ecal_outer_u[index]
		    eocal_v = event.ecal_outer_v[index]
		    eocal_w = event.ecal_outer_w[index]
		    eocal_t = event.ecal_outer_time[index]
		    eocal_l = event.ecal_outer_path[index]
 		    eocal_m2u = event.ecal_outer_m2u[index]
		    eocal_m2v = event.ecal_outer_m2v[index]
		    eocal_m2w = event.ecal_outer_m2w[index]		    
		    eocal_m2 = eocal_m2u + eocal_m2v + eocal_m2w
		}

		if( event.pcal_status.contains(index) ){
		    pcaldep = event.pcal_energy[index]
		    cal_sector = event.pcal_sector[index]-1
		}
		
		def edep = eidep + eodep + pcaldep
		def sf_chi = KinTool.getSamplingFractionFromMomentum(cal_sector, p , edep/p)

		value.eachWithIndex{ cut_res, cut_ind -> cut_res ? histos.computeIfAbsent(p_name+"_"+"passResults",histoBuilders.passrate).fill(cut_ind) : null }
		
		value.eachWithIndex{ cut_result, cut_index ->		
		    
		    if(cut_result && cut_index==0){
			if(nphe != null ) {histos.computeIfAbsent(p_name+"_"+"nphe_cut$cut_index",histoBuilders.nphe).fill(nphe)}
			if(nphe != null ) histos.computeIfAbsent(p_name+"_"+"nphe_s"+cherenkov_sector+"_cut$cut_index",histoBuilders.nphe).fill(nphe)

	 		if( e_pcal != null ) histos.computeIfAbsent(p_name+"_"+"pcalsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p) 

 			if( e_eical != null ) histos.computeIfAbsent(p_name+"_"+"eicalsf_cut$cut_index"+"_s$eical_sect",histoBuilders.ecsf).fill(p, e_eical/p) 

			if( eocal_sect >= 0 ){
			    histos.computeIfAbsent(p_name+"_"+"eocalsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.ecsf).fill(p, e_eocal/p) 
			}

 			if(cal_sector >= 0 ){
			    histos.computeIfAbsent(p_name+"_"+"calsfchi_cut$cut_index",histoBuilders.ecsf_chi).fill(pidchi2) 	
			    histos.computeIfAbsent(p_name+"_"+"calsfchi_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(pidchi2) 	
		    	    
			    histos.computeIfAbsent(p_name+"_"+"calsf_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
			    histos.computeIfAbsent(p_name+"_"+"calsf_cut$cut_index",histoBuilders.ecsf).fill(p, edep/p) 	
			    histos.computeIfAbsent(p_name+"_"+"pcalvsecal_cut$cut_index"+"_s$cal_sector",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 		    histos.computeIfAbsent(p_name+"_"+"pcalvsecal_cut$cut_index",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
			    
			    
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    

			    histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_eivspcal_cut$cut_index",histoBuilders.calsf).fill(eidep/p, pcaldep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    
			    
			    for( int sig_region = 1; sig_region <= 5; sig_region++ ){			
				if( Math.abs(pidchi2) <= sig_region ){
				    histos.computeIfAbsent(p_name+"_"+"calsf_cut$cut_index"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p) 	
				    histos.computeIfAbsent(p_name+"_"+"calsf_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p) 	
				    histos.computeIfAbsent(p_name+"_"+"pcalvsecal_cut$cut_index"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 			    histos.computeIfAbsent(p_name+"_"+"pcalvsecal_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
				    
				    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
 				    histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(eidep/p, eodep/p)
				    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				    histos.computeIfAbsent(p_name+"_"+"calsf_eivspcal_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(eidep/p, pcaldep/p)
 				    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
				    
				}
			    }
 			    
 			    if( p > 4.0 ){
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
				histos.computeIfAbsent(p_name+"_"+"calsfchi_bigp_cut$cut_index",histoBuilders.ecsf_chi).fill(pidchi2) 	
				histos.computeIfAbsent(p_name+"_"+"calsfchi_bigp_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(pidchi2) 	
				
 				histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_bigp_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    }
			    else if ( p < 4.0 ){
				histos.computeIfAbsent(p_name+"_"+"calsfchi_smallp_cut$cut_index",histoBuilders.ecsf_chi).fill(pidchi2) 	
				histos.computeIfAbsent(p_name+"_"+"calsfchi_smallp_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(pidchi2) 	
			
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
				//if( eodep > 0 && eidep > 0 )
				histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_lowp_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
				//if( pcaldep > 0 && eidep > 0 )
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				//if( pcaldep > 0 && eodep > 0 )
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    }			    
			}

			if( ei_sect == eo_sect && ei_sect >= 0 ){			
			    histos.computeIfAbsent(p_name+"_"+"eieo_cut$cut_index",histoBuilders.eieo).fill(e_eical, e_eocal)			    
			    histos.computeIfAbsent(p_name+"_"+"eieo_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)			    
			    histos.computeIfAbsent(p_name+"_"+"eieo_small_cut$cut_index",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
			    histos.computeIfAbsent(p_name+"_"+"eieo_small_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
			    //if( e_eocal > 0 && e_eical > 0 ) histos.computeIfAbsent(p_name+"_"+"ecalsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
			    histos.computeIfAbsent(p_name+"_"+"ecalsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
			}

			histos.computeIfAbsent(p_name+"_"+"vz_cut$cut_index",histoBuilders.vz).fill(vz)
			histos.computeIfAbsent(p_name+"_"+"ptheta_cut$cut_index",histoBuilders.ptheta).fill(p, theta)
                        histos.computeIfAbsent(p_name+"_"+"p_cut$cut_index",histoBuilders.hp).fill(p)
			histos.computeIfAbsent(p_name+"_"+"theta_cut$cut_index",histoBuilders.htheta).fill(theta)
			histos.computeIfAbsent(p_name+"_"+"phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
			histos.computeIfAbsent(p_name+"_"+"vztheta_cut$cut_index",histoBuilders.vztheta).fill(vz, theta)
			histos.computeIfAbsent(p_name+"_"+"vxy_cut$cut_index",histoBuilders.vxy).fill(vx,vy)
			histos.computeIfAbsent(p_name+"_"+"vzphi_cut$cut_index",histoBuilders.vzphi).fill(phi,vz)
		    }
		}
		
		def res_p = lv.mag() - mclv.p()
		def res_theta = lv.theta() - mclv.theta()
		def res_phi = lv.phi() - mclv.phi()
		if(!value.contains(false)){		    
  		    histos.computeIfAbsent(p_name+"_"+"res_p",histoBuilders.res_p).fill(mclv.p(), res_p)
		    histos.computeIfAbsent(p_name+"_"+"res_theta",histoBuilders.res_theta).fill(mclv.theta()*180.0/3.1415, res_theta*180.0/3.1415)		
		    histos.computeIfAbsent(p_name+"_"+"res_phi",histoBuilders.res_phi).fill(mclv.phi()*180.0/3.1415, res_phi*180.0/3.1415)
		}
		if( value.contains(false)){
		    histos.computeIfAbsent(p_name+"_"+"res_p_fail_at_least1",histoBuilders.res_p).fill(mclv.p(), res_p)
		    histos.computeIfAbsent(p_name+"_"+"res_theta_fail_at_least1",histoBuilders.res_theta).fill(mclv.theta()*180.0/3.1415, res_theta*180.0/3.1415)		
		    histos.computeIfAbsent(p_name+"_"+"res_phi_fail_at_least1",histoBuilders.res_phi).fill(mclv.phi()*180.0/3.1415, res_phi*180.0/3.1415)
		}		
	    	
		//plot nphe for tracks passing all cuts
		if( !value.contains(false) ){
		    if( nphe != null ) histos.computeIfAbsent(p_name+"_"+"nphe_passall",histoBuilders.nphe).fill(nphe)
		    if( nphe != null ) histos.computeIfAbsent(p_name+"_"+"nphe_passall_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  
		}		    
		if( value.contains(false) ){
		    if(nphe != null ) histos.computeIfAbsent(p_name+"_"+"nphe_fail_at_least1",histoBuilders.nphe).fill(nphe)
		    if(nphe != null ) histos.computeIfAbsent(p_name+"_"+"nphe_fail_at_least1_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  
		}		    
		

		if(!value.contains(false)){
		    if( e_pcal != null )  histos.computeIfAbsent(p_name+"_"+"pcalsf_passall"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p)
		}		
 		else if(value.contains(false)){
		    if( e_pcal != null ) histos.computeIfAbsent(p_name+"_"+"pcalsf_fail_at_least1"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p)
		}		
		
		if(!value.contains(false)){
		    if( e_eical != null ) histos.computeIfAbsent(p_name+"_"+"eicalsf_passall_s$eical_sect",histoBuilders.ecsf).fill(p, e_eical/p) 
		}
		if(value.contains(false)){
		    if( e_eical != null ) histos.computeIfAbsent(p_name+"_"+"eicalsf_fail_at_least1_s$eical_sect",histoBuilders.ecsf).fill(p, e_eical/p)
		}

		if( eocal_sect >= 0 ){
 		    if(!value.contains(false)){
			histos.computeIfAbsent(p_name+"_"+"eocalsf_passall_s$eocal_sect",histoBuilders.ecsf).fill(p, e_eocal/p) 
		    }
 		    if(value.contains(false)){
			histos.computeIfAbsent(p_name+"_"+"eocalsf_fail_at_least1_s$eocal_sect",histoBuilders.ecsf).fill(p, e_eocal/p) 
		    }
		}

 		if(true){//cal_sector >= 0 ){
		    if(!value.contains(false)){		    
			histos.computeIfAbsent(p_name+"_"+"calsf_passall"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
			histos.computeIfAbsent(p_name+"_"+"calsf_passall",histoBuilders.ecsf).fill(p, edep/p) 	
			//histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall"+"_s$cal_sector",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
			//histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall"+"_s$cal_sector",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 		histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall",histoBuilders.eieo).fill(pcaldep,eidep + eodep)
			histos.computeIfAbsent(p_name+"_"+"calsfchi_passall",histoBuilders.ecsf_chi).fill(pidchi2) 	
			histos.computeIfAbsent(p_name+"_"+"calsfchi_passall_s$cal_sector",histoBuilders.ecsf_chi).fill(pidchi2) 	

			if( pcal_m2 != null ) {
			    histos.computeIfAbsent(p_name+"_"+"pcal_m2u_passall_s$cal_sector",histoBuilders.second_moments).fill(pcal_m2)
			    histos.computeIfAbsent(p_name+"_"+"pcal_m2u_passall",histoBuilders.second_moments).fill(pcal_m2)
			}
			if( eical_m2 != null ){
			    histos.computeIfAbsent(p_name+"_"+"eical_m2u_passall_s$cal_sector",histoBuilders.second_moments).fill(eical_m2)
			    histos.computeIfAbsent(p_name+"_"+"eical_m2u_passall",histoBuilders.second_moments).fill(eical_m2)
			}
			
			if( eocal_m2 != null ){
			    histos.computeIfAbsent(p_name+"_"+"eocal_m2u_passall_s$cal_sector",histoBuilders.second_moments).fill(eocal_m2)
 			    histos.computeIfAbsent(p_name+"_"+"eocal_m2u_passall",histoBuilders.second_moments).fill(eocal_m2)
			}
			//histos.computeIfAbsent(p_name+"_"+"calsf_modsquare_pcalvsecal_passall",histoBuilders.calsf).fill(Math.pow(pcaldep/p,2), Math.pow((eidep + eodep)/p,2))	
			///histos.computeIfAbsent(p_name+"_"+"calsf_modsquare_pcalvsei_passall",histoBuilders.calsf).fill(Math.pow(pcaldep/p,2), Math.pow(eidep/p,2))
			//histos.computeIfAbsent(p_name+"_"+"calsf_modsquare_pcalvseiminuseo_passall",histoBuilders.calsf).fill(Math.pow(pcaldep/p,2), Math.pow(eidep-eodep/p,2))

						
			histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)

			boolean pass_rej_pion_cut = -pcaldep/p + 0.2 > eidep/p 
			if( pass_rej_pion_cut ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall_rej_pions",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_passall_rej_pions",histoBuilders.ecsf).fill(p, edep/p)
			    histos.computeIfAbsent(p_name+"_"+"p_passall_rej_pions",histoBuilders.hp).fill(p)
			}
			else if( !pass_rej_pion_cut ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall_pass_pions",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_passall_pass_pions",histoBuilders.ecsf).fill(p, edep/p)
			    histos.computeIfAbsent(p_name+"_"+"p_passall_pass_pions",histoBuilders.hp).fill(p)
			}

			def pass_good_electron_rej_pion = 0.2 < ((eidep+eodep)/p + pcaldep/p)
			if( pass_good_electron_rej_pion ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_passall_remove_pions",histoBuilders.ecsf).fill(p, edep/p)			    
			}
						    
			histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			histos.computeIfAbsent(p_name+"_"+"calsfmod_pcalvsei_passall",histoBuilders.calsf).fill(pcaldep/p, eidep)

			for( int sig_region = 1; sig_region <= 6; sig_region++ ){			
			    if( Math.abs(pidchi2) <= sig_region ){
				histos.computeIfAbsent(p_name+"_"+"calsf_passall"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p) 	
				histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 			histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)

	 			//if( pcal_m2 != null ) {	histos.computeIfAbsent(p_name+"_"+"pcal_m2_passall"+"_pidchi2_$sig_region"+"sig",histoBuilders.second_moments).fill(pcal_m2) }
	 			//if( eical_m2 != null ) { histos.computeIfAbsent(p_name+"_"+"eical_m2_passall"+"_pidchi2_$sig_region"+"sig",histoBuilders.second_moments).fill(eical_m2) }
	 			//if( eocal_m2 != null ) {histos.computeIfAbsent(p_name+"_"+"eocal_m2_passall"+"_pidchi2_$sig_region"+"sig",histoBuilders.second_moments).fill(eocal_m2) }

				if( pass_rej_pion_cut ){
				    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall_rej_pions_"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				    histos.computeIfAbsent(p_name+"_"+"calsf_passall_rej_pions_"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p)
				    histos.computeIfAbsent(p_name+"_"+"p_passall_rej_pions_"+"pidchi2_$sig_region"+"sig",histoBuilders.hp).fill(p)
				}
				else if( !pass_rej_pion_cut ){
				    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall_pass_pions_"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				    histos.computeIfAbsent(p_name+"_"+"calsf_passall_pass_pions_"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p)
				    histos.computeIfAbsent(p_name+"_"+"p_passall_pass_pions_"+"pidchi2_$sig_region"+"sig",histoBuilders.hp).fill(p)
				}

				//histos.computeIfAbsent(p_name+"_final_"+"calsfchi_prange"+p_range+"_passall_pidchi2"+sig_region,histoBuilders.ecsf_chi).fill(pidchi2) 	

				
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
				histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(eidep/p, eodep/p)
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				histos.computeIfAbsent(p_name+"_"+"calsfmod_pcalvsei_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep/p, eidep)
				histos.computeIfAbsent(p_name+"_"+"mod_pcalvsei_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep)
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
				histos.computeIfAbsent(p_name+"_"+"calsfmod_pcalvseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep/p, eodep)
				histos.computeIfAbsent(p_name+"_"+"mod_pcalvseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eodep)
			    }
			    else{
				histos.computeIfAbsent(p_name+"_"+"calsf_passall"+"_s$cal_sector"+"failpidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p) 	
				histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall"+"_s$cal_sector"+"failpidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 			histos.computeIfAbsent(p_name+"_"+"pcalvsecal_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
				
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
				histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(eidep/p, eodep/p)
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				histos.computeIfAbsent(p_name+"_"+"calsfmod_pcalvsei_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep/p, eidep)
				histos.computeIfAbsent(p_name+"_"+"mod_pcalvsei_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep)
				histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
				histos.computeIfAbsent(p_name+"_"+"calsfmod_pcalvseo_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep/p, eodep)
				histos.computeIfAbsent(p_name+"_"+"mod_pcalvseo_passall"+"failpidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eodep)
			    }
			}

			
			if( p > 4.0 ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_bigp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}
			else if ( p < 4.0 ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    //if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_lowp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    //if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    //if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_lowp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}
		    }
		    else{//if(value.contains(false)){
			histos.computeIfAbsent(p_name+"_"+"calsf_fail_at_least1"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
			histos.computeIfAbsent(p_name+"_"+"calsf_fail_at_least1",histoBuilders.ecsf).fill(p, edep/p) 	
			histos.computeIfAbsent(p_name+"_"+"pcalvsecal_fail_at_least1"+"_s$cal_sector",histoBuilders.eieo).fill( eidep + eodep, pcaldep)
	 		histos.computeIfAbsent(p_name+"_"+"pcalvsecal_fail_at_least1",histoBuilders.eieo).fill(eidep + eodep, pcaldep)
			histos.computeIfAbsent(p_name+"_"+"calsfchi_fail_at_least1",histoBuilders.ecsf_chi).fill(pidchi2) 	
			histos.computeIfAbsent(p_name+"_"+"calsfchi_fail_at_least1_s$cal_sector",histoBuilders.ecsf_chi).fill(pidchi2) 	

			histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			//if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			//if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			//if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)			

			if( p > 4.0 ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    //if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_bigp_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    ///if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    ///if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_bigp_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}
			else if ( p < 4.0 ){
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsecal_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    //if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_lowp_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    //if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    //if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_eivseo_lowp_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvsei_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    histos.computeIfAbsent(p_name+"_"+"calsf_pcalvseo_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}
		    }
		}


		if( ei_sect == eo_sect && ei_sect >= 0 ){			
		    if(!value.contains(false)){
			histos.computeIfAbsent(p_name+"_"+"eieo_passall",histoBuilders.eieo).fill(e_eical, e_eocal)
			histos.computeIfAbsent(p_name+"_"+"eieo_passall"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
			histos.computeIfAbsent(p_name+"_"+"eieo_small_passall",histoBuilders.eieo_small).fill(e_eical, e_eocal)
			histos.computeIfAbsent(p_name+"_"+"eieo_small_passall"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
			//if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent(p_name+"_"+"ecalsf_eivseo_passall",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
			histos.computeIfAbsent(p_name+"_"+"ecalsf_eivseo_passall",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
		    }		   
		    if(value.contains(false)){	    
			histos.computeIfAbsent(p_name+"_"+"eieo_fail_at_least1",histoBuilders.eieo).fill(e_eical, e_eocal)
			histos.computeIfAbsent(p_name+"_"+"eieo_fail_at_least1"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
			histos.computeIfAbsent(p_name+"_"+"eieo_small_fail_at_least1",histoBuilders.eieo_small).fill(e_eical, e_eocal)
			histos.computeIfAbsent(p_name+"_"+"eieo_small_fail_at_least1"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
			//if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent(p_name+"_"+"ecalsf_eivseo_fail_at_least1",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
			histos.computeIfAbsent(p_name+"_"+"ecalsf_eivseo_fail_at_least1",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
		    }		   
		}

		if(!value.contains(false)){		    
		    histos.computeIfAbsent(p_name+"_"+"vz_passall",histoBuilders.vz).fill(vz)
 		    histos.computeIfAbsent(p_name+"_"+"ptheta_passall",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent(p_name+"_"+"p_passall",histoBuilders.hp).fill(p)

		    for( int sig_region = 1; sig_region <= 5; sig_region++ ){			
			if( Math.abs(pidchi2) <= sig_region ){
			    histos.computeIfAbsent(p_name+"_"+"p_passall_pidchi2$sig_region",histoBuilders.hp).fill(p)			    
			}
		    }			   

		    histos.computeIfAbsent(p_name+"_"+"theta_passall",histoBuilders.htheta).fill(theta)		    
		    histos.computeIfAbsent(p_name+"_"+"phitheta_passall",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent(p_name+"_"+"vztheta_passall",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent(p_name+"_"+"vzphi_passall",histoBuilders.vzphi).fill(phi,vz)
		}
		if(value.contains(false)){		    
		    histos.computeIfAbsent(p_name+"_"+"vz_fail_at_least1",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent(p_name+"_"+"ptheta_fail_at_least1",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent(p_name+"_"+"phitheta_fail_at_least1",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent(p_name+"_"+"p_fail_at_least1",histoBuilders.hp).fill(p)
		    histos.computeIfAbsent(p_name+"_"+"theta_fail_at_least1",histoBuilders.htheta).fill(theta)
		    histos.computeIfAbsent(p_name+"_"+"vztheta_fail_at_least1",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent(p_name+"_"+"vzphi_fail_at_least1",histoBuilders.vzphi).fill(phi,vz)
		}


		if( !value.contains(false) ){
		    for( int sig_region = 1; sig_region <= 5; sig_region++ ){			
			if( Math.abs(pidchi2) <= sig_region ){
			    histos.computeIfAbsent(p_name+"_"+"ptheta_pidchi2lvl$sig_region",histoBuilders.ptheta).fill(p, theta)			
			    histos.computeIfAbsent(p_name+"_"+"p_pidchi2lvl$sig_region",histoBuilders.hp).fill(p)
			    histos.computeIfAbsent(p_name+"_"+"theta_pidchi2lvl$sig_region",histoBuilders.htheta).fill(theta)
			}									
		    }
		    histos.computeIfAbsent(p_name+"_"+"allparicles_ptheta",histoBuilders.ptheta).fill(p, theta)
		}	
	    }
	}
	 */		
    }
}
    
