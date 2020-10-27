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


class ElasticMon{

    def beam = 10.604
    def target_position = -3.0 //cm
    def speed_of_light = 29.9792458 // cm/ns
    def lv_beam = LorentzVector.withPID(11,0,0,10.604)
    def target = LorentzVector.withPID(2212,0,0,0)
    
    def el_pid = 11
    def pim_pid = -211
    def el_cut_strictness=[]

    //call electron cut constructor
    def field_config = ''
    def electron = new ElectronFromEvent();
    
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

    def getGeneratedSector(phi){
	return Math.ceil(phi / 60) + 3
    }
    

    def getKin(beam, target, electron) {

	def missing = new Particle(beam)
	missing.combine(target, 1)
	missing.combine(electron, -1)
	def w = missing.mass()

	def q = new Particle(beam)
	q.combine(electron, -1)
	def q2 = -1 * q.mass2()

	def nu = beam.e() - electron.e()
	def y = nu / beam.e()
	def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

	return [x: x, y: y, w: w, nu: nu, q2: q2]
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



    def histos = new ConcurrentHashMap()
    def histoBuilders = [
	
	particle_counter : { title -> new H1F("$title","$title", 9.0, 1.0, 10.0 ) },
	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
	hp : {title -> new H1F("$title","$title", 100, 0, 10.00) },
	htheta : {title -> new H1F("$title","$title", 900, 0, 65.0) },
	htheta_small : {title -> new H1F("$title","$title", 40, 5.0, 15.0) },
	htheta_med : {title -> new H1F("$title","$title", 20, 5.0, 15.0) },
	phitheta : {title -> new H2F("$title","$title", 900, 180, -180, 900, 0.0, 65.0 ) },
	hw : { title -> new H1F("$title","$title", 100, 0.0, 1.50)},
	hw_small : { title -> new H1F("$title","$title", 150, 0.5, 1.50)},
	wphi : { title -> new H2F("$title","$title", 600, -180, 180, 600, 0.0, 1.50 ) },
	wtheta : {title -> new H2F("$title","$title", 20, 5.0, 15.0, 600, 0.0, 1.50 ) },
	hq2x : { title -> new H2F("$title","$title", 500, 0.0, 1.1, 500, 0.0, beam)},
	hq2 : { title -> new H1F("$title","$title", 500, 0.0, beam)},
	hx : { title -> new H1F("$title","$title", 500, 0.0, 1.1)},
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
	fcup_gated : { title -> new H1F("$title","$title", 10, 0, 10 )},
	fcup_ungated : { title -> new H1F("$title","$title", 10, 0, 10 )},
	fcup_tot : { title -> new H1F("$title","$title", 10, 0, 10 )},
		      

    ]

        def processEvent( event ){
	
 	///println("[ElectronMon::processEvent] - processing event now")
	def run_number = event.run_number 
	def fc_ungated = event.fcup_ungated
	def fc_gated = event.fcup_gated
	def event_number = event.event_number
	//&& passPCALFiducialCut(event,it)
	def el_cand = (0..<event.npart).findAll{ event.charge[it]<0 && event.pid[it] == 11 && Math.abs(event.vz[it]+3.0) < 5  && (Math.abs(event.status[it]) >= 2000 && Math.abs(event.status[it]) < 4000) && passPCALFiducialCut(event,it) }.collect{ ii -> [ii, [(event.pid[ii] == el_pid )] ] }.collectEntries()	
	
	//println("[ElectronMon::processEvent] - Found Good Electrons")
	
	el_cand.each{ index, value ->
	    def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
	    
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
	    
	    def cal_sector = -1
	    def pcaldep = null
	    if( event.pcal_status.contains(index) ){
 		pcaldep = event.pcal_energy[index]
		cal_sector = event.pcal_sector[index]-1
	    }
	    				
	    def is_elastic = eX.mass() < 1.2
	    if(!value.contains(false) && is_elastic ){
	
		histos.computeIfAbsent('hw_el_pass_all',histoBuilders.hw).fill(eX.mass())
		histos.computeIfAbsent('hw_small_el_pass_all',histoBuilders.hw_small).fill(eX.mass())
		histos.computeIfAbsent('hq2_pass_all',histoBuilders.hq2).fill(q2)
		histos.computeIfAbsent('hx_pass_all',histoBuilders.hx).fill(xb)
		histos.computeIfAbsent('hq2x_pass_all',histoBuilders.hq2x).fill(xb,q2)
		histos.computeIfAbsent('htheta_small_el_pass_all',histoBuilders.htheta_small).fill(theta)
		histos.computeIfAbsent('htheta_med_el_pass_all',histoBuilders.htheta_med).fill(theta)
		histos.computeIfAbsent('htheta_el_pass_all',histoBuilders.htheta).fill(theta)
		histos.computeIfAbsent('hptheta_el_pass_all',histoBuilders.ptheta).fill(p,theta)
		histos.computeIfAbsent('hphitheta_el_pass_all',histoBuilders.phitheta).fill(phi,theta)
		histos.computeIfAbsent('hphivz_el_pass_all',histoBuilders.vzphi).fill(phi,vz)
		histos.computeIfAbsent('hvztheta_el_pass_all',histoBuilders.vztheta).fill(vz,theta)
		histos.computeIfAbsent('hphiw_el_pass_all',histoBuilders.wphi).fill(phi,eX.mass())
		histos.computeIfAbsent('hwtheta_el_pass_all',histoBuilders.wtheta).fill(theta,eX.mass())

		histos.computeIfAbsent('hw_el_pass_all_' + cal_sector,histoBuilders.hw).fill(eX.mass())
		histos.computeIfAbsent('hw_small_el_pass_all_' + cal_sector,histoBuilders.hw_small).fill(eX.mass())
		histos.computeIfAbsent('htheta_small_el_pass_all' + cal_sector,histoBuilders.htheta_small).fill(theta)
		histos.computeIfAbsent('htheta_med_el_pass_all' + cal_sector,histoBuilders.htheta_med).fill(theta)
		histos.computeIfAbsent('htheta_el_pass_all' + cal_sector,histoBuilders.htheta).fill(theta)
		histos.computeIfAbsent('hptheta_el_pass_all' + cal_sector,histoBuilders.ptheta).fill(p,theta)
		histos.computeIfAbsent('hphitheta_el_pass_all' + cal_sector,histoBuilders.phitheta).fill(phi,theta)
		histos.computeIfAbsent('hphivz_el_pass_all' + cal_sector,histoBuilders.vzphi).fill(phi,vz)
		histos.computeIfAbsent('hvztheta_el_pass_all' + cal_sector,histoBuilders.vztheta).fill(vz,theta)
		histos.computeIfAbsent('hphiw_el_pass_all' + cal_sector,histoBuilders.wphi).fill(phi,eX.mass())
		histos.computeIfAbsent('hwtheta_el_pass_all_' + cal_sector,histoBuilders.wtheta).fill(theta,eX.mass())

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// per run
		histos.computeIfAbsent('hw_el_pass_all_r'+run_number,histoBuilders.hw).fill(eX.mass())
		histos.computeIfAbsent('hw_small_el_pass_all_r'+run_number,histoBuilders.hw_small).fill(eX.mass())
		histos.computeIfAbsent('hq2_pass_all_r'+run_number,histoBuilders.hq2).fill(q2)
		histos.computeIfAbsent('hx_pass_all_r'+run_number,histoBuilders.hx).fill(xb)
		histos.computeIfAbsent('hq2x_pass_all_r'+run_number,histoBuilders.hq2x).fill(xb,q2)
		histos.computeIfAbsent('htheta_small_el_pass_all_r'+run_number,histoBuilders.htheta_small).fill(theta)
		histos.computeIfAbsent('htheta_med_el_pass_all_r'+run_number,histoBuilders.htheta_med).fill(theta)
		histos.computeIfAbsent('htheta_el_pass_all_r'+run_number,histoBuilders.htheta).fill(theta)
		histos.computeIfAbsent('hptheta_el_pass_all_r'+run_number,histoBuilders.ptheta).fill(p,theta)
		histos.computeIfAbsent('hphitheta_el_pass_all_r'+run_number,histoBuilders.phitheta).fill(phi,theta)
		histos.computeIfAbsent('hphivz_el_pass_all_r'+run_number,histoBuilders.vzphi).fill(phi,vz)
		histos.computeIfAbsent('hvztheta_el_pass_all_r'+run_number,histoBuilders.vztheta).fill(vz,theta)
		histos.computeIfAbsent('hphiw_el_pass_all_r'+run_number,histoBuilders.wphi).fill(phi,eX.mass())
		histos.computeIfAbsent('hwtheta_el_pass_all_' + run_number,histoBuilders.wtheta).fill(theta,eX.mass())

		histos.computeIfAbsent('hw_el_pass_all_' + cal_sector + '_r'+run_number,histoBuilders.hw).fill(eX.mass())
		histos.computeIfAbsent('hw_small_el_pass_all_' + cal_sector + '_r'+run_number,histoBuilders.hw_small).fill(eX.mass())
		histos.computeIfAbsent('htheta_small_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.htheta_small).fill(theta)
		histos.computeIfAbsent('htheta_med_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.htheta_med).fill(theta)
		histos.computeIfAbsent('htheta_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.htheta).fill(theta)
		histos.computeIfAbsent('hptheta_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.ptheta).fill(p,theta)
		histos.computeIfAbsent('hphitheta_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.phitheta).fill(phi,theta)
		histos.computeIfAbsent('hphivz_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.vzphi).fill(phi,vz)
		histos.computeIfAbsent('hvztheta_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.vztheta).fill(vz,theta)
		histos.computeIfAbsent('hphiw_el_pass_all' + cal_sector + '_r'+run_number,histoBuilders.wphi).fill(phi,eX.mass())
		histos.computeIfAbsent('hwtheta_el_pass_all_' + cal_sector + '_r'+run_number,histoBuilders.wtheta).fill(theta,eX.mass())
 		
	    }
	}

	(0 ..< event.mc_npart).find{
	    event.mc_pid[it] == 11
	}?.each{ idx -> 
	    
	    def ele = LorentzVector.withPID(11, event.mc_px[idx], event.mc_py[idx], event.mc_pz[idx])


	    def p = ele.p()
	    def vz = event.mc_vz[idx]
	    def theta = Math.toDegrees(ele.theta())


	   


	    def phi = Math.toDegrees(ele.phi())
	    def pidchi2 = event.chi2pid[idx]
	    def vx = event.mc_vx[idx]
	    def vy = event.mc_vy[idx]
	    def gen_sector = getGeneratedSector(phi)

	    def eX = lv_beam+target-ele
	    def q = lv_beam - ele
	    def q2 = -q.mass2()
	    def xb = q2/(2*target.mass()*q.e())
	    def nu = q2/(2*target.mass()*xb)
	    def y = nu/lv_beam.e()	    
		    
	    histos.computeIfAbsent('h_gen_w_el_pass_all',histoBuilders.hw).fill(eX.mass())
	    histos.computeIfAbsent('h_gen_w_small_el_pass_all',histoBuilders.hw_small).fill(eX.mass())
	    histos.computeIfAbsent('h_gen_q2_pass_all',histoBuilders.hq2).fill(q2)
	    histos.computeIfAbsent('h_gen_x_pass_all',histoBuilders.hx).fill(xb)
	    histos.computeIfAbsent('h_gen_q2x_pass_all',histoBuilders.hq2x).fill(xb,q2)
	    histos.computeIfAbsent('h_gen_theta_small_el_pass_all',histoBuilders.htheta_small).fill(theta)
	    histos.computeIfAbsent('h_gen_theta_med_el_pass_all',histoBuilders.htheta_med).fill(theta)
	    histos.computeIfAbsent('h_gen_theta_el_pass_all',histoBuilders.htheta).fill(theta)
	    histos.computeIfAbsent('h_gen_ptheta_el_pass_all',histoBuilders.ptheta).fill(p,theta)
	    histos.computeIfAbsent('h_gen_phitheta_el_pass_all',histoBuilders.phitheta).fill(phi,theta)
	    histos.computeIfAbsent('h_gen_phivz_el_pass_all',histoBuilders.vzphi).fill(phi,vz)
	    histos.computeIfAbsent('h_gen_vztheta_el_pass_all',histoBuilders.vztheta).fill(vz,theta)
	    histos.computeIfAbsent('h_gen_phiw_el_pass_all',histoBuilders.wphi).fill(phi,eX.mass())
	    histos.computeIfAbsent('h_gen_wtheta_el_pass_all',histoBuilders.wtheta).fill(theta,eX.mass())
	    
	    histos.computeIfAbsent('h_gen_w_el_pass_all_' + gen_sector,histoBuilders.hw).fill(eX.mass())
	    histos.computeIfAbsent('h_gen_w_small_el_pass_all_' + gen_sector,histoBuilders.hw_small).fill(eX.mass())
	    histos.computeIfAbsent('h_gen_theta_small_el_pass_all' + gen_sector,histoBuilders.htheta_small).fill(theta)
	    histos.computeIfAbsent('h_gen_theta_med_el_pass_all' + gen_sector,histoBuilders.htheta_med).fill(theta)
	    histos.computeIfAbsent('h_gen_theta_el_pass_all' + gen_sector,histoBuilders.htheta).fill(theta)
	    histos.computeIfAbsent('h_gen_ptheta_el_pass_all' + gen_sector,histoBuilders.ptheta).fill(p,theta)
	    histos.computeIfAbsent('h_gen_phitheta_el_pass_all' + gen_sector,histoBuilders.phitheta).fill(phi,theta)
	    histos.computeIfAbsent('h_gen_phivz_el_pass_all' + gen_sector,histoBuilders.vzphi).fill(phi,vz)
	    histos.computeIfAbsent('h_gen_vztheta_el_pass_all' + gen_sector,histoBuilders.vztheta).fill(vz,theta)
	    histos.computeIfAbsent('h_gen_phiw_el_pass_all' + gen_sector,histoBuilders.wphi).fill(phi,eX.mass())
	    histos.computeIfAbsent('h_gen_wtheta_el_pass_all_' + gen_sector,histoBuilders.wtheta).fill(theta,eX.mass())
	    	   
	}

    }
}


