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

class HadronMonV2{

    def beam = 10.604
    def rfPeriod = 4.008
    def rf_large_integer = 1000
    //call electron cut constructor
    def electron = new ElectronFromEvent();

    def myElectronCutStrategies = [
	electron.passElectronStatus,
	electron.passElectronChargeCut,
	electron.passElectronEBPIDCut,
	//electron.passElectronNpheCut,
	//electron.passElectronVertexCut
	//electron.passElectronPCALFiducialCut,
	//electron.passElectronEIEOCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3
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
	deltat : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -2, 2  ) },
	deltab : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, -1, 1  ) },
	betap : {title -> new H2F("$title","$title", 500, 0.0, beam, 500, 0.0, 1.4  ) },
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
	iron_cross : { title -> new H2F("$title","$title", 500, -20.0, 20.0, 500, -20.0, 20.0 ) }

	

    ]

    def processEvent( event, hevent ){

	//println( ' processing event ' )
	def start_time = event.start_time
	def rf_time = event.rf_time
	if( start_time != null && rf_time != null ){
	    histos.computeIfAbsent("start_time",histoBuilders.time).fill(start_time)
	    histos.computeIfAbsent("rf_time",histoBuilders.time).fill(rf_time)
	}
		
	def electron_candidates = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	  
	
	def ev_sum_tof_time = 0
	def ev_sum_tof_time_pion_region = 0
	def n_tracks = 0
	def n_tracks_pion_region = 0
 	/*def time_hadrons = electron_candidates.findResults{el_candidate -> ! el_candidate.value.contains(false) ? el_candidate.key : null }.findResults{ indx ->
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    	    
	    return [indx, ele]
	}.collectMany{ indx, ele -> 
	    //look at all charged hadrons
	    (0..<event.npart).findAll{ event.charge[it] != 0 && ( it != indx )}.findResults{index->		
		def vt = event.vt[index]
		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		def p = lv.mag()
		def beta = event.beta[index]

		event.tof[index].eachWithIndex{ hh, ii ->  
		
		    def ftof_sector = hh.sector
		    def ftof_t = hh.time
 		    def ftof_p = hh.path
		    def ftof_paddle = hh.paddle
		    def ftof_layer = hh.layer
		    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)
		    def delta_t = ftof_t - vt
		    if( ftof_layer == 2){
			ev_sum_tof_time += ftof_t - vt
			n_tracks+=1
		    }		    
		    if( ftof_layer == 2 && p < 1.0 && beta > 0.95 ){
			ev_sum_tof_time_pion_region += ftof_t - vt
			n_tracks_pion_region+=1			
		    }
		    if(  ftof_layer == 2 && p < 1.0 && beta > 0.95 && delta_t > 20.28 ){
			histos.computeIfAbsent("ftof_RECbetap_pion_region_deltat_cut",histoBuilders.betap).fill(p, beta ) 			
		    }

		}				
	    }
	}

	if( n_tracks != 0 ){
	    histos.computeIfAbsent("ev_sum_tof_time",histoBuilders.sum_time).fill(ev_sum_tof_time)
	    histos.computeIfAbsent("ev_sum_tof_time_averaged_over_ntracks",histoBuilders.sum_time).fill(ev_sum_tof_time/n_tracks)
	    histos.computeIfAbsent("ev_sum_tof_time_vs_ntracks",histoBuilders.sum_time2d).fill(ev_sum_tof_time, n_tracks)
	}

	if( n_tracks_pion_region != 0 ){
	    histos.computeIfAbsent("ev_sum_tof_time_pion_region",histoBuilders.sum_time).fill(ev_sum_tof_time_pion_region)
	    histos.computeIfAbsent("ev_sum_tof_time_pion_region_averaged_over_ntracks",histoBuilders.sum_time).fill(ev_sum_tof_time_pion_region/n_tracks_pion_region)
	}
	if ( true) {//ev_sum_tof_time < 86.616 + 5 && ev_sum_tof_time > 86.616 - 5){
	    //histos.computeIfAbsent("ev_sum_tof_time_cut",histoBuilders.sum_time).fill(ev_sum_tof_time)
	 */

 	def hadrons = electron_candidates.findResults{el_candidate -> ! el_candidate.value.contains(false) ? el_candidate.key : null }.findResults{ indx ->
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    	    
	    return [indx, ele]
	}.collectMany{ indx, ele -> 
	    //look at positive charged hadrons

	    (0..<event.npart).findAll{event.charge[it] > 0 && ( it != indx )}.findResults{index->		
		def lv = new Vector3(event.px[index], event.py[index], event.pz[index])
		def pid = event.pid[index]
		def p = lv.mag()
		def theta = Math.toDegrees(lv.theta())
		def phi = Math.toDegrees(lv.phi())
		def beta = event.beta[index]
		def vz = event.vz[index]
		def start_time_rf_corrected = event.vt[index] // start time with vertex correction
		def elvz = event.vz[indx]
		def vt = event.vt[index]
 		def run_number = event.run_number
		def delta_vz_abs = Math.abs(vz - elvz)
		def delta_vz = (vz - elvz)
		
		if( event.dc1_status.contains(index) ){		    
		    def hits = event.dc1[index]
		    def hit = hits.find{it.layer==12}.each{
 			histos.computeIfAbsent("dcr1hit_pos",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}

		if( event.dc2_status.contains(index) ){		    
		    def hits = event.dc2[index]
		    def hit = hits.find{it.layer==24}.each{
 			histos.computeIfAbsent("dcr2hit_pos",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }
		}

		if( event.dc3_status.contains(index) ){		    
		    def hits = event.dc3[index]
		    def hit = hits.find{it.layer==36}.each{
 			histos.computeIfAbsent("dcr3hit_pos",histoBuilders.hitxy).fill(it.x, it.y) 		    		
		    }		    
		}
		
		if (event.dc3_status.contains(index) ){		    
 		    def dc_sector = event.dc_sector[index]		    
		    histos.computeIfAbsent("vz_p_pos_s"+dc_sector,histoBuilders.vz_p).fill(p,vz)
		    histos.computeIfAbsent("vz_theta_pos_s"+dc_sector,histoBuilders.vz_theta).fill(theta,vz)
		    histos.computeIfAbsent("vz_phi_pos_s"+dc_sector,histoBuilders.vz_phi).fill(phi,vz)		    
		    
		    histos.computeIfAbsent("delta_vz_p_pos_s"+dc_sector,histoBuilders.delta_vz_p).fill(p,delta_vz)
		    histos.computeIfAbsent("delta_vz_theta_pos_s"+dc_sector,histoBuilders.delta_vz_theta).fill(theta,delta_vz)
		    histos.computeIfAbsent("delta_vz_phi_pos_s"+dc_sector,histoBuilders.delta_vz_phi).fill(phi,delta_vz)	    
		    histos.computeIfAbsent("delta_vz_beta_pos_s"+dc_sector,histoBuilders.delta_vz_phi).fill(delta_vz,beta)	    
		    
		    if(pid == 2212){
			histos.computeIfAbsent("vz_p_pro_s"+dc_sector,histoBuilders.vz_p).fill(p,vz)
			histos.computeIfAbsent("vz_theta_pro_s"+dc_sector,histoBuilders.vz_theta).fill(theta,vz)
			histos.computeIfAbsent("vz_phi_pro_s"+dc_sector,histoBuilders.vz_phi).fill(phi,vz)		    

			histos.computeIfAbsent("delta_vz_p_pos_s"+dc_sector,histoBuilders.delta_vz_p).fill(p,delta_vz)
			histos.computeIfAbsent("delta_vz_theta_pos_s"+dc_sector,histoBuilders.delta_vz_theta).fill(theta,delta_vz)
			histos.computeIfAbsent("delta_vz_phi_pos_s"+dc_sector,histoBuilders.delta_vz_phi).fill(phi,delta_vz)    
		    }
		    else if( pid == 321 ){
			histos.computeIfAbsent("vz_p_kp_s"+dc_sector,histoBuilders.vz_p).fill(p,vz)
			histos.computeIfAbsent("vz_theta_kp_s"+dc_sector,histoBuilders.vz_theta).fill(theta,vz)
			histos.computeIfAbsent("vz_phi_kp_s"+dc_sector,histoBuilders.vz_phi).fill(phi,vz)		    

			histos.computeIfAbsent("delta_vz_p_kp_s"+dc_sector,histoBuilders.delta_vz_p).fill(p,delta_vz)    
			histos.computeIfAbsent("delta_vz_theta_kp_s"+dc_sector,histoBuilders.delta_vz_theta).fill(theta,delta_vz)    
			histos.computeIfAbsent("delta_vz_phi_kp_s"+dc_sector,histoBuilders.delta_vz_phi).fill(phi,delta_vz)    
		    }
		    else if( pid == 211 ){
			histos.computeIfAbsent("vz_p_pip_s"+dc_sector,histoBuilders.vz_p).fill(p,vz)
			histos.computeIfAbsent("vz_theta_pip_s"+dc_sector,histoBuilders.vz_theta).fill(theta,vz)
			histos.computeIfAbsent("vz_phi_pip_s"+dc_sector,histoBuilders.vz_phi).fill(phi,vz)		    

			histos.computeIfAbsent("delta_vz_p_pip_s"+dc_sector,histoBuilders.delta_vz_p).fill(p,delta_vz)
			histos.computeIfAbsent("delta_vz_theta_pip_s"+dc_sector,histoBuilders.delta_vz_theta).fill(theta,delta_vz)
			histos.computeIfAbsent("delta_vz_phi_pip_s"+dc_sector,histoBuilders.delta_vz_phi).fill(phi,delta_vz)		    
		    }
		}
		
		//get electron tof info
		if(event.tof_status.contains(index) ){		  
		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[index].layer
		    def my_layer = []
		    
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }
		    if( my_layer.contains('2') && !my_layer.contains('3') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }		    
		    event.tof[index].eachWithIndex{ hh, ii ->  
			if(layer_to_use!=-1 && layer_to_use == hh.layer && hh.layer!=3 && hh.layer!=0 ){// hh.layer == 2 ){
			    //println('sector ' + hh.sector + ' paddle ' + hh.paddle + ' layer ' + hh.layer + ' time  ' + hh.time + ' path  ' + hh.path + ' energy ' + hh.energy)
			    def measured_beta = beta// ftof_p/(ftof_t - start_time)/(29.9792458)
			    def shiftp= 0
			    def shift=0
			    if(hh.layer == 2 ){  
				shiftp=0
				shift=0 // -4
			    }
			    def ftof_sector = hh.sector
 			    def ftof_t = hh.time
			    def ftof_p = hh.path+shiftp
			    def ftof_paddle = hh.paddle
			    def ftof_layer = hh.layer
			    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    		   
			    def measured_t  = ftof_t - vt


 			    def calc_beta_pr = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(2212)**2)
			    def calc_time_pr = ftof_p/(calc_beta_pr*(29.9792458))
			    
			    def calc_beta_kp = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(321)**2)
			    def calc_time_kp = ftof_p/(calc_beta_kp*(29.9792458))
			    
			    def calc_beta_pip = p/Math.sqrt(p**2 + PDGDatabase.getParticleMass(211)**2)
			    def calc_time_pip = ftof_p/(calc_beta_pip*(29.9792458))
			    
			    def delta_beta_pr = calc_beta_pr - measured_beta
			    def delta_t_pr = calc_time_pr - measured_t
			    
			    def delta_beta_kp = calc_beta_kp - measured_beta
			    def delta_t_kp = calc_time_kp - measured_t
			    
			    def delta_beta_pip = calc_beta_pip - measured_beta
			    def delta_t_pip = calc_time_pip - measured_t

			    //////////////////////////////////
			    // delta t 
			    //
			    histos.computeIfAbsent("ftof_deltat_pr",histoBuilders.deltat).fill(p, delta_t_pr)
			    histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pr)
			    histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_pr)
			    histos.computeIfAbsent("ftof_deltat_pr_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pr)

			    histos.computeIfAbsent("ftof_deltat_kp",histoBuilders.deltat).fill(p, delta_t_kp)
			    histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_kp)
			    histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_kp)
			    histos.computeIfAbsent("ftof_deltat_kp_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_kp)

			    histos.computeIfAbsent("ftof_deltat_pip",histoBuilders.deltat).fill(p, delta_t_pip)
			    histos.computeIfAbsent("ftof_deltat_pip_s$ftof_sector",histoBuilders.deltat).fill(p, delta_t_pip)
			    histos.computeIfAbsent("ftof_deltat_pip_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltat).fill(p, delta_t_pip)
			    histos.computeIfAbsent("ftof_deltat_pip_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltat).fill(p, delta_t_pip)


			    //delta beta
			    histos.computeIfAbsent("ftof_deltab_pr",histoBuilders.deltab).fill(p, delta_beta_pr)
			    histos.computeIfAbsent("ftof_deltab_pr_s$ftof_sector",histoBuilders.deltab).fill(p, delta_beta_pr)
			    histos.computeIfAbsent("ftof_deltab_pr_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltab).fill(p, delta_beta_pr)
			    histos.computeIfAbsent("ftof_deltab_pr_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pr)

			    histos.computeIfAbsent("ftof_deltab_kp",histoBuilders.deltab).fill(p, delta_beta_kp)
			    histos.computeIfAbsent("ftof_deltab_kp_s$ftof_sector",histoBuilders.deltab).fill(p, delta_beta_kp)
			    histos.computeIfAbsent("ftof_deltab_kp_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltab).fill(p, delta_beta_kp)
			    histos.computeIfAbsent("ftof_deltab_kp_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_kp)

			    histos.computeIfAbsent("ftof_deltab_pip",histoBuilders.deltab).fill(p, delta_beta_pip)
			    histos.computeIfAbsent("ftof_deltab_pip_s$ftof_sector",histoBuilders.deltab).fill(p, delta_beta_pip)
			    histos.computeIfAbsent("ftof_deltab_pip_s$ftof_sector"+"_layer$ftof_layer",histoBuilders.deltab).fill(p, delta_beta_pip)
			    histos.computeIfAbsent("ftof_deltab_pip_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle",histoBuilders.deltab).fill(p, delta_beta_pip)

			
			    histos.computeIfAbsent("ftof_rec_betap_pos",histoBuilders.betap).fill(p, measured_beta) 			    	
  			    histos.computeIfAbsent("ftof_rec_betap_pos"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 			    histos.computeIfAbsent("ftof_rec_betap_pos"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta)
 			    histos.computeIfAbsent("ftof_rec_betap_pos_layer$layer_to_use",histoBuilders.betap).fill(p, measured_beta)
 			    histos.computeIfAbsent("ftof_rec_betap_pos"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
 			    histos.computeIfAbsent("ftof_rec_betap_pos"+"_layer$layer_to_use"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
			    histos.computeIfAbsent("ftof_rec_betap_pos"+"_s$ftof_sector"+"_layer$layer_to_use",histoBuilders.betap).fill(p, measured_beta)
  			    histos.computeIfAbsent("vz_vt_pos"+"_s$ftof_sector"+"_layer$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  			    histos.computeIfAbsent("rf_vt_pos"+"_s$ftof_sector"+"_layer$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)
  			    histos.computeIfAbsent("vz_vt_pos"+"_s$ftof_sector"+"_layer$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  			    histos.computeIfAbsent("rf_vt_pos"+"_s$ftof_sector"+"_layer$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)
			    

			    if(pid == 2212){
 				def pr_FTOF_vt = ftof_t - ftof_p / ( calc_beta_pr * 29.98f ) - start_time
				def thisTime = ftof_t -  ftof_p / ( calc_beta_pr * 29.98f ) - rf_time
				thisTime = (thisTime+(rf_large_integer+0.5f)*rfPeriod) % rfPeriod
				thisTime = thisTime - rfPeriod/2
 				histos.computeIfAbsent("mon_time_proton_l$ftof_layer",histoBuilders.mon_time).fill(thisTime,pr_FTOF_vt)
			    }
			    if(pid == 321){
 				def kp_FTOF_vt = ftof_t - ftof_p / ( calc_beta_kp * 29.98f ) - start_time
				def thisTime = ftof_t -  ftof_p / ( calc_beta_kp * 29.98f ) - rf_time
				thisTime = (thisTime+(rf_large_integer+0.5f)*rfPeriod) % rfPeriod
				thisTime = thisTime - rfPeriod/2
 				histos.computeIfAbsent("mon_time_kaonPlus_l$ftof_layer",histoBuilders.mon_time).fill(thisTime,kp_FTOF_vt)
			    }
			    if(pid == 211){
 				def pip_FTOF_vt = ftof_t - ftof_p / ( calc_beta_pip * 29.98f ) - start_time
				def thisTime = ftof_t -  ftof_p / ( calc_beta_pip * 29.98f ) - rf_time
				thisTime = (thisTime+(rf_large_integer+0.5f)*rfPeriod) % rfPeriod
				thisTime = thisTime - rfPeriod/2
 				histos.computeIfAbsent("mon_time_pionPlus_l$ftof_layer",histoBuilders.mon_time).fill(thisTime,pip_FTOF_vt)
			    }
			    
			    if( vt > 80.0 && vt < 94 ){
 				histos.computeIfAbsent("ftof_rec_betap_pos_vtcut_layer$layer_to_use"+"_elvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_rec_betap_pos_elvzcut_vtcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_rec_betap_pos_elvzcut_vtcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
  				histos.computeIfAbsent("vz_vt_vtcut_layer$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  				histos.computeIfAbsent("rf_vt_vtcut_layer$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)				
			    }
			    if( vt > 80.0 && vt < 94 && vz < 10 && vz > -12 ){
 				histos.computeIfAbsent("ftof_rec_betap_pos_vtcut_vzcut_layer$layer_to_use"+"_elvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_rec_betap_pos_elvzcut_vtcut_vzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_RECbetap_pos_elvzcut_vtcut_vzcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
  				histos.computeIfAbsent("vz_vt_vtcut_vzcut_layer$layer_to_use",histoBuilders.vz_vt).fill(vz, start_time_rf_corrected)
  				histos.computeIfAbsent("rf_vt_vtcut_vzcut_layer$layer_to_use",histoBuilders.rf_vt).fill(rf_time, start_time_rf_corrected)				
			    }			
			    else{
				histos.computeIfAbsent("ftof_rec_betap_pos_layer$layer_to_use"+"_BADelvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_rec_betap_pos_BADelvzcut",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_rec_betap_pos_BADelvzcut"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta)
			    }				
 			
			    histos.computeIfAbsent("ftof_betap_pos_layer$layer_to_use",histoBuilders.betap).fill(p, measured_beta)
 			    histos.computeIfAbsent("ftof_betap_pos",histoBuilders.betap).fill(p, measured_beta)
 			    
			
			    
			    
 			    
			    /*def measured_t = ftof_t - start_time 
			    def measured_beta2 = ftof_p/(ftof_t - start_time_rf_corrected + shift)/(29.9792458)
			    def measured_beta2shift = ftof_p/(ftof_t - start_time - 2)/(29.9792458)
			     if( true ){// p < 1.0 && measured_beta2 > 0.91 ){				
  				if( measured_beta > 1.1 && measured_beta < 1.23 && p < 2.0 && Math.abs(measured_beta - measured_beta2) > 0.1 && Math.abs(vz-elvz) < 5.0 ){
				}				
   				histos.computeIfAbsent("ftof_mybetap_pos"+"_s$ftof_sector"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos"+"_s$ftof_sector",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
 				histos.computeIfAbsent("ftof_mybetap_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)

 				if( Math.abs(vz-elvz) < 5.0 ){
				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_elvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_elvzcut"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				    histos.computeIfAbsent("ftof_mybetap_pos_elvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)

				}
				else{
				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_BADelvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_BADelvzcut"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
  				    histos.computeIfAbsent("ftof_mybetap_pos_BADelvzcut"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betap).fill(p, measured_beta2)
				    histos.computeIfAbsent("ftof_betap_vs_mybeta_pos_BADelvzcut"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)

				}
				histos.computeIfAbsent("ftof_betap_vs_mybeta_pos"+"_layer$ftof_layer"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)
 				histos.computeIfAbsent("ftof_betap_vs_mybeta_pos"+"_r$run_number",histoBuilders.betacomp).fill(measured_beta, measured_beta2)

				    
			 }
			    
			     */
						
			}
		    }
		}
	    }
	}


	def timing_hadrons = electron_candidates.findResults{el_candidate -> ! el_candidate.value.contains(false) ? el_candidate.key : null }.findResults{ indx ->
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    	    
	    return [indx, ele]
	}.collectMany{ indx, ele -> 
	    (0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipr->
		def delta_t_pr=-1000
		def lv = new Vector3(event.px[ipr], event.py[ipr], event.pz[ipr])
		def pid = event.pid[ipr]
		def p_pr = lv.mag()
		def vt = event.vt[ipr]
		def beta_pr=event.beta[ipr]

		if(event.tof_status.contains(ipr) ){		  
		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[ipr].layer
		    def my_layer = []
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }

		    if( my_layer.contains('2') && !my_layer.contains('3') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }		    
		    

		    event.tof[ipr].eachWithIndex{ hh, ii ->  
			//println('proton tof ' )
			//println(ii)
			//println('layer is -')
			//println(hh.layer)
			//println('layer to use os ')
			//println(layer_to_use)
			if(layer_to_use!=-1 && layer_to_use == hh.layer && hh.layer!=3 && hh.layer!=0 ){
			    //println('PASSSEEDDDDD')
			    def ftof_sector = hh.sector
			    def ftof_t = hh.time
			    def ftof_p = hh.path
			    def ftof_paddle = hh.paddle
			    def ftof_layer = hh.layer
			    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    	
			    def calc_beta_pr = p_pr/Math.sqrt(p_pr**2 + PDGDatabase.getParticleMass(2212)**2)
			    def measured_t  = ftof_t - vt
			    def calc_time_pr = ftof_p/(calc_beta_pr*(29.9792458))
			    delta_t_pr = calc_time_pr - measured_t
			    //println(delta_t_pr)
			}
		    }
		}
		//println('final t pr ')
		//println(delta_t_pr)
		return [delta_t_pr, beta_pr, p_pr]
	    }
	}.collectMany{delta_t_pr, beta_pr, p_pr->
	    (0..<event.npart).findAll{event.pid[it]==211}.findResults{ikp->
		def delta_t_kp=-1000	 
		def lv = new Vector3(event.px[ikp], event.py[ikp], event.pz[ikp])
		def pid = event.pid[ikp]
		def p_kp = lv.mag()
		def vt = event.vt[ikp]
		def beta_kp = event.beta[ikp]

		if(event.tof_status.contains(ikp) ){		  
		    def layer_to_use = -1
		    //println(' layers available ' + event.tof[index].layer)
		    def layer_avail = event.tof[ikp].layer
		    def my_layer = []
		    for( ii in layer_avail){
			my_layer.add(Integer.toString(ii))
		    }


		    if( my_layer.contains('2') && !my_layer.contains('3') ){
			//println('use layer 2 for timing ' )
			layer_to_use = 2
		    }
		    else if ( !my_layer.contains('2') && my_layer.contains('1') && !my_layer.contains('3')  ){
			//println('use layer 1 for timing ' )
			layer_to_use = 1
		    }		    

		    event.tof[ikp].eachWithIndex{ hh, ii ->  
			if(layer_to_use!=-1 && layer_to_use == hh.layer && hh.layer!=3 && hh.layer!=0 ){
			    def ftof_sector = hh.sector
			    def ftof_t = hh.time
			    def ftof_p = hh.path
			    def ftof_paddle = hh.paddle
			    def ftof_layer = hh.layer
			    def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    	
			    def calc_beta_kp = p_kp/Math.sqrt(p_kp**2 + PDGDatabase.getParticleMass(211)**2)
			    def measured_t  = ftof_t - vt
			    def calc_time_kp = ftof_p/(calc_beta_kp*(29.9792458))
			    delta_t_kp = calc_time_kp - measured_t
			}
		    }
		}
		return [delta_t_pr,delta_t_kp, beta_pr, beta_kp, p_pr, p_kp]
	    }
	}.findResults{delta_t_pr, delta_t_kp, beta_pr, beta_kp, p_pr, p_kp->
	    return [delta_t_pr, delta_t_kp, beta_pr, beta_kp, p_pr, p_kp]
	}

	
	
	if( timing_hadrons.size() > 0 ){
	    for(th in timing_hadrons){
		//println(th)
		def (delta_t_pr,delta_t_kp, beta_pr, beta_kp, p_pr, p_kp)=th
		//println(delta_t_pr)
		//println(delta_t_kp)
		//println(beta_pr)
		//println(beta_kp)
		if( delta_t_pr != -1000 && delta_t_kp != -1000 ){
		    histos.computeIfAbsent("iron_cross_pipx_pry_deltat",histoBuilders.iron_cross).fill(delta_t_kp,delta_t_pr)
		    histos.computeIfAbsent("iron_cross_pip_rec_betap",histoBuilders.betap).fill(p_kp, beta_kp) 			    	
		    histos.computeIfAbsent("iron_cross_pr_rec_betap",histoBuilders.betap).fill(p_pr, beta_pr)
 		    if( delta_t_kp > 3 ){
			histos.computeIfAbsent("iron_cross_pip_rec_betap_shifted",histoBuilders.betap).fill(p_kp, beta_kp) 			    	
			histos.computeIfAbsent("iron_cross_pr_rec_betap_shifted",histoBuilders.betap).fill(p_pr, beta_pr)
		    }
		}
	    }	    
	}
	   
    }
}
