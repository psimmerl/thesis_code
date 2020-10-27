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

    //def rf_offset = [5164:
    //5032: 
    
    //call electron cut constructor
    def field_config = ''
    def electron = new ElectronFromEvent();

    def setCutParameters(){
	println(' [epkpkm::setCutParameters] -> field configuration is ' + field_config)
	println(' [epkpkm::setCutParameters] -> setting cut values for field ' + field_config)
	electron.setElectronCutParameters(field_config)
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
    
    
    def myElectronCutStrategies = [
	electron.passElectronChargeCut,
	//electron.passElectronMinMomentum,
	//electron.passElectronEIEOCut,
	//electron.passElectronNpheCut,
	electron.passElectronVertexCut,
	//electron.passElectronPCALFiducialCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3,
	//electron.passElectronSamplingFractionCut,
	electron.passElectronStatus,
	//electron.passElectronTrackQualityCut
	electron.passElectronEBPIDCut
    ]

    def charge_cut_index=0
    def min_mntm_cut_index=1
    def eieo_cut_index=2
    def nphe_cut_index=3
    def vz_cut_index=4
    def pcalfid_cut_index=5
    //def dcr1_cut_index=6
    //def dcr2_cut_index=7
    //def dcr3_cut_index=8
    def ecsf_cut_index=6
    def status_cut_index=7
    def eb_el_pid_cut_index=9

    /*
     def trck_qual_cut_index=6
    def ecsf_cut_index=11
     */


    
    def histos = new ConcurrentHashMap()
    def histoBuilders = [

	particle_counter : { title -> new H1F("$title","$title", 9.0, 1.0, 10.0 ) },
	ptheta : {title -> new H2F("$title","$title", 900, 0, beam, 900, 0.0, 65.0 ) },
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
	betap : { title -> new H2F("$title","$title", 200, 0.0, 8.0, 200, 0.0, 1.20 ) }

    ]

    def processEvent( event ){
	
	//println("[ElectronMon::processEvent] - processing event now")
	def run_number = event.run_number
	def my_el_cuts = (0..<event.npart).findAll{event.charge[it]<0}.collect{ ii -> [ii, myElectronCutStrategies.collect{ el_test -> el_test(event,ii) } ] }.collectEntries()	
	//println("[ElectronMon::processEvent] - Found Good Electrons")
	//println(my_el_cuts)
	
	my_el_cuts.each{ index, value ->
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

	    }
	    if( event.pcal_status.contains(index) ){
		pcaldep = event.pcal_energy[index]
		cal_sector = event.pcal_sector[index]-1
	    }
	    
	    def edep = eidep + eodep + pcaldep
	    def sf_chi = KinTool.getSamplingFractionFromMomentum(cal_sector, p , edep/p)

	    value.eachWithIndex{ cut_res, cut_ind -> cut_res ? histos.computeIfAbsent("passResults",histoBuilders.passrate).fill(cut_ind) : null }
	    
	    value.eachWithIndex{ cut_result, cut_index ->		
		
		if(cut_result && cut_index==0){
		    if(nphe != null ) {histos.computeIfAbsent("nphe_cut$cut_index",histoBuilders.nphe).fill(nphe)}
		    if(nphe != null ) histos.computeIfAbsent("nphe_s"+cherenkov_sector+"_cut$cut_index",histoBuilders.nphe).fill(nphe)

	 	    if( e_pcal != null ) histos.computeIfAbsent("pcalsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p) 
		    if( pcal_u != null && e_pcal != null ) histos.computeIfAbsent("pcalusf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		    if( pcal_v != null && e_pcal != null ) histos.computeIfAbsent("pcalvsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		    if( pcal_w != null && e_pcal != null ) histos.computeIfAbsent("pcalwsf_cut$cut_index"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)

 		    if( e_eical != null ) histos.computeIfAbsent("eicalsf_cut$cut_index"+"_s$eical_sect",histoBuilders.ecsf).fill(p, e_eical/p) 
		    if( eical_u != null && e_eical != null ) histos.computeIfAbsent("eical_usf_cut$cut_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_u, e_eical/p)
		    if( eical_v != null && e_eical != null ) histos.computeIfAbsent("eical_vsf_cut$cut_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_v, e_eical/p)
		    if( eical_w != null && e_eical != null ) histos.computeIfAbsent("eical_wsf_cut$cut_index"+"_s$eical_sect",histoBuilders.uvwsf).fill(eical_w, e_eical/p)
		    
		    if( eical_u != null && eical_t != null ) histos.computeIfAbsent("eical_ut_cut$cut_index"+"_s$eical_sect",histoBuilders.time).fill(eical_u, eical_t)
		    if( eical_v != null && eical_t != null ) histos.computeIfAbsent("eical_vt_cut$cut_index"+"_s$eical_sect",histoBuilders.time).fill(eical_v, eical_t)
		    if( eical_w != null && eical_t != null ) histos.computeIfAbsent("eical_wt_cut$cut_index"+"_s$eical_sect",histoBuilders.time).fill(eical_w, eical_t)

		    if( eocal_sect >= 0 ){
			histos.computeIfAbsent("eocalsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.ecsf).fill(p, e_eocal/p) 
 			histos.computeIfAbsent("eocal_usf_cut$cut_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_u, e_eocal/p)
			histos.computeIfAbsent("eocal_vsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_v, e_eocal/p)
			histos.computeIfAbsent("eocal_wsf_cut$cut_index"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_w, e_eocal/p)
			
			histos.computeIfAbsent("eocal_ut_cut$cut_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_u, eocal_t)
			histos.computeIfAbsent("eocal_vt_cut$cut_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_v, eocal_t)
			histos.computeIfAbsent("eocal_wt_cut$cut_index"+"_s$eocal_sect",histoBuilders.time).fill(eocal_w, eocal_t)
		    }

 		    if(cal_sector >= 0 ){
			histos.computeIfAbsent("calsfchi_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
			histos.computeIfAbsent("calsfchi_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	
		    			
			for( int sig_region = 1; sig_region < 5; sig_region++ ){			
			    if( Math.abs(pidchi2) < sig_region ){
				histos.computeIfAbsent("calsf_cut$cut_index"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p) 	
				histos.computeIfAbsent("pcalvsecal_cut$cut_index"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 			histos.computeIfAbsent("pcalvsecal_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
			
				histos.computeIfAbsent("calsf_pcalvsecal_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
				if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(eidep/p, eodep/p)
				if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
				if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_cut$cut_index"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    }
			}
			
			histos.computeIfAbsent("calsf_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
			histos.computeIfAbsent("pcalvsecal_cut$cut_index"+"_s$cal_sector",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 		histos.computeIfAbsent("pcalvsecal_cut$cut_index",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
		    
		    
			histos.computeIfAbsent("calsf_pcalvsecal_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)

			def delta_p = 1.0
			def p_min = 0
			def p_max = 1
			def p_range = -1
			for( int p_bin = 1; p_bin < 12; p_bin++ ){
			    p_max = delta_p*p_bin
			    if( p < p_max && p > p_min ){
				p_range = p_range
			    }
			    p_min = p_max
			}

			if( p_range >= 0 ){
			    histos.computeIfAbsent("calsf_pcalvsecal_prange"+"$p_range"+"_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    histos.computeIfAbsent("calsfchi_prange"+"$p_range"+"_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
			    histos.computeIfAbsent("calsfchi_prange"+"$p_range"+"_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	
			    
			    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_prange"+"$p_range"+"_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_prange"+"$p_range"+"_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_prange"+"$p_range"+"_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}
			
 			if( p > 4.0 ){
			    histos.computeIfAbsent("calsf_pcalvsecal_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    histos.computeIfAbsent("calsfchi_bigp_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
			    histos.computeIfAbsent("calsfchi_bigp_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	
			    
			    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_bigp_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_bigp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}
			else if ( p < 4.0 ){
			    histos.computeIfAbsent("calsfchi_smallp_cut$cut_index",histoBuilders.ecsf_chi).fill(sf_chi) 	
			    histos.computeIfAbsent("calsfchi_smallp_cut$cut_index"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	
			    
			    histos.computeIfAbsent("calsf_pcalvsecal_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_lowp_cut$cut_index",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_lowp_cut$cut_index",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			}			    
		    }

		    if( ei_sect == eo_sect && ei_sect >= 0 ){			
			histos.computeIfAbsent("eieo_cut$cut_index",histoBuilders.eieo).fill(e_eical, e_eocal)			    
			histos.computeIfAbsent("eieo_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)			    
			histos.computeIfAbsent("eieo_small_cut$cut_index",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
			histos.computeIfAbsent("eieo_small_cut$cut_index"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)			    
			if( e_eocal > 0 && e_eical > 0 ) histos.computeIfAbsent("ecalsf_eivseo_cut$cut_index",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
		    }

		    if( dc_sector >= 0 ){
			histos.computeIfAbsent("vz_cut$cut_index"+"_s$dc_sector",histoBuilders.vz).fill(vz)
			histos.computeIfAbsent("ptheta_cut$cut_index"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
			histos.computeIfAbsent("phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
			histos.computeIfAbsent("vztheta_cut$cut_index"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
			histos.computeIfAbsent("vxy_cut$cut_index"+"_s$dc_sector",histoBuilders.vxy).fill(vx,vy)
			histos.computeIfAbsent("vzphi_cut$cut_index"+"_s$dc_sector",histoBuilders.vzphi).fill(phi,vz)
			histos.computeIfAbsent("chi_cut$cut_index"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
			histos.computeIfAbsent("ndf_cut$cut_index"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)
		    }
		    histos.computeIfAbsent("vz_cut$cut_index",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_cut$cut_index",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_cut$cut_index",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_cut$cut_index",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vxy_cut$cut_index",histoBuilders.vxy).fill(vx,vy)
		    histos.computeIfAbsent("vzphi_cut$cut_index",histoBuilders.vzphi).fill(phi,vz)
		}
	    }
	    
	    		
	    //plot nphe for tracks passing all cuts
	    if( !value.contains(false) ){
		histos.computeIfAbsent("nphe_passall",histoBuilders.nphe).fill(nphe)
		histos.computeIfAbsent("nphe_passall_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  
	    }		    
	    if( value.contains(false) ){
		if(nphe != null ) histos.computeIfAbsent("nphe_fail_at_least1",histoBuilders.nphe).fill(nphe)
		if(nphe != null ) histos.computeIfAbsent("nphe_fail_at_least1_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  
	    }		    
	    
	    if( everythingElsePassed(value,nphe_cut_index) ){
		if(nphe != null ) histos.computeIfAbsent("nphe_pass_all_but_nphe",histoBuilders.nphe).fill(nphe)
		if(nphe != null ) histos.computeIfAbsent("nphe_pass_all_but_nphe_s"+cherenkov_sector,histoBuilders.nphe).fill(nphe)		  		    		   
	    }

	    if(!value.contains(false)){
		if( e_pcal != null )  histos.computeIfAbsent("pcalsf_passall"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p)
		if( pcal_u != null && e_pcal != null ) histos.computeIfAbsent("pcalusf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		if( pcal_v != null && e_pcal != null ) histos.computeIfAbsent("pcalvsf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		if( pcal_w != null && e_pcal != null ) histos.computeIfAbsent("pcalwsf_passall"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)
	    }		
	    else if(value.contains(false)){
		if( e_pcal != null ) histos.computeIfAbsent("pcalsf_fail_at_least1"+"_s$pcal_sect",histoBuilders.ecsf).fill(p, e_pcal/p)
		if( pcal_u != null && e_pcal != null ) histos.computeIfAbsent("pcalusf_fail_at_least1"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		if( pcal_v != null && e_pcal != null ) histos.computeIfAbsent("pcalvsf_fail_at_least1"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		if( pcal_w != null && e_pcal != null ) histos.computeIfAbsent("pcalwsf_fail_at_least1"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)

	    }		
	    if( everythingElsePassed(value,pcalfid_cut_index) ){
 		if( pcal_u != null && e_pcal != null )  histos.computeIfAbsent("pcalusf_passall_but_pcalfid"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_u, e_pcal/p)
		if( pcal_u != null && e_pcal != null ) histos.computeIfAbsent("pcalvsf_passall_but_pcalfid"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_v, e_pcal/p)
		if( pcal_v != null && e_pcal != null ) histos.computeIfAbsent("pcalwsf_passall_but_pcalfid"+"_s$pcal_sect",histoBuilders.uvwsf).fill(pcal_w, e_pcal/p)		    
	    }


	    if(!value.contains(false)){
		if( e_eical != null ) histos.computeIfAbsent("eicalsf_passall_s$eical_sect",histoBuilders.ecsf).fill(p, e_eical/p) 
		if( eical_u != null && e_eical != null ) histos.computeIfAbsent("eical_usf_passall_s$eical_sect",histoBuilders.uvwsf).fill(eical_u, e_eical/p)
		if( eical_v != null && e_eical != null ) histos.computeIfAbsent("eical_vsf_passall_s$eical_sect",histoBuilders.uvwsf).fill(eical_v, e_eical/p)
		if( eical_w != null && e_eical != null ) histos.computeIfAbsent("eical_wsf_passall_s$eical_sect",histoBuilders.uvwsf).fill(eical_w, e_eical/p)
		
		if( eical_u != null && eical_t != null ) histos.computeIfAbsent("eical_ut_passall_s$eical_sect",histoBuilders.time).fill(eical_u, eical_t)
		if( eical_v != null && eical_t != null ) histos.computeIfAbsent("eical_vt_passall_s$eical_sect",histoBuilders.time).fill(eical_v, eical_t)
		if( eical_w != null && eical_t != null ) histos.computeIfAbsent("eical_wt_passall_s$eical_sect",histoBuilders.time).fill(eical_w, eical_t)
	    }
	    if(value.contains(false)){
		if( e_eical != null ) histos.computeIfAbsent("eicalsf_fail_at_least1_s$eical_sect",histoBuilders.ecsf).fill(p, e_eical/p)
		if( eical_u != null && e_eical != null )  histos.computeIfAbsent("eical_usf_fail_at_least1_s$eical_sect",histoBuilders.uvwsf).fill(eical_u, e_eical/p)
		if( eical_v != null && e_eical != null ) histos.computeIfAbsent("eical_vsf_fail_at_least1_s$eical_sect",histoBuilders.uvwsf).fill(eical_v, e_eical/p)
		if( eical_w != null && e_eical != null )  histos.computeIfAbsent("eical_wsf_fail_at_least1_s$eical_sect",histoBuilders.uvwsf).fill(eical_w, e_eical/p)
		
		if( eical_u != null && eical_t != null ) histos.computeIfAbsent("eical_ut_fail_at_least1_s$eical_sect",histoBuilders.time).fill(eical_u, eical_t)
		if( eical_v != null && eical_t != null ) histos.computeIfAbsent("eical_vt_fail_at_least1_s$eical_sect",histoBuilders.time).fill(eical_v, eical_t)
		if( eical_w != null && eical_t != null ) histos.computeIfAbsent("eical_wt_fail_at_least1_s$eical_sect",histoBuilders.time).fill(eical_w, eical_t)
	    }

	    if( eocal_sect >= 0 ){
 		if(!value.contains(false)){
		    histos.computeIfAbsent("eocalsf_passall_s$eocal_sect",histoBuilders.ecsf).fill(p, e_eocal/p) 
		    histos.computeIfAbsent("eocal_usf_passall"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_u, e_eocal/p)
		    histos.computeIfAbsent("eocal_vsf_passall"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_v, e_eocal/p)
		    histos.computeIfAbsent("eocal_wsf_passall"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_w, e_eocal/p)
		    
		    histos.computeIfAbsent("eocal_ut_passall"+"_s$eocal_sect",histoBuilders.time).fill(eocal_u, eocal_t)
		    histos.computeIfAbsent("eocal_vt_passall"+"_s$eocal_sect",histoBuilders.time).fill(eocal_v, eocal_t)
		    histos.computeIfAbsent("eocal_wt_passall"+"_s$eocal_sect",histoBuilders.time).fill(eocal_w, eocal_t)
		}
 		if(value.contains(false)){
		    histos.computeIfAbsent("eocalsf_fail_at_least1_s$eocal_sect",histoBuilders.ecsf).fill(p, e_eocal/p) 
		    histos.computeIfAbsent("eocal_usf_fail_at_least1"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_u, e_eocal/p)
		    histos.computeIfAbsent("eocal_vsf_fail_at_least1"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_v, e_eocal/p)
		    histos.computeIfAbsent("eocal_wsf_fail_at_least1"+"_s$eocal_sect",histoBuilders.uvwsf).fill(eocal_w, e_eocal/p)
		    
		    histos.computeIfAbsent("eocal_ut_fail_at_least1"+"_s$eocal_sect",histoBuilders.time).fill(eocal_u, eocal_t)
		    histos.computeIfAbsent("eocal_vt_fail_at_least1"+"_s$eocal_sect",histoBuilders.time).fill(eocal_v, eocal_t)
		    histos.computeIfAbsent("eocal_wt_fail_at_least1"+"_s$eocal_sect",histoBuilders.time).fill(eocal_w, eocal_t)
		}
	    }

	    if(true){//cal_sector >= 0 ){
		if(!value.contains(false)){		    
		    histos.computeIfAbsent("calsf_passall"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
		    histos.computeIfAbsent("pcalvsecal_passall"+"_s$cal_sector",histoBuilders.eieo).fill( eidep + eodep, pcaldep)
	 	    histos.computeIfAbsent("pcalvsecal_passall",histoBuilders.eieo).fill(eidep + eodep, pcaldep)
		    histos.computeIfAbsent("calsfchi_passall",histoBuilders.ecsf_chi).fill(sf_chi) 	
		    histos.computeIfAbsent("calsfchi_passall_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	
		 
   
		    histos.computeIfAbsent("calsf_pcalvsecal_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    
		    for( int sig_region = 1; sig_region < 5; sig_region++ ){			
			if( Math.abs(pidchi2) < sig_region ){
			    histos.computeIfAbsent("calsf_passall"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.ecsf).fill(p, edep/p) 	
			    histos.computeIfAbsent("pcalvsecal_passall"+"_s$cal_sector"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
	 		    histos.computeIfAbsent("pcalvsecal_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep + eodep)
			    
			    histos.computeIfAbsent("calsf_pcalvsecal_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(eidep/p, eodep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsfmod_pcalvsei_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep/p, eidep)
			    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("mod_pcalvsei_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eidep)
			    if( pcaldep > 0 && eodep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
			    if( pcaldep > 0 && eodep > 0 ) histos.computeIfAbsent("calsfmod_pcalvseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep/p, eodep)
			    if( pcaldep > 0 && eodep > 0 ) histos.computeIfAbsent("mod_pcalvseo_passall"+"pidchi2_$sig_region"+"sig",histoBuilders.eieo).fill(pcaldep, eodep)
			}
		    }

		    def delta_p = 1.0
		    def p_min = 0
		    def p_max = 1
		    def p_range = -1
		    for( int p_bin = 1; p_bin < 12; p_bin++ ){
			p_max = delta_p*p_bin
			if( p < p_max && p > p_min ){
			    p_range = p_range
			}
			p_min = p_max
		    }

		    if( p_range >= 0 ){
			histos.computeIfAbsent("calsf_pcalvsecal_prange"+"$p_range"+"_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			histos.computeIfAbsent("calsfchi_prange"+"$p_range"+"_passall",histoBuilders.ecsf_chi).fill(sf_chi) 	
			histos.computeIfAbsent("calsfchi_prange"+"$p_range"+"_passall"+"_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	
			
			if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_prange"+"$p_range"+"_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_prange"+"$p_range"+"_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_prange"+"$p_range"+"_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    }


		    if( p > 4.0 ){
			histos.computeIfAbsent("calsf_pcalvsecal_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_bigp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_bigp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    }
		    else if ( p < 4.0 ){
			histos.computeIfAbsent("calsf_pcalvsecal_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_lowp_passall",histoBuilders.calsf).fill(eidep/p, eodep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_lowp_passall",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    }
		}
		else if(value.contains(false)){
		    histos.computeIfAbsent("calsf_fail_at_least1"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
		    histos.computeIfAbsent("pcalvsecal_fail_at_least1"+"_s$cal_sector",histoBuilders.eieo).fill( eidep + eodep, pcaldep)
	 	    histos.computeIfAbsent("pcalvsecal_fail_at_least1",histoBuilders.eieo).fill(eidep + eodep, pcaldep)
		    histos.computeIfAbsent("calsfchi_fail_at_least1",histoBuilders.ecsf_chi).fill(sf_chi) 	
		    histos.computeIfAbsent("calsfchi_fail_at_least1_s$cal_sector",histoBuilders.ecsf_chi).fill(sf_chi) 	

		    histos.computeIfAbsent("calsf_pcalvsecal_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    
		    if( p > 4.0 ){
			histos.computeIfAbsent("calsf_pcalvsecal_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_bigp_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_bigp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    }
		    else if ( p < 4.0 ){
			histos.computeIfAbsent("calsf_pcalvsecal_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
			if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_lowp_fail_at_least1",histoBuilders.calsf).fill(eidep/p, eodep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
			if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_lowp_fail_at_least1",histoBuilders.calsf).fill(pcaldep/p, eodep/p)
		    }
		}
		if( everythingElsePassed(value,eieo_cut_index) ){
		    histos.computeIfAbsent("calsf_passall_but_eieo"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
		    histos.computeIfAbsent("pcalvsecal_passall_but_eieo_"+"_s$cal_sector",histoBuilders.eieo).fill(eidep + eodep,pcaldep)
	 	    histos.computeIfAbsent("pcalvsecal_passall_but_eieo",histoBuilders.eieo).fill(eidep + eodep, pcaldep)
		    histos.computeIfAbsent("calsf_pcalvsecal_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)

		    histos.computeIfAbsent("calsf_pcalvsecal_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall_but_eieo",histoBuilders.calsf).fill(eidep/p, eodep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall_but_eieo",histoBuilders.calsf).fill(pcaldep/p, eodep/p)

		}


		if( everythingElsePassed(value,ecsf_cut_index) ){
		    histos.computeIfAbsent("calsf_passall_but_ecsf"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
		    histos.computeIfAbsent("pcalvsecal_passall_but_ecsf_"+"_s$cal_sector",histoBuilders.eieo).fill(eidep + eodep,pcaldep)
	 	    histos.computeIfAbsent("pcalvsecal_passall_but_ecsf",histoBuilders.eieo).fill(eidep + eodep, pcaldep)

		    histos.computeIfAbsent("calsf_pcalvsecal_passall_but_ecsf",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall_but_ecsf",histoBuilders.calsf).fill(eidep/p, eodep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall_but_ecsf",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall_but_ecsf",histoBuilders.calsf).fill(pcaldep/p, eodep/p)

		}


		if( everythingElsePassed(value,min_mntm_cut_index) ){
		    histos.computeIfAbsent("cal_sf_passall_but_minmntm"+"_s$cal_sector",histoBuilders.ecsf).fill(p, edep/p) 	
		    histos.computeIfAbsent("pcalvsecal_passall_but_minmntm_"+"_s$cal_sector",histoBuilders.eieo).fill( eidep + eodep, pcaldep)
	 	    histos.computeIfAbsent("pcalvsecal_passall_but_minmntm",histoBuilders.eieo).fill(eidep + eodep, pcaldep)

		    histos.computeIfAbsent("calsf_pcalvsecal_passall_but_minmntm",histoBuilders.calsf).fill(pcaldep/p, (eidep + eodep)/p)			    
		    if( eodep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_eivseo_passall_but_minmntm",histoBuilders.calsf).fill(eidep/p, eodep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvsei_passall_but_minmntm",histoBuilders.calsf).fill(pcaldep/p, eidep/p)
		    if( pcaldep > 0 && eidep > 0 ) histos.computeIfAbsent("calsf_pcalvseo_passall_but_minmntm",histoBuilders.calsf).fill(pcaldep/p, eodep/p)	
		}
	    }


	    if( ei_sect == eo_sect && ei_sect >= 0 ){			
		if(!value.contains(false)){
		    histos.computeIfAbsent("eieo_passall",histoBuilders.eieo).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_passall"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_small_passall",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_small_passall"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		    if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent("ecalsf_eivseo_passall",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
		}		   
		if(value.contains(false)){	    
		    histos.computeIfAbsent("eieo_fail_at_least1",histoBuilders.eieo).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_fail_at_least1"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_small_fail_at_least1",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_small_fail_at_least1"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		    if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent("ecalsf_eivseo_fail_at_least1",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
		}		   
		if( everythingElsePassed(value,eieo_cut_index) ){
		    histos.computeIfAbsent("eieo_passall_but_eieo",histoBuilders.eieo).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_passall_but_eieo"+"_s$ei_sect",histoBuilders.eieo).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_small_passall_but_eieo",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		    histos.computeIfAbsent("eieo_small_passall_but_eieo"+"_s$ei_sect",histoBuilders.eieo_small).fill(e_eical, e_eocal)
		    if(e_eocal > 0 && e_eical > 0) histos.computeIfAbsent("ecalsf_eivseo_passall_but_eieo",histoBuilders.calsf).fill(e_eical/p, e_eocal/p)
		}
	    }

	    if( dc_sector >= 0 ){
		if(!value.contains(false)){		    
		    histos.computeIfAbsent("vz_passall"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_passall"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_passall",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_passall"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_passall",histoBuilders.vzphi).fill(phi,vz)
 		    histos.computeIfAbsent("chi_passall"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 		    histos.computeIfAbsent("ndf_passall"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		
		}
		if(value.contains(false)){		    
		    histos.computeIfAbsent("vz_fail_at_least1"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("ptheta_fail_at_least1"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		    histos.computeIfAbsent("phitheta_fail_at_least1",histoBuilders.phitheta).fill(phi, theta)
		    histos.computeIfAbsent("vztheta_fail_at_least1"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_fail_at_least1",histoBuilders.vzphi).fill(phi,vz)
 		    histos.computeIfAbsent("chi_fail_at_least1"+"_s$dc_sector",histoBuilders.tchi2).fill(dc_chi2)
 		    histos.computeIfAbsent("ndf_fail_at_least1"+"_s$dc_sector",histoBuilders.tndf).fill(dc_ndf)		
 		    
		}

		if( everythingElsePassed(value,vz_cut_index) ){
		    histos.computeIfAbsent("vz_passall_passall_but_vz"+"_s$dc_sector",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("vztheta_passall_but_vz"+"_s$dc_sector",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vz_passall_passall_but_vz",histoBuilders.vz).fill(vz)
		    histos.computeIfAbsent("vztheta_passall_but_vz",histoBuilders.vztheta).fill(vz, theta)
		    histos.computeIfAbsent("vzphi_passall_but_vz",histoBuilders.vzphi).fill(phi,vz)
		}
		
		if( everythingElsePassed(value,min_mntm_cut_index) ){
		    histos.computeIfAbsent("ptheta_passall_but_minmntm"+"_s$dc_sector",histoBuilders.ptheta).fill(p, theta)
		}

		if( !value.contains(false) ){
		    if( pid == 211 || pid == -211 ){
			histos.computeIfAbsent("pions_ptheta",histoBuilders.ptheta).fill(p, theta)						
			for( int sig_region = 1; sig_region < 5; sig_region++ ){			
			    if( Math.abs(pidchi2) < sig_region ){
				histos.computeIfAbsent("pions_ptheta_pidchi2lvl$sig_region",histoBuilders.ptheta).fill(p, theta)			
			    }						
			}
		    }
		    else{
			for( int sig_region = 1; sig_region < 5; sig_region++ ){			
			    if( Math.abs(pidchi2) < sig_region ){
				histos.computeIfAbsent("electrons_ptheta_pidchi2lvl$sig_region",histoBuilders.ptheta).fill(p, theta)			
 			    }									    
			    histos.computeIfAbsent("electrons_ptheta",histoBuilders.ptheta).fill(p, theta)
			}
		    }
		    histos.computeIfAbsent("allparicles_ptheta",histoBuilders.ptheta).fill(p, theta)
		}

	    }
	    
	}
    }
	    
	
	    
    /*
 	    if( event.dc1_status.contains(index) ){		    
		def hits = event.dc1[index]
		def hit = hits.find{it.layer==12}.each{			
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result && cut_index==0){
  			    if( cut_result ){ histos.computeIfAbsent("dcr1hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y) }
			}
		    }		    
		    if(!value.contains(false)){histos.computeIfAbsent("dcr1hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
		}		
	    }
 	    
	    if( event.dc2_status.contains(index) ){		    
		def hits = event.dc2[index]
		def hit = hits.find{it.layer==24}.each{
		    value.eachWithIndex{ cut_result, cut_index ->
			//println(cut_result + ' ' + cut_index)
			if(cut_result && cut_index==0){
			    if( cut_result ){ histos.computeIfAbsent("dcr2hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y)}
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr2hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}		
		    }
		}
	    }

	    if( event.dc3_status.contains(index) ){		    
		def hits = event.dc3[index]
		def hit = hits.find{it.layer==36}.each{
		    value.eachWithIndex{ cut_result, cut_index ->
			if(cut_result && cut_index==0){
			    if( cut_result ){ histos.computeIfAbsent("dcr3hit_cut$cut_index",histoBuilders.hitxy).fill(it.x, it.y) }
			}
			if(!value.contains(false)){histos.computeIfAbsent("dcr3hit_passall",histoBuilders.hitxy).fill(it.x, it.y)}
		    }
		}
	    }
     */

    
}
