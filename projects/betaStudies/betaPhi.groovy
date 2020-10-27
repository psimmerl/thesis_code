package projects.betaStudies
import org.jlab.detector.base.DetectorType
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentLinkedQueue
import pid.electron.ElectronFromEvent
import utils.PhysicsEvent
import event.Event
import org.jlab.jroot.ROOTFile
import org.jlab.jroot.TNtuple
import java.util.concurrent.atomic.AtomicInteger;


class betaPhi {

    def histos = new ConcurrentHashMap()
    def beam = LorentzVector.withPID(11,0,0,10.604)
    def target = LorentzVector.withPID(2212,0,0,0)
   

    def hcombs = {new H1F("$it","$it",20,1,19)}
    def histoBuilders = [
	deltat : {title -> new H2F("$title","$title", 500, 0.0, beam.e(), 500, -2, 2  ) },
	deltab : {title -> new H2F("$title","$title", 500, 0.0, beam.e(), 500, -1, 1  ) },
	betap : {title -> new H2F("$title","$title", 500, 0.0, beam.e(), 500, 0.0, 1.4  ) },
	deltab_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -1, 1  ) },
	deltat_time2d : {title -> new H2F("$title","$title", 100, 0.0, 100, 500, -2, 2  ) }
    ]

    Random random = new Random()

    def electron = new ElectronFromEvent();
    def phyev = new PhysicsEvent()        
    def f1 = null 
    def t1 = null
    def data_type = 0
    //private Integer global_event=0
    private static  AtomicInteger at = new AtomicInteger(0);
    

     def myElectronCutStrategies = [
	 electron.passElectronStatus,
	 electron.passElectronChargeCut,
	 electron.passElectronEBPIDCut,
	 electron.passElectronNpheCut,
	 electron.passElectronVertexCut,
	 electron.passElectronPCALFiducialCut,
	 electron.passElectronEIEOCut,
	 electron.passElectronDCR1,
	 electron.passElectronDCR2,
	 electron.passElectronDCR3	
     ]
 
    def getGlobalEvent(glob_ev){
	global_event=glob_ev
    }

     def createROOTFile(fout){
	 f1 = new ROOTFile(fout)
	 println(' DATA TYPE OF OUTPUT FILE IS ' + data_type )
	 if( this.data_type == 1 ){
	     t1 = f1.makeNtuple('clas12','clas12phi',phyev.genTreeString()+":"+phyev.genMCTreeString()+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb")	    
	 }
	 else{
	     t1 = f1.makeNtuple('clas12','clas12phi',phyev.genTreeString()+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb")	    
	 }

     }

     def fillOutTree(int el_index, int el_status, LorentzVector el, int pr_index, int pr_status, LorentzVector pr, int kp_index, int kp_status, LorentzVector kp, int km_index, int km_status, LorentzVector km, int hel, int glob_ev, int comb_index, double prb, double kpb, double kmb ){
	t1.fill(el_index,
		el_status,
		el.px(),
		el.py(),
		el.pz(),
		el.e(),
		pr_index,
		pr_status,
		pr.px(),
		pr.py(),
		pr.pz(),
		pr.e(),
		kp_index,
		kp_status,
		kp.px(),
		kp.py(),
		kp.pz(),
		kp.e(),
		km_index,
		km_status,
		km.px(),
		km.py(),
		km.pz(),
		km.e(),
		hel,
		glob_ev,
		comb_index,
		prb,
		kpb,
		kmb
	)

    }

    def closeOutTree(){
	t1.write()
	f1.close()
    }

    def processEvents(event) {		
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()

	//get electron proton pair first
	def ep = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    
	    return [index:indx, lv:ele, name:'el', status:event.status[indx] ]
	    //return ele
	}.collectMany{ mele -> 
	    (0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 	
		def mpro = [index:ipro, lv:pro, name:'pro', status:event.status[ipro]]
		return [mele,mpro]
	    }
	}

	if( ep.size() == 1 ){	    
	    def combs = ep.collectMany{mele,mpro->            	
		(0..<event.npart).findAll{event.pid[it]==321}.findResults{ikp->
		    def kp = LorentzVector.withPID(321,event.px[ikp], event.py[ikp], event.pz[ikp]) 		
		    def mkp = [index:ikp, lv:kp,name:'kp',status:event.status[ikp]]
		    return [mele,mpro,mkp]
		}
	    }.collectMany{mele,mpro,mkp->  
		(0..<event.npart).findAll{event.pid[it]==-321}.findResults{ikm->
		    def km = LorentzVector.withPID(-321,event.px[ikm], event.py[ikm], event.pz[ikm]) 						
		    def mkm = [index:ikm, lv:km, name:'km',status:event.status[ikm]]
		    return [mele,mpro,mkp,mkm]
		}
	    }.findResults{ele,pro,kp,km->	    
		//returns a list of map elements
		return [ele,pro,kp,km]
	    }

	    histos.computeIfAbsent('combsize',hcombs).fill(combs.size())
	    if( combs.size() > 0 ){
		def cc=0
		def global_event = at.incrementAndGet()
		for(comb in combs){
		    def (el, pro, kp, km ) = comb		    
		    		    
		    def part_betas = [pro:0,kp:0,km:0]
		    comb.eachWithIndex{it, i->
			if( i!=0 ){
			    def index = it.get('index')
			    def part_lv = it.get('lv')
			    def part_name = it.get('name')
			    def stat_val = (it.get('status'))
			    def part_status = (stat_val >= 2000 && stat_val < 4000 ) ? 'FD' : 'CD'
			    def part_beta = event.beta[index]
			    part_betas.put(part_name,part_beta)
			    def p = part_lv.p()
			
			    def measured_betarec = part_beta			    			    
			    def calc_betarec = p/Math.sqrt(p**2 + part_lv.mass()**2)			    
			    def delta_betarec = calc_betarec - measured_betarec
			    
			    histos.computeIfAbsent("betap_$part_name"+"_det$part_status",histoBuilders.betap).fill(p, measured_betarec)
 			    histos.computeIfAbsent("betap_$part_name"+"_det$part_status",histoBuilders.betap).fill(p, measured_betarec)
			    histos.computeIfAbsent("betap_$part_name"+"_det$part_status",histoBuilders.betap).fill(p, measured_betarec)
			    histos.computeIfAbsent("betap_$part_name"+"_det$part_status",histoBuilders.betap).fill(p, measured_betarec)

			    
			    if( event.tof_status.contains(index) ){		
 				//println( part_status)
				def start_time = event.start_time
				def rf_time = event.rf_time
				
				def ftof_sector = event.tof_sector[index]
				def ftof_t = event.tof_time[index]
				def ftof_p = event.tof_path[index]
				def ftof_paddle = event.tof_paddle[index]
				def ftof_layer = event.tof_layer[index]
				def ftof_t_corr = ftof_t - ftof_p/(29.9792458)		    		   
				
				def measured_beta = ftof_p/(ftof_t - start_time)/(29.9792458)
				def measured_t = ftof_t - start_time

 
				def calc_beta = p/Math.sqrt(p**2 + part_lv.mass()**2)
				def calc_time = ftof_p/(calc_beta*(29.9792458))

				def delta_beta = calc_beta - measured_beta
				def delta_t = calc_time - measured_t

				
 				histos.computeIfAbsent("ftof_betap"+"_layer$ftof_layer"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
				histos.computeIfAbsent("ftof_betap"+"_s$ftof_sector"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
				histos.computeIfAbsent("ftof_betap"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
				
				histos.computeIfAbsent("ftof_betap"+"_s$ftof_sector"+"_layer$ftof_layer"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
				histos.computeIfAbsent("ftof_betap_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_deltab_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"_det$part_status",histoBuilders.deltab).fill(p, delta_beta)
				//histos.computeIfAbsent("ftof_deltat_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"_det$part_status",histoBuilders.deltat).fill(p, delta_t)

 				histos.computeIfAbsent("ftof_deltab_$part_name"+"_s$ftof_sector"+"_det$part_status",histoBuilders.deltab).fill(p, delta_beta)
				//histos.computeIfAbsent("ftof_deltat_$part_name"+"_s$ftof_sector"+"_det$part_status",histoBuilders.deltat).fill(p, delta_t)
			
 				histos.computeIfAbsent("ftof_betap"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_betap_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle"+"_det$part_status",histoBuilders.betap).fill(p, measured_beta)
 				histos.computeIfAbsent("ftof_deltab_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle"+"_det$part_status",histoBuilders.deltab).fill(p, delta_beta)
				//histos.computeIfAbsent("ftof_deltat_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"_paddle$ftof_paddle"+"_det$part_status",histoBuilders.deltat).fill(p, delta_t)

				//////////////// 
				// per paddle 
				histos.computeIfAbsent("ftof_t_deltab_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle"+"_det$part_status",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_beta)
				//histos.computeIfAbsent("ftof_t_deltat_$part_name"+"_s$ftof_sector"+"_layer$ftof_layer"+"per_paddle"+"_det$part_status",histoBuilders.deltab_time2d).fill(ftof_paddle,delta_t)
			    }
			}
		    }

		    def el_index = el.get('index')
		    def el_lv = el.get('lv')
		    def el_stat_val = el.get('status')
		    def pr_index = pro.get('index')
		    def pr_lv = pro.get('lv')
		    def pr_stat_val = pro.get('status')
		    def kp_index = kp.get('index')
		    def kp_lv = kp.get('lv')
		    def kp_stat_val = kp.get('status')
		    def km_index = km.get('index')
		    def km_lv = km.get('lv')
		    def km_stat_val = km.get('status')
		    fillOutTree(el_index,el_stat_val,el_lv,pr_index,pr_stat_val,pr_lv,kp_index,kp_stat_val,kp_lv,km_index,km_stat_val,km_lv,event.helicity,global_event,cc,part_betas.get('pro'),part_betas.get('kp'),part_betas.get('km'))
		    cc+=1


		}
		
	    }
	    
	}
    }


}
