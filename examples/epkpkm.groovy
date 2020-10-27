package examples
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


class epkpkm {
  
    def phyev = new PhysicsEvent()        
    def f1 = null 
    def t1 = null
    def field_config=''
    def el_cut_strictness=[]

    //def beam = 10.604
    def target_position = -3.0 //cm
    def speed_of_light = 29.9792458 // cm/ns
    def rfBucketLength = 4.0 //ns

    //set_data_type -> 0 for DATA and 1 for SIMULATION
    def data_type = -1
    def setDataType(dt){
	data_type = dt
    }

    def hists = new ConcurrentHashMap()
    def beam = LorentzVector.withPID(11,0,0,10.604)
    def target = LorentzVector.withPID(2212,0,0,0)
    def phim = 1.019
    def hmissemm2 = {new H2F("$it","$it",100,-2,5,100,0,5)}
    def hmissephim = {new H2F("$it","$it",100,-2,5,100,0,5)}
    def hlambdaphi = {new H2F("$it","$it",100,1.4,2.0,100,0.95,1.15)}
    def hcpl = {new H1F("$it","$it",100,0.0,50.0)}
    def hpt = {new H1F("$it","$it",100,0.0,2.0)}
    def hphim = {new H1F("$it","$it",75,0.8,1.5)}
    def hcombs = {new H1F("$it","$it",9,1,10)}
    def hW = {new H1F("$it","$it",100,0,5)}
    def hmissmm2phim = {new H2F("$it","$it",100,-1,1,100,0,5)}
    def hmissmm2 = {new H1F("$it","$it",500,-1,1.1)}
    def hpele = {new H1F("$it","$it",75,0.0,10.0)}
    def hptheta = {new H2F("$it","$it",75,0.0,10.0,75,0.0,60)}
    def hbeta = {new H2F("$it","it",100, 0.0, 5.5, 100, 0.0, 1.05)}
    def hang = {new H1F("$it","$it",100,0,60)}
    def hphitheta = {new H2F("$it","$it",200,-180,180,200,0.0,60.0)}
    def hq2x = {new H2F("$it","$it",75,0.0,1.0,75,0.0,beam.e()+1) }
    def ht = {new H1F("$it","$it",50,0.0,5.0) }
    def hq = {new H1F("$it","$it",50,0.0,beam.e()+1) }
    def hx = {new H1F("$it","$it",50,0.0,1.0) }
    def htxb = {new H2F("$it","$it",50,0.0,5.0,50,0.0,1.1)}
    def hq2t = {new H2F("$it","$it",50,0.0,5.0,50,0.0,beam.e()+1)}
    def hcos = {new H1F("$it","$it",15,-1.0,1.0)}
    def hasy = {new H1F("$it","$it",20,0.0,360.0)}
    def hhel = {new H1F("$it","$it",4,-2,2)}
    def htime = {new H1F("$it","$it",150,-50,200)}    
    def hbetap = {new H2F("$it","$it",200,0.0,beam.e(),200,0.0,1.5)}

    Random random = new Random()
    //private Integer global_event=0
    private static  AtomicInteger at = new AtomicInteger(0);

   
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

    def myElectronCutStrategies = [
	electron.passElectronStatus,
	electron.passElectronChargeCut,
	//electron.passElectronTrackQualityCut,
	//electron.passElectronMinMomentum,
	electron.passElectronEBPIDCut,
	//electron.passElectronSamplingFractionCut,
	//electron.passElectronNpheCut,
	//electron.passElectronVertexCut,
	//electron.passElectronPCALFiducialCut,
	//electron.passElectronEIEOCut,
	//electron.passElectronDCR1,
	//electron.passElectronDCR2,
	//electron.passElectronDCR3
	//electron.passElectronAntiPionCut
    ]

    def createROOTFile(fout){
	phyev.setAnalysisTopology("epkpkm")

	def mc_tree_string  = phyev.genMCTreeString()
	def tree_string  = phyev.genTreeString()
	//println( mc_tree_string)

	f1 = new ROOTFile(fout)
	println(' DATA TYPE OF OUTPUT FILE IS ' + data_type )
	if( this.data_type == 1 ){
 	    println(tree_string+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb:"+mc_tree_string)
	    t1 = f1.makeNtuple('clas12','clas12phi',tree_string+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb:"+mc_tree_string)
	}
	else{
	    println(tree_string+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb") 
	    t1 = f1.makeNtuple('clas12','clas12phi',tree_string+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb") 
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

    def fillMCOutTree (int el_index, int el_status, LorentzVector el, int pr_index, int pr_status, LorentzVector pr, int kp_index, int kp_status, LorentzVector kp, int km_index, int km_status, LorentzVector km, int hel, int glob_ev, int comb_index, double prb, double kpb, double kmb, LorentzVector mcel, LorentzVector mcpr, LorentzVector mckp, LorentzVector mckm){
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
		kmb,
		mcel.px(),
		mcel.py(),
		mcel.pz(),
		mcel.e(),
		mcpr.px(),
		mcpr.py(),
		mcpr.pz(),
		mcpr.e(),
		mckp.px(),
		mckp.py(),
		mckp.pz(),
		mckp.e(),
		mckm.px(),
		mckm.py(),
		mckm.pz(),
		mckm.e()	
	)

    }

    def closeOutTree(){
	println('[epkpkm::closeOutTree] -> writing tree')
	t1.write()
	println('[epkpkm::closeOutTree] -> done writing tree')
	println('[epkpkm::closeOutTree] -> closing file now')
	f1.close()
	println('[epkpkm::closeOutTree] -> done closing file')
    }


    def processEventTest(event){
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()

	//get electron proton pair first
	def ep = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
 	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    
	    return [el_index:indx, lv:ele]
	    //return ele
	}.collectMany{ mele -> 
	    (0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 	
		def  mpro = [pr_index:ipro, lv:pro]
		return [mele,mpro]
            }
	}
	
	println(ep)

    }
   
    def processEventComb(event) {		
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()

	def final_el = new LorentzVector(0,0,0,0)
	def final_pr = new LorentzVector(0,0,0,0)
	def final_kp = new LorentzVector(0,0,0,0)
	def final_km = new LorentzVector(0,0,0,0)

	def final_mcel = new LorentzVector(0,0,0,0)
	def final_mcpr = new LorentzVector(0,0,0,0)
	def final_mckp = new LorentzVector(0,0,0,0)
	def final_mckm = new LorentzVector(0,0,0,0)

	def event_rec = 0

	//get MC::Particle Variables 
	if(data_type == 1){
	    if( event.mc_status){
		final_mcel = LorentzVector.withPID(11,event.mc_px[0],event.mc_py[0],event.mc_pz[0])
		final_mcpr = LorentzVector.withPID(2212,event.mc_px[1],event.mc_py[1],event.mc_pz[1])
		final_mckp = LorentzVector.withPID(321,event.mc_px[2],event.mc_py[2],event.mc_pz[2])
		final_mckm = LorentzVector.withPID(-321,event.mc_px[3],event.mc_py[3],event.mc_pz[3])
	    }
	}

	//get electron proton pair first
	def ep = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
 	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    
	    return ele
	}.collectMany{ ele -> 
	    (0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 	
		return [ele,pro]
            }
	}
	
	if( ep.size() == 1 ){	    
	    def combs = ep.collectMany{ele,pro->            	
		(0..<event.npart).findAll{event.pid[it]==321}.findResults{ikp->
 		    def kp = LorentzVector.withPID(321,event.px[ikp], event.py[ikp], event.pz[ikp]) 		
		    return [ele,pro,kp]
		}
            }.collectMany{ele,pro,kp->  
		(0..<event.npart).findAll{event.pid[it]==-321}.findResults{ikm->
		    def km = LorentzVector.withPID(-321,event.px[ikm], event.py[ikm], event.pz[ikm]) 						
		    return [ele,pro,kp,km]
		}
            }.findResults{ele,pro,kp,km->	    
		def epkpkmX = beam+target-ele-pro-kp-km	    

		def ekpkmX = beam+target-ele-kp-km
		def epkpX = beam+target-ele-pro-kp
		def epkmX = beam+target-ele-pro-km
		
		def cplpro = ekpkmX.vect().theta(pro.vect())
		def cplkm = epkpX.vect().theta(km.vect())
		def cplkp = epkmX.vect().theta(kp.vect())
		def pt = Math.sqrt(epkpkmX.px()**2 + epkpkmX.py()**2)

		def eX = beam+target-ele
		def epX = beam+target-ele-pro
		def lambda = pro+kp
		def vphi = kp+km
		def q = beam-ele
		def q2 = -q.mass2()	    
		def xb = q2/(2*target.mass()*q.e())
		def t = 2*target.mass()*(pro.e()-target.mass())
		def nu = q2/(2*target.mass()*xb)
		def y = nu/beam.e()

 		hists.computeIfAbsent('cplpro',hcpl).fill(cplpro)
		hists.computeIfAbsent('cplkp',hcpl).fill(cplkp)
		hists.computeIfAbsent('cplkm',hcpl).fill(cplkm)
		hists.computeIfAbsent('pt',hpt).fill(pt)
		hists.computeIfAbsent('epkpkmXmm2',hmissmm2).fill(epkpkmX.mass2())
		hists.computeIfAbsent('epkpkmX',hmissmm2).fill(ekpkmX.mass2())
		hists.computeIfAbsent('epkpX',hmissmm2).fill(epkpX.mass2())
		hists.computeIfAbsent('epkmX',hmissmm2).fill(epkmX.mass2())	    	    	   	    

		hists.computeIfAbsent('W',hW).fill(eX.mass())
		hists.computeIfAbsent('missemm2',hmissemm2).fill(epkpkmX.e(), epX.mass2())
		hists.computeIfAbsent('phim',hphim).fill(vphi.mass())
		hists.computeIfAbsent('mm2phim',hmissephim).fill(epX.mass2(), vphi.mass())
		hists.computeIfAbsent('missephim',hmissephim).fill(epkpkmX.e(), vphi.mass())
		hists.computeIfAbsent('mismm2phim',hmissmm2phim).fill(epkpkmX.mass2(), vphi.mass())
		hists.computeIfAbsent('missmm2',hmissmm2).fill(epkpkmX.mass2())
		hists.computeIfAbsent('lambdaphi',hlambdaphi).fill(lambda.mass(),vphi.mass())

		hists.computeIfAbsent('el_ptheta',hptheta).fill(ele.p(), Math.toDegrees(ele.theta()))
		hists.computeIfAbsent('pro_ptheta',hptheta).fill(pro.p(), Math.toDegrees(pro.theta()))
		hists.computeIfAbsent('kp_ptheta',hptheta).fill(kp.p(), Math.toDegrees(kp.theta()))
		hists.computeIfAbsent('km_ptheta',hptheta).fill(km.p(), Math.toDegrees(km.theta()))
		hists.computeIfAbsent('el_phitheta',hphitheta).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta() ) )
		hists.computeIfAbsent('pro_phitheta',hphitheta).fill(Math.toDegrees(pro.phi()), Math.toDegrees(pro.theta()))
		hists.computeIfAbsent('kp_phitheta',hphitheta).fill(Math.toDegrees(kp.phi()), Math.toDegrees(kp.theta()))
		hists.computeIfAbsent('km_phitheta',hphitheta).fill(Math.toDegrees(km.phi()), Math.toDegrees(km.theta()))
		

		return [ele,pro,kp,km,epkpkmX] 
	    }
	    
	    hists.computeIfAbsent('combsize',hcombs).fill(combs.size())
	    if( combs.size() > 0 ){
		//fill the mass histograms by randomly choosing possible events	    
		//def rand_comb=random.nextInt(combs.size())
		def global_event = at.incrementAndGet()
		for( comb in combs){
		    def (ele,pro,kp,km,epkpkmX) = combs[rand_comb]
		    final_el = ele
		    final_pr = pro
		    final_kp = kp
		    final_km = km
		    def hel = event.helicity
		    def rftime = event.rf_time
		    
		    //set event_rec for when simulation is needed
		    event_rec = 1
		    
		    if(data_type == 0 ){
 			fillMCOutTree(event_rec,el_index,el_stat_val,el_lv,pr_index,pr_stat_val,pr_lv,kp_index,kp_stat_val,kp_lv,km_index,km_stat_val,km_lv,event.helicity,global_event,cc)
		    }
		}
	    }
	}
    

	if(data_type == 1){
	    
	    //fillMCOutTree(event_rec,el_index,el_stat_val,el_lv,pr_index,pr_stat_val,pr_lv,kp_index,kp_stat_val,kp_lv,km_index,km_stat_val,km_lv,event.helicity,global_event,cc,final_mcel,final_mcpr,final_mckp,final_mckm,-1000)
	}
	

    }

    def fillElectronHist(ele, event, index){
	def el_sect = event.pcal_sector[index]
	hists.computeIfAbsent('val_el_p',hpele).fill(ele.p())
    	hists.computeIfAbsent('val_el_ptheta',hptheta).fill(ele.p(), Math.toDegrees(ele.theta()))
	hists.computeIfAbsent('val_el_phitheta',hptheta).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta()))

	hists.computeIfAbsent('val_el_p_S'+el_sect,hpele).fill(ele.p())
    	hists.computeIfAbsent('val_el_ptheta_S'+el_sect,hptheta).fill(ele.p(), Math.toDegrees(ele.theta()))
	hists.computeIfAbsent('val_el_phitheta_S'+el_sect,hptheta).fill(Math.toDegrees(ele.phi()), Math.toDegrees(ele.theta()))


    }

    def processEventDetailed(event){
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()

	//get electron proton pair first
	def ep = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    
	    fillElectronHist(ele,event, indx)
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
	    
	    
	    hists.computeIfAbsent('combsize',hcombs).fill(combs.size())
	    if( combs.size() > 0 ){
		def cc=0
		def global_event = at.incrementAndGet()
		for(comb in combs){
		    def (el, pro, kp, km ) = comb		    
		    
		    //write out betas
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


    def processEventDetailedMC(event){

	///////////////////////////////////////////////////////////////
	//set the MC information for each of the final state particles
	def final_mcel = new LorentzVector(0,0,0,0)
	def final_mcpr = new LorentzVector(0,0,0,0)
	def final_mckp = new LorentzVector(0,0,0,0)
	def final_mckm = new LorentzVector(0,0,0,0)
	def event_rec = 0

	//get MC::Particle Variables 
	if(data_type == 1){
	    if( event.mc_status){
		final_mcel = LorentzVector.withPID(11,event.mc_px[0],event.mc_py[0],event.mc_pz[0])
		final_mcpr = LorentzVector.withPID(2212,event.mc_px[1],event.mc_py[1],event.mc_pz[1])
		final_mckp = LorentzVector.withPID(321,event.mc_px[2],event.mc_py[2],event.mc_pz[2])
		final_mckm = LorentzVector.withPID(-321,event.mc_px[3],event.mc_py[3],event.mc_pz[3])
	    }
	}
	//////////////////////////////////////////////////////////////

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
	    
	    
	    hists.computeIfAbsent('combsize',hcombs).fill(combs.size())
	    if( combs.size() > 0 ){
		def cc=0
		def global_event = at.incrementAndGet()
		for(comb in combs){
		    def (el, pro, kp, km ) = comb		    
		    
		    //write out betas
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
		    
		    fillMCOutTree(el_index,el_stat_val,el_lv,pr_index,pr_stat_val,pr_lv,kp_index,kp_stat_val,kp_lv,km_index,km_stat_val,km_lv,event.helicity,global_event,cc,part_betas.get('pro'),part_betas.get('kp'),part_betas.get('km'),final_mcel,final_mcpr,final_mckp,final_mckm)
		    cc+=1
		    
		}
	    }
	    
	    
	}
    }


    def processEventMethod(event, method){
	if( method == 0 ){
	    processEvent(event)
	}
	else if( method == 1){
	    processEventComb(event)
	}
	else if( method == 2){
	    processEventDetailed(event)
	}
	else if( method == 3){
	    //println('here')
	    processEventDetailedMC(event)
	}

    }

}


