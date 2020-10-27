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

class epkpkm {
  
    def phyev = new PhysicsEvent()        
    def f1 = null 
    def t1 = null
    
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
    def hW = {new H1F("$it","$it",100,1,4)}
    def hmissmm2phim = {new H2F("$it","$it",100,-1,1,100,0,5)}
    def hmissmm2 = {new H1F("$it","$it",500,-1,1.1)}
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

    def electron = new ElectronFromEvent();

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

    def createROOTFile(fout){
	f1 = new ROOTFile(fout)
	println(' DATA TYPE OF OUTPUT FILE IS ' + data_type )
	if( this.data_type == 1 ){
	     t1 = f1.makeNtuple('clas12','clas12phi',phyev.genTreeString()+":"+phyev.genMCTreeString()+":global_event"+":comb_index:"+"prb:"+"kpb:"+"kmb")	    
	}
	else{
	    t1 = f1.makeNtuple('clas12','clas12phi',phyev.genTreeString())
	}

    }
    
    def fillOutTree(int el_index, int el_status, LorentzVector el, int pr_index, int pr_status, LorentzVector pr, int kp_index, int kp_status, LorentzVector kp, int km_index, int km_status, LorentzVector km, int hel){
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
		hel
	)

    }

    def fillMCOutTree(int event_rec, int el_status, LorentzVector el, int pr_index, int pr_status, LorentzVector pr, int kp_index, int kp_status, LorentzVector kp, int km_index, int km_status, LorentzVector km, int hel, int glob_ev, int comb_index, LorentzVector mcel, LorentzVector mcpr, LorentzVector mckp, LorentzVector mckm, int mchel){
 	t1.fill(event_rec,
		el_index,
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
		mckm.e(),
		mchel
	)

    }

    def closeOutTree(){
	t1.write()
	f1.close()
    }

    def processEvent(event) {	
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()
//	println('here ')
	def combs = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
 	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    	    
	    //println(' got electron - finding proton next ')
	    return ele
	}.collectMany{ ele -> 
	    (0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		//println(' found proton - setting lv ' )
		def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 	
		//println(pro.e())
		//return (pro.p()>0.5 && pro.p()<4) ? [ele, pro] : null
		return [ele,pro]
		//return pro
            }//.take(1).collect{pro -> [ele,pro] }
	}.collectMany{ele,pro->            
	    (0..<event.npart).findAll{event.pid[it]==321}.findResults{ikp->
		//println(' found kaon - setting lv ' )
		//println(' el e ' + ele.e() + ' pro e ' + pro.e() )
		def kp = LorentzVector.withPID(321,event.px[ikp], event.py[ikp], event.pz[ikp]) 		
		//return (kp.p()<2.5 && Math.toDegrees(kp.theta())<35) ? [ele,pro,kp] : null
		return [ele,pro,kp]
            }
        }.collectMany{ele,pro,kp->  
	    (0..<event.npart).findAll{event.pid[it]==-321}.findResults{ikm->
		def km = LorentzVector.withPID(-321,event.px[ikm], event.py[ikm], event.pz[ikm]) 						
		//return (km.p()<2.5 && Math.toDegrees(km.theta())<40) ? [ele,pro,kp,km] : null
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
	    
	    //hists.computeIfAbsent('cplpro',hcpl).fill(cplpro)
	    //hists.computeIfAbsent('cplkp',hcpl).fill(cplkp)
	    //hists.computeIfAbsent('cplkm',hcpl).fill(cplkm)
	    //hists.computeIfAbsent('pt',hpt).fill(pt)
	    //hists.computeIfAbsent('epkpkmXmm2',hmissmm2).fill(epkpkmX.mass2())
	    //hists.computeIfAbsent('epkpkmX',hmissmm2).fill(ekpkmX.mass2())
	    //hists.computeIfAbsent('epkpX',hmissmm2).fill(epkpX.mass2())
	    //hists.computeIfAbsent('epkmX',hmissmm2).fill(epkmX.mass2())	    	    	   	    
	    return [ele,pro,kp,km,epkpkmX] 
	}
		
	hists.computeIfAbsent('combsize',hcombs).fill(combs.size())
	//next we sort the possible combinations based on missing mass2 and return combination with smallest value
	def min_mm2_event = combs.sort{ it.get(4).mass2().abs() }.take(1).collect{min_mm2_event->; return min_mm2_event}
 
	//if( combs.size() > 0 ){
	    //println(' final event is this ' )
	    //println(' number of combinations is ' + combs.size() )
	    //println( combs.each{ it.each{  println(it.e() ) } } )
	    //println(' getting event with small missing mass ' )
	    //println(min_mm2_event)
	    

	    //fill the mass histograms by randomly choosing possible events
	    
	  //  def rand_comb=random.nextInt(combs.size())
	   // def (ele,pro,kp,km,epkpkmX) = combs[rand_comb]
	   // def hel = event.helicity
	    //def rftime = event.rf_time
	    
	   /// fillOutTree(ele,pro,kp,km,hel)
	    //println('combo size is ' + combs.size() + ' picking random number ' + rand_event)
	    
	//}
	
	if(min_mm2_event.size()==1) {            
	    def (ele,pro,kp,km,epkpkmX) = min_mm2_event[0]
	    def hel = event.helicity
	    def rftime = event.rf_time
	    
	    //fillOutTree(ele,pro,kp,km,hel)
	}	    
	

	    /*
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
	    def e1 = Math.pow(y*xb*target.mass(),2)/q2
	    def eps = (1-y-e1)/(1-y+ y*y/2 + e1)
	    def w2 = target.mass()*target.mass() + 2*target.mass()*nu - q2
	    def w = Math.sqrt(w2)
	    def wth = target.mass() + phim 
	    def e1cm = target.mass()*(target.mass() + nu)/w
	    def p1cm = target.mass()*Math.sqrt(nu**2 + q2)/w
	    def e2cm = (w2 + target.mass()**2 - phim*phim)/(2.0*w)
	    def p2cm = Math.sqrt(e2cm*e2cm - target.mass()**2)
	    def e3cm = (w2 - phim**2 + target.mass()**2)/(2*w)
	    def p3cm = Math.sqrt(e3cm**2 - target.mass()**2)
	    def tmin = 2*(target.mass()**2 - e1cm*e2cm - p1cm*p2cm)
	    def tmax = 2*(target.mass()**2 - e1cm*e2cm + p1cm*p2cm)
	    
	    def boost_phicm = vphi.boostVector()
	    boost_phicm.negative()
	    def boost_kp = new LorentzVector(kp)
	    def boost_km = new LorentzVector(km)
	    boost_kp.boost(boost_phicm)//.x(),boost_phicm.y(),boost_phicm.z())
	    boost_km.boost(boost_phicm)
	    def cos_theta_cm_kp = Math.cos(boost_kp.vect().theta())

	    def v3l = beam.vect().cross(ele.vect())
	    def v3h = pro.vect().cross((beam-ele).vect())
	    def trento = Math.toDegrees(Math.acos(v3l.dot(v3h)/(v3l.mag()*v3h.mag())))

	    if(v3l.dot(pro.vect())<0){trento=360-trento}

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
	    hists.computeIfAbsent('q2x',hq2x).fill(xb,q2)
	    hists.computeIfAbsent('txb',htxb).fill(t,xb)
	    hists.computeIfAbsent('tq2',hq2t).fill(t,q2)
	    hists.computeIfAbsent('xb',hx).fill(xb)
	    hists.computeIfAbsent('t',ht).fill(t)
	    hists.computeIfAbsent('q2',hq).fill(q2)
	    hists.computeIfAbsent('cos_kp',hcos).fill(cos_theta_cm_kp)

	    //event based information
	    //def hel = banks.ev.getByte('helicity',0)
	    //def rftime = banks.ev.getFloat('RFTime',0)
	    def hel = event.helicity
	    def rftime = event.rf_time

	    hists.computeIfAbsent('helicity',hhel).fill(hel)
	    hists.computeIfAbsent('rftime',htime).fill(rftime)
	    
	    if( hel > 0 ){hists.computeIfAbsent('pos_asy',hasy).fill(trento)}
	    else if( hel < 0 ){hists.computeIfAbsent('neg_asy',hasy).fill(trento)}

	    if( vphi.mass().abs() < 1.05 ){
		if( hel > 0 ){hists.computeIfAbsent('pos_asy_phicut',hasy).fill(trento)}
		else if( hel < 0 ){hists.computeIfAbsent('neg_asy_phicut',hasy).fill(trento)}
	    }

	    
            hists.computeIfAbsent('W',hW).fill(eX.mass())
            hists.computeIfAbsent('missemm2',hmissemm2).fill(epkpkmX.e(), epX.mass2())
            hists.computeIfAbsent('phim',hphim).fill(vphi.mass())
            hists.computeIfAbsent('mm2phim',hmissephim).fill(epX.mass2(), vphi.mass())
            hists.computeIfAbsent('missephim',hmissephim).fill(epkpkmX.e(), vphi.mass())
            hists.computeIfAbsent('mismm2phim',hmissmm2phim).fill(epkpkmX.mass2(), vphi.mass())
            hists.computeIfAbsent('missmm2',hmissmm2).fill(epkpkmX.mass2())

     }
	 */
    }



    def processEventMax(event){
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()

	def combMax = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }.findResults{ indx -> 
	    def ele = LorentzVector.withPID(11, event.px[indx], event.py[indx], event.pz[indx]) 	    
	    return ele
	}.collectMany{ ele -> 
	    (0..<event.npart).findAll{event.pid[it]==2212}.findResults{ipro->
		def pro = LorentzVector.withPID(2212,event.px[ipro], event.py[ipro], event.pz[ipro]) 
		return pro
            }.sort{-it.e()}.take(1).collect{ pro -> [ele,pro] } //  return proton with sort 
	}.collectMany{ele,pro->
	    (0..<event.npart).findAll{event.pid[it]==321}.findResults{ikp->
		def kp = LorentzVector.withPID(321,event.px[ikp], event.py[ikp], event.pz[ikp]) 		
		return kp
            }sort{-it.e()}.take(1).collect{ kp -> [ele,pro,kp] }
        }.collectMany{ele,pro,kp->  
	    (0..<event.npart).findAll{event.pid[it]==-321}.findResults{ikm->
		def km = LorentzVector.withPID(-321,event.px[ikm], event.py[ikm], event.pz[ikm]) 				
		return km
            }.sort{-it.e()}.take(1).collect{ km-> [ele,pro,kp,km] } // pick last element (-1) . take[-1] also works 
        }.findResults{ele,pro,kp,km->	    
	    def epkpkmX = beam+target-ele-pro-kp-km
	    
	    def ekpkmX = beam+target-ele-kp-km
	    def epkpX = beam+target-ele-pro-kp
	    def epkmX = beam+target-ele-pro-km
	    
	    def cplpro = ekpkmX.vect().theta(pro.vect())
	    def cplkm = epkpX.vect().theta(km.vect())
	    def cplkp = epkmX.vect().theta(kp.vect())
	    def pt = Math.sqrt(epkpkmX.px()**2 + epkpkmX.py()**2)

	    //hists.computeIfAbsent('cplpro',hcpl).fill(cplpro)
	    //hists.computeIfAbsent('cplkp',hcpl).fill(cplkp)
	    //hists.computeIfAbsent('cplkm',hcpl).fill(cplkm)
	    //hists.computeIfAbsent('pt',hpt).fill(pt)
	    //hists.computeIfAbsent('epkpkmXmm2',hmissmm2).fill(epkpkmX.mass2())
	    //hists.computeIfAbsent('ekpkmXmm2',hmissmm2).fill(ekpkmX.mass2())
	    //hists.computeIfAbsent('epkpXmm2',hmissmm2).fill(epkpX.mass2())
	    //hists.computeIfAbsent('epkmXmm2',hmissmm2).fill(epkmX.mass2())

	    //return (epkpkmX.mass2().abs()<0.06 && cplpro<30 && cplkm<20 && cplkp<20 && pt<0.5) ? [ele,pro,kp,km,epkpkmX,ekpkmX,epkpX,epkmX] : null
	    return [ele,pro,kp,km]
	}

	hists.computeIfAbsent('combsize',hcombs).fill(combMax.size())
	if(combMax.size()==1) {
	    def (ele,pro,kp,km) = combMax[0]
	    def hel = event.helicity
	    def rftime = event.rf_time

	    //fillOutTree(ele,pro,kp,km,hel)
	}
	/*
	if(combMax.size()==1) {
            def (ele,pro,kp,km,epkpkmX,ekpkmX,epkpX,epkmX) = combMax[0]
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
	    def e1 = Math.pow(y*xb*target.mass(),2)/q2
	    def eps = (1-y-e1)/(1-y+ y*y/2 + e1)
	    def w2 = target.mass()*target.mass() + 2*target.mass()*nu - q2
	    def w = Math.sqrt(w2)
	    def wth = target.mass() + phim 
	    def e1cm = target.mass()*(target.mass() + nu)/w
	    def p1cm = target.mass()*Math.sqrt(nu**2 + q2)/w
	    def e2cm = (w2 + target.mass()**2 - phim*phim)/(2.0*w)
	    def p2cm = Math.sqrt(e2cm*e2cm - target.mass()**2)
	    def e3cm = (w2 - phim**2 + target.mass()**2)/(2*w)
	    def p3cm = Math.sqrt(e3cm**2 - target.mass()**2)
	    def tmin = 2*(target.mass()**2 - e1cm*e2cm - p1cm*p2cm)
	    def tmax = 2*(target.mass()**2 - e1cm*e2cm + p1cm*p2cm)
	    
	    def boost_phicm = vphi.boostVector()
	    boost_phicm.negative()
	    def boost_kp = new LorentzVector(kp)
	    def boost_km = new LorentzVector(km)
	    boost_kp.boost(boost_phicm)//.x(),boost_phicm.y(),boost_phicm.z())
	    boost_km.boost(boost_phicm)
	    def cos_theta_cm_kp = Math.cos(boost_kp.vect().theta())

	    def v3l = beam.vect().cross(ele.vect())
	    def v3h = pro.vect().cross((beam-ele).vect())
	    def trento = Math.toDegrees(Math.acos(v3l.dot(v3h)/(v3l.mag()*v3h.mag())))

	    if(v3l.dot(pro.vect())<0){trento=360-trento}

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
	    hists.computeIfAbsent('q2x',hq2x).fill(xb,q2)
	    hists.computeIfAbsent('txb',htxb).fill(t,xb)
	    hists.computeIfAbsent('tq2',hq2t).fill(t,q2)
	    hists.computeIfAbsent('xb',hx).fill(xb)
	    hists.computeIfAbsent('t',ht).fill(t)
	    hists.computeIfAbsent('q2',hq).fill(q2)
	    hists.computeIfAbsent('cos_kp',hcos).fill(cos_theta_cm_kp)

	    //event based information
	    def hel = event.helicity
	    def rftime = event.rf_time

	    hists.computeIfAbsent('helicity',hhel).fill(hel)
	    hists.computeIfAbsent('rftime',htime).fill(rftime)
	   
	    if( hel > 0 ){hists.computeIfAbsent('pos_asy',hasy).fill(trento)}
	    else if( hel < 0 ){hists.computeIfAbsent('neg_asy',hasy).fill(trento)}

	    if( vphi.mass().abs() < 1.05 ){
		if( hel > 0 ){hists.computeIfAbsent('pos_asy_phicut',hasy).fill(trento)}
		else if( hel < 0 ){hists.computeIfAbsent('neg_asy_phicut',hasy).fill(trento)}
	    }


	}
	 */
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
   
    def processEventRandomComb(event) {		
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
			//fillOutTree(ele,pro,kp,km,hel)	    
		    }
		}
	    }
	}
    

	if(data_type == 1){
	    
	    //	    fillMCOutTree(event_rec,el_index,el_stat_val,el_lv,pr_index,pr_stat_val,pr_lv,kp_index,kp_stat_val,kp_lv,km_index,km_stat_val,km_lv,event.helicity,global_event,cc,final_mcel,final_mcpr,final_mckp,final_mckm,-1000)
	}
	

    }



    def processEventMethod(event, method){
	if( method == 0 ){
	    processEvent(event)
	}
	else if( method == 1){
	    processEventMax(event)
	}
	else if( method == 2){
	    processEventRandomComb(event)
	}
    }
    
}


