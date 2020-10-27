package projects.multiplicity

import org.jlab.io.hipo.HipoDataEvent
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F

import pid.electron.ElectronFromEvent
import utils.PhysicsEvent
import event.Event
import org.jlab.jroot.ROOTFile
import org.jlab.jroot.TNtuple

import java.util.concurrent.ConcurrentHashMap


class PhiMultiplicity{


    def hists = new ConcurrentHashMap()   
    def hmultiplicity = {new H1F("$it","$it",19,1,20)}


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

    def getEventMultiplicity(event){
	
	def e_cand = (0..<event.npart)?.findAll{idx->event.charge[idx]<0}.collect{ii -> [ii, myElectronCutStrategies.collect{ el_cut -> el_cut(event,ii)} ] }.collectEntries()
	
	def ele = e_cand.findResults{cand -> !cand.value.contains(false) ? cand.key : null }
	if(ele.size() == 1){	    	    
	    //def map_multiplicity = [pro:[],kp:[],km:[]]
	    def pr=[]
	    def kp=[]
	    def km=[]
	    //def n_hadrons = (0..<event.npart).findAll{event.pid[it]==2212 || event.pid[it]==321 || event.pid[it]==-321 }	    	   
	    def n_hadrons = (0..<event.npart).each{
		if( event.pid[it]==2212 ){		    
		    //map_multiplicity.get('pro').add(it)
		    pr.add(it)
		}
		else if(event.pid[it]==321){		    
		    //map_multiplicity.get('kp').add(it)
		    kp.add(it)
		}
		else if(event.pid[it]==-321){		    
		    //map_multiplicity.get('km').add(it)
		    km.add(it)
		}
		//println(map_multiplicity)
	    }
	    //look at the total number of hadrons present regardless of conditions on the multiplicity
 	    hists.computeIfAbsent('h_multiplicity_all_hads',hmultiplicity).fill(pr.size() + kp.size() + km.size() )
 
	    //select events with at least one proton + at least one kaonPlus + at least one kaonMinus
	    if(pr.size() >= 1 && kp.size() >= 1 && km.size() >= 1 ){
 		hists.computeIfAbsent('h_multiplicity_integrated',hmultiplicity).fill(pr.size() + kp.size() + km.size() )

 		hists.computeIfAbsent('h_multiplicity_pr_integrated_allkaons',hmultiplicity).fill(pr.size())
 		hists.computeIfAbsent('h_multiplicity_kp_integrated_allprkm',hmultiplicity).fill(kp.size())
 		hists.computeIfAbsent('h_multiplicity_km_integrated_allprkp',hmultiplicity).fill(km.size())						
		
	    }

	    //select events with only proton + at least one kaonPlus + at least one kaonMinus
	    if(pr.size() == 1 && kp.size() >= 1 && km.size() >= 1 ){
 		hists.computeIfAbsent('h_multiplicity_singlepr_integratedkaons',hmultiplicity).fill(pr.size() + kp.size() + km.size() )

 		hists.computeIfAbsent('h_multiplicity_singlepr_pr_integrated_allkaons',hmultiplicity).fill(pr.size())
 		hists.computeIfAbsent('h_multiplicity_singlepr_kp_integrated_allprkm',hmultiplicity).fill(kp.size())
 		hists.computeIfAbsent('h_multiplicity_singlepr_km_integrated_allprkp',hmultiplicity).fill(km.size())						
		
	    }
	    
	    //select events with at least one proton + at least one kaonPlus
	    if(pr.size() >= 1 && kp.size() >= 1){
 		hists.computeIfAbsent('h_multiplicity_prkp',hmultiplicity).fill(pr.size() + kp.size())

 		hists.computeIfAbsent('h_multiplicity_pr_integrated_topprkp',hmultiplicity).fill(pr.size())
 		hists.computeIfAbsent('h_multiplicity_kp_integrated_topprkp',hmultiplicity).fill(kp.size())
	    }

	    //select events with at least one proton + at least one kaonPlus + at least one kaonPlus
	    if(pr.size() >= 1 && km.size() >= 1){
 		hists.computeIfAbsent('h_multiplicity_prkm',hmultiplicity).fill(pr.size() + km.size())

 		hists.computeIfAbsent('h_multiplicity_pr_integrated_topprkm',hmultiplicity).fill(pr.size())
 		hists.computeIfAbsent('h_multiplicity_km_integrated_topprkm',hmultiplicity).fill(km.size())
	    }

	    //select events with at least one kaonPlus + at least one kaonMinus
	    if(kp.size() >= 1 && km.size() >= 1){
		hists.computeIfAbsent('h_multiplicity_kpkm',hmultiplicity).fill(kp.size() + km.size())

 		hists.computeIfAbsent('h_multiplicity_kp_integrated_topkpkm',hmultiplicity).fill(kp.size())
 		hists.computeIfAbsent('h_multiplicity_km_integrated_topkpkm',hmultiplicity).fill(km.size())
	    }

	    
	    
	    
	    

	}

    }

}
