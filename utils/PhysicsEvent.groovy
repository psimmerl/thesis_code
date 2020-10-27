package utils
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.pdg.PDGDatabase

class PhysicsEvent {
    
    def e_p = 'el_indx:el_status:el_px:el_py:el_pz:el_e'
    def pr_p = 'pr_indx:pr_status:pr_px:pr_py:pr_pz:pr_e'
    def kp_p = 'kp_indx:kp_status:kp_px:kp_py:kp_pz:kp_e'
    def km_p = 'km_indx:km_status:km_px:km_py:km_pz:km_e'
    def hel = 'hel'


    //MC::Particle Info
    def event_rec = 'ev_rec'
    def mc_e_p = 'mc_el_px:mc_el_py:mc_el_pz:mc_el_e'
    def mc_pr_p = 'mc_pr_px:mc_pr_py:mc_pr_pz:mc_pr_e'
    def mc_kp_p = 'mc_kp_px:mc_kp_py:mc_kp_pz:mc_kp_e'
    def mc_km_p = 'mc_km_px:mc_km_py:mc_km_pz:mc_km_e'    
    def mc_hel = 'mc_hel'


    def part_kin_string =  null//[e_p,pr_p,kp_p,km_p,hel]   
    def mc_part_kin_string = [mc_e_p,mc_pr_p,mc_kp_p,mc_km_p]//,mc_hel]
    def physics_event=''
    def mc_physics_event=''

    //def analysisTopology = null
    

    def setAnalysisTopology(topology){
	if( topology == 'epkpkm' ){
	    part_kin_string = [e_p,pr_p,kp_p,km_p,hel]   
	}
	else if( topology == 'epkpX' ){
	    part_kin_string = [e_p,pr_p,kp_p,hel]   
	}
	else if( topology == 'epkmX' ){
	    part_kin_string = [e_p,pr_p,km_p,hel]   
	}
    }

    def genTreeString(){
	for( ks in part_kin_string){
	    if( ks != part_kin_string.last() ){
		physics_event+=ks+':'
	    }
	    else if( ks == part_kin_string.last() ){
		physics_event+=ks
	    }
	}
	return physics_event
    }

    def genMCTreeString(){	
	for( ks in mc_part_kin_string){
	    if( ks != mc_part_kin_string.last() ){
		mc_physics_event+=ks+':'
	    }
	    else if( ks == mc_part_kin_string.last() ){
		mc_physics_event+=ks
	    }
	}
	return mc_physics_event
    }


}
