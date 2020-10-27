package pid.kaons

import org.jlab.clas.pdg.PDGDatabase
import org.jlab.detector.base.DetectorType
import event.Event


class KaonPID{

    def ebPID = 0
    def charge = 0

    def kaonPlusPDG = 321
    def kaonMinusPDG = -321

    //dcr1,2,3 fiducial cut parameters    
    def sect_angle_coverage = 60
    // height is lower for outbending runs - this parameter is field config. / sector dependent
    def heightR1 = 16
    def radiusR1 = 25

    def heightR2 = 30
    def radiusR2 = 41

    def heightR3 = 52
    def radiusR3 = 71

    //these need to be stored elsewhere- need to plan out a cut DB procedure
    def kp_mean_p0 = [1.00236, 0.999663, 0.999094, 0.999417, 1.00021, 1.00034]
    def kp_mean_p1 = [0.933781, 0.919325, 0.908676, 0.919357, 0.90765, 0.908293]
    def kp_sig_p0 = [0.0100589, 0.010606, 0.0101618, 0.0113726, 0.0098664, 0.00972719]
    def kp_sig_p1 = [-0.0119585, -0.0125292, -0.0115819, -0.0140509, -0.0117452, -0.0109518]
    def kp_sig_p2 = [0.0208288, 0.020747, 0.0202862, 0.0216464, 0.0207485, 0.0202303]

    def k_mean_p0 = null
    def k_mean_p1 = null
    def k_sig_p0 = null
    def k_sig_p1 = null
    def k_sig_p2 = null

    
    def initCutParameters(charge){
	if( charge > 0 ){
	    ebPID = kaonPlusPDG
	    k_mean_p0 = kp_mean_p0
	    k_mean_p1 = kp_mean_p1
	    k_sig_p0 = kp_mean_p0
	    k_sig_p1 = kp_mean_p1
	    k_sig_p2 = kp_mean_p2
	}
	else if( charnge < 0 ){
	    ebPID = kaonMinusPDG

	}	

    }

    def pr_sig_range = 3.0


    def passProtonEBPIDCut = { event, index ->
	return (event.pid[index] == ebPID)
    }


    // FIDUCIAL CUTS ON DC
    //detector layer r1-12, r2-24, r3-36
    //rotate hit position based on sector
    def rotateDCHitPosition(hit, sec) {
        def ang = Math.toRadians(sec * sect_angle_coverage)
        def x1_rot = hit.get(1) * Math.sin(ang) + hit.get(0) * Math.cos(ang)
        def y1_rot = hit.get(1) * Math.cos(ang) - hit.get(0) * Math.sin(ang)
        return [x1_rot, y1_rot]
    }

    //define left right 
    def borderDCHitPosition(y_rot, height) {
        def slope = 1 / Math.tan(Math.toRadians(0.5 * sect_angle_coverage))
        def left = (height - slope * y_rot)
        def right = (height + slope * y_rot)
        return [left, right]
    }

    def passProtonDCR1 = {event, index ->
	if (event.dc1_status.contains(index)){
	    def sec = event.dc_sector[index]
	    def hit = event.dc1.get(index).find{ hit -> hit.layer == 12}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec-1)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR1)
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1))
            } else {
                return false
            }
	}
	return false
    }
    
    def passProtonDCR2 = { event, index ->
	if (event.dc2_status.contains(index)){
	    def sec = event.dc_sector[index]
	    def hit = event.dc2.get(index).find{ hit -> hit.layer == 24}
            if (hit){
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec-1)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR2)
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1))
            } else {
                return false
            }
	}
	return false
    }

    def passProtonDCR3 = { event, index ->
	if (event.dc3_status.contains(index)){
	    def sec = event.dc_sector[index]
	    def hit = event.dc3.get(index).find{ hit -> hit.layer == 36}
	    if (hit) {
                def hit_rotate = rotateDCHitPosition([hit.x, hit.y], sec-1)
                def left_right = borderDCHitPosition(hit_rotate.get(1), heightR3)

                // from lower line: && x1_rot**2 > radius2_DCr1){ return true }
                return (hit_rotate.get(0) > left_right.get(0) && hit_rotate.get(0) > left_right.get(1))
            } else {
                return false
            }
	}
	return false
    }


    def passProtonBetaFitCut = {event, index ->

	if( event.tof_status.contains(index) ){
	    def ftof_sector = event.tof_sector[index]
	    def particle_beta =  event.beta[index]
	    def p = event.p[index]
	    def mean = prot_mean_p0[sect] * p / Math.sqrt(p**2 + prot_mean_p1[sect]) 
	    def sigma = prot_sig_p0[sect] + prot_sig_p1[sect]/Math.sqrt(p) + prot_sig_p2[sect]/(p**2)
	    // CHANGE 
	    def beta_upper = mean + pr_sig_range*sigma
	    def beta_lower = mean - pr_sig_range*sigma

	    if( particle_beta <= beta_upper && particle_beta >= beta_lower ){
		return true
	    }	    	    
	}
	return false

	//add central beta cut in here to
    }

    
	//PDGDatabase.getParticleMass(2212)**2)


}
