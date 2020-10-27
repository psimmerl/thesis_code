package utils
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.pdg.PDGDatabase


class KinTool{

    //move somplace else
    static def p0mean = [0.105631, 0.11551, 0.112799, 0.109937, 0.116249, 0.119057]
    static  def p1mean = [-0.153951, -0.0253273, -0.125718, 0.165414, 0.0768411, 0.0555026]
    static  def p2mean = [0.00860091, 0.00706291, 0.00908884, 0.00499666, 0.00448701, 0.00558927]
    static  def p3mean = [-0.000821675, -0.000711488, -0.000930922, -0.000298311, -0.000455716, -0.000657084]
    static  def p0sigma = [0.0149613, 0.0115116, 0.00580737, 0.0106817, 0.012667, 0.00553471]
    static def  p1sigma = [0.00700773, 0.0116193, 0.0202375, 0.0126958, 0.00892239, 0.0216206]

    

    static def Q2(){
	return 0
    }

    static def delta_meas_energy( Double beam,  LorentzVector measured_el ){
	def calc_e = beam/(1+ (beam/PDGDatabase.getParticleMass(2212))*(1-Math.cos(measured_el.theta())) );
	return (calc_e  - measured_el.e() )
    }

    static def delta_meas_theta( Double beam,  LorentzVector measured_el ){
	def calc_theta = Math.toDegrees(Math.acos( 1 + (PDGDatabase.getParticleMass(2212)/beam)*( 1 - beam/measured_el.e()) ))
	return Math.toDegrees((calc_theta  - Math.toDegrees(measured_el.theta()) ))
    }

    static def getSamplingFractionFromMomentum(Integer sector, Double p, Double sf_meas){
	def sf_calc = p0mean[sector] * (1 + p/Math.sqrt(p*p + p1mean[sector])) + p2mean[sector]*p + p3mean[sector]*p*p
	def sf_sigma = p0sigma[sector] + p1sigma[sector]/Math.sqrt(p)
	def sf_chi = (sf_meas - sf_calc)/sf_sigma
	return sf_chi
    }

}


