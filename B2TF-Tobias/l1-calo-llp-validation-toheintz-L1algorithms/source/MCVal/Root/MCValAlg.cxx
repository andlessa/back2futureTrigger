#include <AsgMessaging/MessageCheck.h>
#include <MCVal/MCValAlg.h>
#include <TruthUtils/HepMCHelpers.h>
#include <bits/stdc++.h>

MCValAlg :: MCValAlg (const std::string& name,
		ISvcLocator *pSvcLocator)
	: EL::AnaAlgorithm (name, pSvcLocator)
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  This is also where you
	// declare all properties for your algorithm.  Note that things like
	// resetting statistics variables or booking histograms should
	// rather go into the initialize() function.

	declareProperty( "pdgIdBSM", m_pdgIdBSM = 55,  "BSM Higgs/Scalar PDGId" );
	declareProperty( "pdgIdHS", m_pdgIdHS = 4000023,  "Chi1 PDGId" );
	declareProperty( "pdgIdchi0", m_pdgIdchi0 = 4000022,  "Chi0 PDGId" );

}

std::vector<float> extrapolate (const xAOD::TruthParticle* part) {

	// This function returns the position of a truth particle within the calorimeters
	// N.B. It is intended to be used for detector-stabel truth particles
	// It returns a vector with (x,y,y,t,readout time). All coordinates have units of mm
	// The exact return values depend on the production vertex position:
	// a) Particles with production vertex within the ECAL or HCAL: x,y,z,t are just the coordinates of the production vertes. The readout time is calculated by subtracting r=sqrt(x*x+y*y+z*z) from t.
	// b) Particles with production vertex within the Tracker: Starting from the production vertex (x0,y0,z0,t0) and knowing the momentum (and energy) of the particle, a linear extrapolation to the calorimeter entry (i.e. the inner boundary of the ECAL) is done to get the (x,y,z,t) of the particle within the camorimeter. The readout time is given by by correcting t by the "travelled distance", i.e. the distance from the IA to the production vertex and the distance from the production vertex to the calorimeter entry
	// c) Particles produced outside the HCAL are currently not considered, therefore the time and readout time are set to 1e10 //FIXME need to think about the case that the momentum could bring the particle back into the CALO [should be rare]
	// d) particles without production vertex are not considered
	
	//Inner boundaries of the ECAL
	float rhoCalo = 1400;
        float zCalo = 3700;

	//Outer boundaries of the HCAL
	float rhoCaloMax = 3500;
	float zCaloMax = 5500;

	float x0 = 0 , y0 = 0 , z0 = 0 , t0 = 1e10;
	//Get production vertex
        if (part->hasProdVtx()){
		x0 = part->prodVtx()->x(); //mm
        	y0 = part->prodVtx()->y(); //mm
        	z0 = part->prodVtx()->z(); //mm
		t0 = part->prodVtx()->t(); //mm
	}
	float r0 = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0); //mm
	float readout0 = t0 - r0; //mm

	// case d) returns (0,0,0,1e10,1e10), case a) returns (x0,y0,z0,t0,t0-r0)
	std::vector<float> result {x0 , y0 , z0 , t0 , readout0};

	// case c) returns (x0,y0,y0,1e10,1e10) //We might want to return (0,0,0,1e10,1e10) instead?
	if (part->hasProdVtx() and (part->prodVtx()->perp() > rhoCaloMax or part->prodVtx()->z() > zCaloMax)){
		result.at(3) = 1e10;  
		result.at(4) = 1e10; 
	}

	// In case b) we extrapolate the calorimeter entry:
        if (part->hasProdVtx() and (part->prodVtx()->perp() < rhoCalo and part->prodVtx()->z() < zCalo)){
		//N.B. we're extrapolating the production vertex linearly to the CALO entry (assuming that the particles are energetic enough such that they are moving on straight paths rather than circular motion)
	        float px = part->px();
	        float py = part->py();
	        float pz = part->pz();
		float e = part->e();
		float vx = px/e;
		float vy = py/e;
		float vz = pz/e;

		//approximate x(t) = x0 + vx * t etc. (N.B. for simplicity we start t=0 at the production vertex, rather than accounting for the time offset t0)
		//solve rho(t) = sqrt(x(t)^2 + y(t)^2) = rhoCalo and z(t) = zCalo
		
		//vector to store the solutions of the equations
                std::vector<float> t;
		//(i) for rho, we solve quadratic equation t_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a} with:
		float a = vx * vx  + vy * vy;
	        float b = 2 * (vx * x0 + vy * y0);
        	float c = x0 * x0 + y0 * y0 - rhoCalo * rhoCalo;

        	if (b * b >= 4 * a * c and a != 0) {
        		float discriminant = std::sqrt(b * b - 4 * a * c);
	                float t1 = (-b + discriminant) / (2 * a);
        	        float t2 = (-b - discriminant) / (2 * a);
                	if (t1 >= 0) t.push_back(t1); //we're only interested in positive time solutions
	                if (t2 >= 0) t.push_back(t2); //we're only interested in positive time solutions
        	}

		//(ii) in z direction the equation is linear
		if (vz != 0){
			float tZ = (zCalo - z0) / vz;
			if (tZ >= 0) t.push_back(tZ); //we're only interested in positive time solutions
		}

		//take the smallest solution as the time of the calorimeter entry. Calculate the correspionding x,y,z values, and the readout time
		if (!t.empty()) {
			float tEntry = *std::min_element(t.begin(), t.end()); //take the smallest solution
			float xEntry = vx * tEntry;
                        float yEntry = vy * tEntry;
                        float zEntry = vz * tEntry;
                        float rEntry = std::sqrt(xEntry * xEntry + yEntry * yEntry + zEntry * zEntry);

                        result.at(0) = x0 + xEntry;
                        result.at(1) = y0 + yEntry;
                        result.at(2) = z0 + zEntry;
                        result.at(3) = t0 + tEntry;
                        result.at(4) = (t0 + tEntry) - (r0 + rEntry);
		}
		else { //extrapolation failed, e.g. vx = vy = vz = 0
			result.at(3) = 1e10;
                        result.at(4) = 1e10;
		}
	} //end of case b)
	// case b) returns either (xEntry,yEntry,zEntry,tEntry,tEmtryReadout) or in rare occasions (x0,y0,y0,1e10,1e10) //We might want to return (0,0,0,1e10,1e10) instead?
	
	return result;
}

float remainingInteractionLength(float eta0 , float rho0 , float z0){

	//This is a function to caclculate the remaining interaction length from a jet originating in the ECAL from its "production vertex" at (x0,y0,z0) to the ECAL exit

	//std::abs(eta0) is binned from 0 to 3.3 in 0.1 steps, i.e. first bin from 0 to 0.1 and last bin from 3.2 to 3.3
	int etaBin = int(10 * std::abs(eta0));

	//Material budget in interaction length (values taken from figure 5.2 in https://iopscience.iop.org/article/10.1088/1748-0221/3/08/S08003):
	const float LAMBDA_EXIT[33] = {2.25 , 2.25 , 2.20 , 2.30 , 2.35 , 2.50 , 2.60 , 2.85 , 2.92 , 3.25 , 3.70 , 3.95 , 3.95 , 3.85 , 4.20 , 2.95 , 2.78 , 2.50 , 2.35 , 2.20 , 2.30 , 2.30 , 2.30 , 2.30 , 2.40 , 2.10 , 2.10 , 2.20 , 2.15 , 2.20 , 2.20 , 2.30 , 1.80};
	const float LAMBDA_ENTRY[33] = {0.60 , 0.65 , 0.55 , 0.60 , 0.62 , 0.65 , 0.70 , 0.80 , 0.85 , 0.90 , 0.98 , 1.00 , 1.10 , 1.25 , 1.40 , 1.55 , 1.52 , 1.10 , 0.82 , 0.68 , 0.72 , 0.70 , 0.65 , 0.60 , 0.60 , 0.65 , 0.65 , 0.70 , 0.65 , 0.65  , 0.60 , 0.55 , 0.50};

	float lambda_exit = LAMBDA_EXIT[etaBin];
	float lambda_entry = LAMBDA_ENTRY[etaBin];

	//We'll scale the lambda difference between the ECAL entry and exit according to the ratio of the ECAL thickness and the "remaing" distance from the production vertex to the ECAL exit
	//in the Barrel region:
	float rho_exit = 2000; 
	float rho_entry = 1400; 
	float ratio_Barrel = (rho_exit - rho0) / (rho_exit - rho_entry); //ratio remaining distance : ECAL thickness

        //in the forward region:
        float z_exit = 4300; 
        float z_entry = 3700; 
	float ratio_Endcap = (z_exit - z0) / (z_exit - z_entry); //ratio remaining distance : ECAL thickness

	//in the transition region:
	//the jet is travelling through the Barrel as well as the endcap ECAL
	//for the distance (and thickness) calculations, we need to account for the Barrel (as above) and the endcap (as above). 
	//Based on eta, we linear combine both distances:
	float alpha = (std::abs(eta0) - 1.37) / (1.70 - 1.37); //or 1.375 to 1.475 instead? In our "simplified" geometry, the transition region goes from 1.37 to 1.70, in the actual geometry from 1.375 to 1.475
	//Due to the linear combination, the sin(theta) and cos(theta) factors do not cancel out in the ratio
	float theta = 2 * std::atan(std::exp(-eta0));
	float ratio_Transition = ((1 - alpha) * ((rho_exit - rho0) / std::sin(theta)) + alpha * ((z_exit - z0) / std::cos(theta))) / ((1 - alpha) * ((rho_exit - rho_entry) / std::sin(theta)) + alpha * ((z_exit - z_entry) / std::cos(theta)));

	//based on eta, we grab the according ratio
	float ratio = 0;
	if (std::abs(eta0) < 1.37)		ratio = ratio_Barrel; //or 1.375 instead?
	else if (std::abs(eta0) > 1.70) 	ratio = ratio_Endcap; //or 1.475 instead? 
	else 					ratio = ratio_Transition;
	
	//We'll scale the lambda difference according to the ratio
	return (lambda_exit - lambda_entry) * ratio;
}

StatusCode MCValAlg :: initialize ()
{

	ATH_MSG_ALWAYS("THE PDGID OF THE SCALAR IS " << m_pdgIdBSM);
	ATH_MSG_ALWAYS("THE PDGID OF THE CHI1 IS " << m_pdgIdHS);
	ATH_MSG_ALWAYS("THE PDGID OF THE CHI0 IS " << m_pdgIdchi0);
	//N.B. the PDGIDs are actually set by https://gitlab.cern.ch/toheintz/l1-calo-llp-validation/-/blob/toheintz/L1algorithms/source/MCVal/share/MCVal_eljob.py?ref_type=heads#L48-50
	
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.

	// Create output tree
	ANA_CHECK (book (TTree ("analysis", "My analysis ntuple")));
	TTree* mytree = tree ("analysis");

	// Declare branches
	mytree->Branch ("RunNumber", &m_runNumber);
	mytree->Branch ("EventNumber", &m_eventNumber);

	//L1Calo MET in BC N-1
        mytree->Branch ("L1_MET_Nminus1_met", &m_L1_MET_Nminus1_met);
        mytree->Branch ("L1_MET_Nminus1_mpt", &m_L1_MET_Nminus1_mpt);
        mytree->Branch ("L1_MET_Nminus1_phi", &m_L1_MET_Nminus1_phi);
        mytree->Branch ("L1_MET_Nminus1_eta", &m_L1_MET_Nminus1_eta);
        mytree->Branch ("L1_MET_Nminus1_px", &m_L1_MET_Nminus1_px);
        mytree->Branch ("L1_MET_Nminus1_py", &m_L1_MET_Nminus1_py);
        mytree->Branch ("L1_MET_Nminus1_pz", &m_L1_MET_Nminus1_pz);

	//Truth MET in BC N-1
	mytree->Branch ("truth_MET_Nminus1_mpt", &m_truth_MET_Nminus1_mpt);
	mytree->Branch ("truth_MET_Nminus1_met", &m_truth_MET_Nminus1_met);

        //jFex smallR jets
        m_L1_jet_N_e = new std::vector<float>();
        m_L1_jet_N_et = new std::vector<float>();
        m_L1_jet_N_eta = new std::vector<float>();
        m_L1_jet_N_phi = new std::vector<float>();

	m_L1_delPhi = new std::vector<float>();

	mytree->Branch ("L1_jet_N_e", &m_L1_jet_N_e);
	mytree->Branch ("L1_jet_N_et", &m_L1_jet_N_et);
        mytree->Branch ("L1_jet_N_eta", &m_L1_jet_N_eta);
        mytree->Branch ("L1_jet_N_phi", &m_L1_jet_N_phi);
        mytree->Branch ("L1_jet_N_N", &m_L1_jet_N_N);
	
        mytree->Branch ("L1_delPhi", &m_L1_delPhi);

        mytree->Branch ("L1_leadingJet_N_e", &m_L1_leadingJet_N_e);
        mytree->Branch ("L1_leadingJet_N_et", &m_L1_leadingJet_N_et);
        mytree->Branch ("L1_leadingJet_N_eta", &m_L1_leadingJet_N_eta);
        mytree->Branch ("L1_leadingJet_N_phi", &m_L1_leadingJet_N_phi);
        mytree->Branch ("L1_leading_delPhi", &m_L1_leading_delPhi);

	//offline jets (HLT)
        m_offline_jet_N_e = new std::vector<float>();
        m_offline_jet_N_et = new std::vector<float>();
        m_offline_jet_N_eta = new std::vector<float>();
        m_offline_jet_N_phi = new std::vector<float>();

        mytree->Branch ("offline_jet_N_e", &m_offline_jet_N_e);
        mytree->Branch ("offline_jet_N_et", &m_offline_jet_N_et);
        mytree->Branch ("offline_jet_N_eta", &m_offline_jet_N_eta);
        mytree->Branch ("offline_jet_N_phi", &m_offline_jet_N_phi);
        mytree->Branch ("offline_jet_N_N", &m_offline_jet_N_N);

	//Scalar Particle
	m_scalar_pt = new std::vector<float>();
	m_scalar_eta = new std::vector<float>();
	m_scalar_phi = new std::vector<float>();
	m_scalar_m = new std::vector<float>();
	m_scalar_pdgid = new std::vector<int>();
	m_scalar_status = new std::vector<int>();

	mytree->Branch ("scalar_pt", &m_scalar_pt);
	mytree->Branch ("scalar_eta", &m_scalar_eta);
	mytree->Branch ("scalar_phi", &m_scalar_phi);
	mytree->Branch ("scalar_m", &m_scalar_m);
	mytree->Branch ("scalar_pdgid", &m_scalar_pdgid);
	mytree->Branch ("scalar_N", &m_scalar_N);
	mytree->Branch ("scalar_status", &m_scalar_status);

	//reconstruct scalar particle from its decay products:
	mytree->Branch ("scalarReco_pt", &m_scalarReco_pt);
	mytree->Branch ("scalarReco_pz", &m_scalarReco_pz);
	mytree->Branch ("scalarReco_m", &m_scalarReco_m);

	//Higgs Boson (only on-shell higgs bosons are stored in this vector)
	m_higgs_pt = new std::vector<float>();
	m_higgs_eta = new std::vector<float>();
	m_higgs_phi = new std::vector<float>();
	m_higgs_m = new std::vector<float>();
	m_higgs_pdgid = new std::vector<int>();
	m_higgs_status = new std::vector<int>();

	mytree->Branch ("higgs_pt", &m_higgs_pt);
	mytree->Branch ("higgs_eta", &m_higgs_eta);
	mytree->Branch ("higgs_phi", &m_higgs_phi);
	mytree->Branch ("higgs_m", &m_higgs_m);
	mytree->Branch ("higgs_pdgid", &m_higgs_pdgid);
	mytree->Branch ("higgs_N", &m_higgs_N);
	mytree->Branch ("higgs_status", &m_higgs_status);

	//reconstruct Higgs bosons from their decay products (b quarks)
	//on Time Higgs (reconstructed from b quarks)
	mytree->Branch ("onTimeHiggsReco_et", &m_onTimeHiggsReco_et);
	mytree->Branch ("onTimeHiggsReco_pt", &m_onTimeHiggsReco_pt);
	mytree->Branch ("onTimeHiggsReco_pz", &m_onTimeHiggsReco_pz);
	mytree->Branch ("onTimeHiggsReco_eta", &m_onTimeHiggsReco_eta);
	mytree->Branch ("onTimeHiggsReco_phi", &m_onTimeHiggsReco_phi);
	mytree->Branch ("onTimeHiggsReco_m", &m_onTimeHiggsReco_m);
	mytree->Branch ("onTimeHiggsReco_ProdVtx_t", &m_onTimeHiggsReco_ProdVtx_t);
	mytree->Branch ("onTimeHiggsReco_ProdVtx_Lxy", &m_onTimeHiggsReco_ProdVtx_Lxy);
	mytree->Branch ("onTimeHiggsReco_ProdVtx_Lz", &m_onTimeHiggsReco_ProdVtx_Lz);

	//out of time Higgs (reconstructed from b quarks)
	mytree->Branch ("outOfTimeHiggsReco_et", &m_outOfTimeHiggsReco_et);
	mytree->Branch ("outOfTimeHiggsReco_pt", &m_outOfTimeHiggsReco_pt);
	mytree->Branch ("outOfTimeHiggsReco_pz", &m_outOfTimeHiggsReco_pz);
	mytree->Branch ("outOfTimeHiggsReco_eta", &m_outOfTimeHiggsReco_eta);
	mytree->Branch ("outOfTimeHiggsReco_phi", &m_outOfTimeHiggsReco_phi);
	mytree->Branch ("outOfTimeHiggsReco_m", &m_outOfTimeHiggsReco_m);
	mytree->Branch ("outOfTimeHiggsReco_ProdVtx_t", &m_outOfTimeHiggsReco_ProdVtx_t);
	mytree->Branch ("outOfTimeHiggsReco_ProdVtx_Lxy", &m_outOfTimeHiggsReco_ProdVtx_Lxy);
	mytree->Branch ("outOfTimeHiggsReco_ProdVtx_Lz", &m_outOfTimeHiggsReco_ProdVtx_Lz);

	//delta phi between reconstructed Higgs	
	mytree->Branch ("HiggsReco_del_phi" , &m_HiggsReco_del_phi);	
	
	//clustered Higgs bosons from b jets //FIXME previously used for jet and MET estimation. Only used for comaprison currently.
        //on Time Higgs (reclustered from b jets)
	mytree->Branch ("onTimeHiggsPseudoApproach2_et" , &m_onTimeHiggsPseudoApproach2_et);
	mytree->Branch ("onTimeHiggsPseudoApproach2_pt" , &m_onTimeHiggsPseudoApproach2_pt);
	mytree->Branch ("onTimeHiggsPseudoApproach2_pz" , &m_onTimeHiggsPseudoApproach2_pz);	
	mytree->Branch ("onTimeHiggsPseudoApproach2_eta" , &m_onTimeHiggsPseudoApproach2_eta);
	mytree->Branch ("onTimeHiggsPseudoApproach2_phi" , &m_onTimeHiggsPseudoApproach2_phi);
	mytree->Branch ("onTimeHiggsPseudoApproach2_m" , &m_onTimeHiggsPseudoApproach2_m);
	mytree->Branch ("onTimeHiggsPseudoApproach2_ProdVtx_t" , &m_onTimeHiggsPseudoApproach2_ProdVtx_t);
	mytree->Branch ("onTimeHiggsPseudoApproach2_ProdVtx_Lxy" , &m_onTimeHiggsPseudoApproach2_ProdVtx_Lxy);
	mytree->Branch ("onTimeHiggsPseudoApproach2_ProdVtx_Lz" , &m_onTimeHiggsPseudoApproach2_ProdVtx_Lz);

	//out of Time Higgs (reclustered from b jets)
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_et" , &m_outOfTimeHiggsPseudoApproach2_et);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_pt" , &m_outOfTimeHiggsPseudoApproach2_pt);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_pz" , &m_outOfTimeHiggsPseudoApproach2_pz);	
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_eta" , &m_outOfTimeHiggsPseudoApproach2_eta);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_phi" , &m_outOfTimeHiggsPseudoApproach2_phi);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_m" , &m_outOfTimeHiggsPseudoApproach2_m);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_ProdVtx_t" , &m_outOfTimeHiggsPseudoApproach2_ProdVtx_t);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_ProdVtx_Lxy" , &m_outOfTimeHiggsPseudoApproach2_ProdVtx_Lxy);
	mytree->Branch ("outOfTimeHiggsPseudoApproach2_ProdVtx_Lz" , &m_outOfTimeHiggsPseudoApproach2_ProdVtx_Lz);
        
	//delta phi between clustered Higgs
	mytree->Branch ("HiggsPseudoApproach2_del_phi" , &m_HiggsPseudoApproach2_del_phi);	

	//BSM Chi1 particle
	m_bsm_pt = new std::vector<float>();
	m_bsm_eta = new std::vector<float>();
	m_bsm_phi = new std::vector<float>();
	m_bsm_m = new std::vector<float>();
	m_bsm_pdgid = new std::vector<int>();
	m_bsm_status = new std::vector<int>();

	mytree->Branch ("bsm_pt", &m_bsm_pt);
	mytree->Branch ("bsm_eta", &m_bsm_eta);
	mytree->Branch ("bsm_phi", &m_bsm_phi);
	mytree->Branch ("bsm_m", &m_bsm_m);
	mytree->Branch ("bsm_pdgid", &m_bsm_pdgid);
	mytree->Branch ("bsm_N", &m_bsm_N);
	mytree->Branch ("bsm_status", &m_bsm_status);

	mytree->Branch("llp_del_phi",&m_llp_del_phi);
	mytree->Branch("llp_del_phi_eff",&m_llp_del_phi_eff);
	mytree->Branch("llp_del_R",&m_llp_del_R);

	//out of time Chi1
	mytree->Branch("nonpromptLLPpt",&m_llp1_pt);
	mytree->Branch("nonpromptLLPphi",&m_llp1_phi);
	mytree->Branch("nonpromptLLPpx",&m_llp1_px);
	mytree->Branch("nonpromptLLPpy",&m_llp1_py);
	mytree->Branch("nonpromptLLPe",&m_llp1_e);
	mytree->Branch("nonpromptLLPpz",&m_llp1_pz);
	mytree->Branch("nonpromptLLPLxy",&m_llp1_Lxy);
	mytree->Branch("nonpromptLLPLz",&m_llp1_Lz);
	mytree->Branch("nonpromptLLPt",&m_llp1_t);
	mytree->Branch("nonpromptLLPt_eff",&m_llp1_t_eff);
	mytree->Branch("nonpromptLLPctau",&m_llp1_ctau);
	mytree->Branch("nonpromptLLPboost",&m_llp1_boost);

	//on-time Chi1
	mytree->Branch("promptLLPpt",&m_llp2_pt);
	mytree->Branch("promptLLPe",&m_llp2_e);
	mytree->Branch("promptLLPphi",&m_llp2_phi);
	mytree->Branch("promptLLPpx",&m_llp2_px);
	mytree->Branch("promptLLPpy",&m_llp2_py);
	mytree->Branch("promptLLPLxy",&m_llp2_Lxy);
	mytree->Branch("promptLLPLz",&m_llp2_Lz);
	mytree->Branch("promptLLPt",&m_llp2_t);
	mytree->Branch("promptLLPt_eff",&m_llp2_t_eff);
	mytree->Branch("promptLLPctau",&m_llp2_ctau);
	mytree->Branch("promptLLPboost",&m_llp2_boost);

	//all chi0
	m_chi0_status = new std::vector<int>();
	mytree->Branch ("chi0_status", &m_chi0_status);

	//chi1-matched chi0
	m_chi0_llp_status = new std::vector<int>();
	mytree->Branch ("chi0_llp_status", &m_chi0_llp_status);
	mytree->Branch("chi0_llp_N", &m_chi0_llp_N); //Number of chi1->gamma+chi0 vertices

	mytree->Branch("chi0_llp_del_phi",&m_chi0_llp_del_phi);
	mytree->Branch("chi0_llp_del_R",&m_chi0_llp_del_R);
	mytree->Branch("chi0_llp_del_t",&m_chi0_llp_del_t);

	//delayed-chi1 matched chi0
	mytree->Branch("chi0_llp1_pt",&m_chi0_llp1_pt);
	mytree->Branch("chi0_llp1_e",&m_chi0_llp1_e);
	mytree->Branch("chi0_llp1_px",&m_chi0_llp1_px);
	mytree->Branch("chi0_llp1_py",&m_chi0_llp1_py);
	mytree->Branch("chi0_llp1_eta",&m_chi0_llp1_eta);
	mytree->Branch("chi0_llp1_phi", &m_chi0_llp1_phi);

	//on-time-chi1-matched chi0
	mytree->Branch("chi0_llp2_pt", &m_chi0_llp2_pt);
	mytree->Branch("chi0_llp2_e", &m_chi0_llp2_e);
	mytree->Branch("chi0_llp2_px", &m_chi0_llp2_px);
	mytree->Branch("chi0_llp2_py", &m_chi0_llp2_py);
	mytree->Branch("chi0_llp2_pz", &m_chi0_llp2_pz);
	mytree->Branch("chi0_llp2_eta", &m_chi0_llp2_eta);
	mytree->Branch("chi0_llp2_phi", &m_chi0_llp2_phi);
	
        //neutrinos (all flavours)
        m_nu_pt = new std::vector<float>();
        m_nu_e = new std::vector<float>();
        m_nu_px = new std::vector<float>();
        m_nu_py = new std::vector<float>();
        m_nu_eta = new std::vector<float>();
        m_nu_phi = new std::vector<float>();
        m_nu_m = new std::vector<float>();
        m_nu_pdgid = new std::vector<int>();
        m_nu_status = new std::vector<int>();

        mytree->Branch ("nu_pt", &m_nu_pt);
        mytree->Branch ("nu_e", &m_nu_e);
        mytree->Branch ("nu_px", &m_nu_px);
        mytree->Branch ("nu_py", &m_nu_py);
        mytree->Branch ("nu_eta", &m_nu_eta);
        mytree->Branch ("nu_phi", &m_nu_phi);
        mytree->Branch ("nu_m", &m_nu_m);
        mytree->Branch ("nu_pdgid", &m_nu_pdgid);
        mytree->Branch ("nu_N", &m_nu_N);
        mytree->Branch ("nu_status", &m_nu_status);

	//leptons (all flavours)
	m_lep_pt = new std::vector<float>();
	m_lep_eta = new std::vector<float>();
	m_lep_phi = new std::vector<float>();
	m_lep_m = new std::vector<float>();
	m_lep_pdgid = new std::vector<int>();
	m_lep_status = new std::vector<int>();

	mytree->Branch ("lep_pt", &m_lep_pt);
	mytree->Branch ("lep_eta", &m_lep_eta);
	mytree->Branch ("lep_phi", &m_lep_phi);
	mytree->Branch ("lep_m", &m_lep_m);
	mytree->Branch ("lep_pdgid", &m_lep_pdgid);
	mytree->Branch ("lep_N", &m_lep_N);
	mytree->Branch ("lep_status", &m_lep_status);

	//electrons
	m_el_pt = new std::vector<float>();
	m_el_eta = new std::vector<float>();
	m_el_phi = new std::vector<float>();
	m_el_m = new std::vector<float>();
	m_el_status = new std::vector<int>();

	mytree->Branch ("el_pt", &m_el_pt);
	mytree->Branch ("el_eta", &m_el_eta);
	mytree->Branch ("el_phi", &m_el_phi);
	mytree->Branch ("el_m", &m_el_m);
	mytree->Branch ("el_N", &m_el_N);
	mytree->Branch ("el_status", &m_el_status);

	//muons
	m_mu_pt = new std::vector<float>();
	m_mu_eta = new std::vector<float>();
	m_mu_phi = new std::vector<float>();
	m_mu_m = new std::vector<float>();
	m_mu_status = new std::vector<int>();

	mytree->Branch ("mu_pt", &m_mu_pt);
	mytree->Branch ("mu_eta", &m_mu_eta);
	mytree->Branch ("mu_phi", &m_mu_phi);
	mytree->Branch ("mu_m", &m_mu_m);
	mytree->Branch ("mu_N", &m_mu_N);
	mytree->Branch ("mu_status", &m_mu_status);

	//quarks (all flavours)
	m_q_pt = new std::vector<float>();
	m_q_eta = new std::vector<float>();
	m_q_phi = new std::vector<float>();
	m_q_m = new std::vector<float>();
	m_q_pdgid = new std::vector<int>();
	m_q_status = new std::vector<int>();

	mytree->Branch ("q_pt", &m_q_pt);
	mytree->Branch ("q_eta", &m_q_eta);
	mytree->Branch ("q_phi", &m_q_phi);
	mytree->Branch ("q_m", &m_q_m);
	mytree->Branch ("q_pdgid", &m_q_pdgid);
	mytree->Branch ("q_N", &m_q_N);
	mytree->Branch ("q_status", &m_q_status);

	//b quarks
	m_bq_charge = new std::vector<double>();
	m_bq_pt = new std::vector<float>();
	m_bq_eta = new std::vector<float>();
	m_bq_phi = new std::vector<float>();
	m_bq_m = new std::vector<float>();
    	m_bq_status = new std::vector<int>();

	m_bq_prodVtx_perp = new std::vector<float>;
	m_bq_prodVtx_z = new std::vector<float>;
	m_bq_prodVtx_t = new std::vector<float>;
	m_bq_prodVtx_eta = new std::vector<float>;
	m_bq_prodVtx_phi = new std::vector<float>;
	m_bq_prodVtx_nIn = new std::vector<int>;
	m_bq_prodVtx_nOut = new std::vector<int>;
	m_bq_prodVtx_ID	= new std::vector<int>;
	m_bq_prodVtx_barcode = new std::vector<int>;

	mytree->Branch ("bq_charge", &m_bq_charge);
	mytree->Branch ("bq_pt", &m_bq_pt);
	mytree->Branch ("bq_eta", &m_bq_eta);
	mytree->Branch ("bq_phi", &m_bq_phi);
	mytree->Branch ("bq_m", &m_bq_m);
    	mytree->Branch ("bq_status", &m_bq_status);
	mytree->Branch ("bq_N", &m_bq_N);

	mytree->Branch ("bq_prodVtx_perp", &m_bq_prodVtx_perp);
	mytree->Branch ("bq_prodVtx_z", &m_bq_prodVtx_z);
	mytree->Branch ("bq_prodVtx_t", &m_bq_prodVtx_t);
	mytree->Branch ("bq_prodVtx_eta", &m_bq_prodVtx_eta);
	mytree->Branch ("bq_prodVtx_phi", &m_bq_prodVtx_phi);
	mytree->Branch ("bq_prodVtx_nIn", &m_bq_prodVtx_nIn);
	mytree->Branch ("bq_prodVtx_nOut", &m_bq_prodVtx_nOut);
	mytree->Branch ("bq_prodVtx_ID", &m_bq_prodVtx_ID);
	mytree->Branch ("bq_prodVtx_barcode", &m_bq_prodVtx_barcode);

	//jets (all flavours)
        m_j_pt = new std::vector<float>();
        m_j_eta = new std::vector<float>();
        m_j_phi = new std::vector<float>();
        m_j_m = new std::vector<float>();
        m_j_status = new std::vector<int>();
        m_j_truthcode = new std::vector<int>();

        mytree->Branch ("j_pt", &m_j_pt);
        mytree->Branch ("j_eta", &m_j_eta);
        mytree->Branch ("j_phi", &m_j_phi);
        mytree->Branch ("j_m", &m_j_m);
        mytree->Branch ("j_status", &m_j_status);
        mytree->Branch ("j_N", &m_j_N);
        mytree->Branch ("j_truthcode", &m_j_truthcode);

	//b jets
	m_bj_pt = new std::vector<float>();
	m_bj_eta = new std::vector<float>();
	m_bj_phi = new std::vector<float>();
	m_bj_m = new std::vector<float>();
	m_bj_status = new std::vector<int>();
	m_bj_truthcode = new std::vector<int>();

	mytree->Branch ("bj_pt", &m_bj_pt);
	mytree->Branch ("bj_eta", &m_bj_eta);
	mytree->Branch ("bj_phi", &m_bj_phi);
	mytree->Branch ("bj_m", &m_bj_m);
	mytree->Branch ("bj_status", &m_bj_status);
	mytree->Branch ("bj_N", &m_bj_N);
	mytree->Branch ("bj_truthcode", &m_bj_truthcode);

	return StatusCode::SUCCESS;
}

//function to sort an std::vector<const xAOD::TruthParticle*> based on the particles' pT
bool comparePt(const xAOD::TruthParticle* part1, const xAOD::TruthParticle* part2)
{
	return (part1->pt()>part2->pt());
}

//function to sort an std::vector<const xAOD::TruthParticle*> based on the particles' decay time
bool compareTdecay(const xAOD::TruthParticle* part1, const xAOD::TruthParticle* part2)
{
	if(part1->hasDecayVtx() && part2-> hasDecayVtx())		return (part1->decayVtx()->t()>part2->decayVtx()->t());
	else return true;
}

//function to sort an std::vector<const xAOD::TruthParticle*> based on the particles' production time
bool compareTprod(const xAOD::TruthParticle* part1, const xAOD::TruthParticle* part2)
{
	if(part1->hasProdVtx() && part2-> hasProdVtx())		return (part1->prodVtx()->t()>part2->prodVtx()->t());
	else return true;
}


StatusCode MCValAlg :: execute ()
{
	
	// Here you do everything that needs to be done on every single
	// events, e.g. read input variables, apply cuts, and fill
	// histograms and trees.  This is where most of your actual analysis
	// code will go.

	// Truth particles
	const xAOD::TruthParticleContainer* truth_particles = nullptr;
	ANA_CHECK ( evtStore()->retrieve(truth_particles, "TruthParticles" ) );
	
	// Truth Jets
	const xAOD::JetContainer* truth_jets = nullptr;
	ANA_CHECK ( evtStore()->retrieve(truth_jets, "AntiKt4TruthDressedWZJets" ) );
	//we rather want to have "AntiKt4TruthJets" to also include prompt leptons from Higgs decays
	//see e.g. https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TruthDAOD
	//but it seems they are not used anymore (and it is recommended to use DressedWZJets instead, e.g. https://atlas-talk.web.cern.ch/t/non-wz-antik4-truth-jet-collection-in-daod-phys-for-rel22/30731/4)
	//ANA_CHECK ( evtStore()->retrieve(truth_jets, "AntiKt4TruthJets" ) ); //returns error, as AntiKt4TruthJets are not stored in TRUTH0, TRUTH1, and TRUTH3
	const xAOD::JetContainer* truth_largeR_jets = nullptr;
	ANA_CHECK ( evtStore()->retrieve(truth_largeR_jets, "AntiKt10TruthTrimmedPtFrac5SmallR20Jets" ) );
	
	//truth vertices >> currently not used, though might be useful
	const xAOD::TruthVertexContainer *truth_vertices = nullptr;
	ANA_CHECK ( evtStore()->retrieve(truth_vertices, "TruthVertices") );

	//create vectors to store the different particle types
	std::vector<const xAOD::TruthParticle*> vec_finalStateParticles;
	std::vector<const xAOD::TruthParticle*> vec_scalar;
	std::vector<const xAOD::TruthParticle*> vec_higgs;
	std::vector<const xAOD::TruthParticle*> vec_llp;
	std::vector<const xAOD::TruthParticle*> vec_chi0;
	std::vector<const xAOD::TruthParticle*> vec_chi0_llp;
	std::vector<const xAOD::TruthParticle*> vec_nu;
	std::vector<const xAOD::TruthParticle*> vec_lep;
	std::vector<const xAOD::TruthParticle*> vec_el;
	std::vector<const xAOD::TruthParticle*> vec_mu;
	std::vector<const xAOD::TruthParticle*> vec_q;
	std::vector<const xAOD::TruthParticle*> vec_bquark;
	
	//create vectors to store the different jet types
	std::vector<const xAOD::Jet_v1*> vec_jet;
	std::vector<const xAOD::Jet_v1*> vec_bjet;
	std::vector<const xAOD::Jet_v1*> vec_largeRjet;

	//create vectors to store different vertex types
	std::vector<const xAOD::TruthVertex_v1*> vec_vertices;

	//Fill the jet vectors:
	//https://atlas-sw-doxygen.web.cern.ch/atlas-sw-doxygen/atlas_main--Doxygen/docs/html/db/d8c/classxAOD_1_1Jet__v1.html
	for(auto jet : *truth_jets){
		vec_jet.push_back(jet);
		int truthcode;
		// Default is actually here 
		// https://ftag.docs.cern.ch/algorithms/labelling/jet_labels/#hadron-based-labelling
		jet->getAttribute("HadronConeExclTruthLabelID",truthcode); 	
		if(abs(truthcode) == 5){
			vec_bjet.push_back(jet);
		}
	}	
	
	for(auto jet : *truth_largeR_jets){
		vec_largeRjet.push_back(jet);
	}	

	//Fill the particle vectors:
	int idx = 0;
	for(auto tp : *truth_particles){
		idx +=1;

 		if (std::abs(tp->status()) == 1){
                	vec_finalStateParticles.push_back(tp);
                }

		// Skip particles that come from hadron decays
		if(!Truth::notFromHadron(tp)) continue;

		// Skip gluons
		if(MC::isGluon(tp->pdgId())) continue;

		//scalar particle
		if(std::abs(tp->pdgId()) == m_pdgIdBSM) {
			if(tp->parent()->pdgId() == tp->pdgId()) continue;
			vec_scalar.push_back(tp);
		}

		//Higgs bosons
		if(tp->isHiggs()) {
			vec_higgs.push_back(tp);
		}

		// BSM particles
		if(std::abs(tp->pdgId()) == m_pdgIdHS) {
			// Skip if it is a self-decay
			if(tp->parent()->pdgId() == tp->pdgId()) continue;
			// Store BSM particle
			vec_llp.push_back(tp);
		}

		/*
		if(Truth::isFromParticle(tp,m_pdgIdHS)) {
			vec_child.push_back(tp);
		}
		*/

		//Pick any chi0, then LLP chi0
		if(std::abs(tp->pdgId()) == m_pdgIdchi0){
			vec_chi0.push_back(tp);
			if(std::abs(tp->parent()->pdgId())==m_pdgIdHS){	
				vec_chi0_llp.push_back(tp);
			}
		}

		//neutrinos
		if(std::abs(tp->pdgId()) == 12 or std::abs(tp->pdgId()) == 14 or std::abs(tp->pdgId()) == 16){
			if (std::abs(tp->status()) == 1)	vec_nu.push_back(tp);
		}

		//Electrons
		if(Truth::isFinalElectron(tp)) {
			vec_el.push_back(tp);
			vec_lep.push_back(tp);
		}

		//Muons
		if(Truth::isFinalMuon(tp)) {
			vec_mu.push_back(tp);
			vec_lep.push_back(tp);
		}

		//Quarks
		if(Truth::isFinalQuark(tp) && !Truth::isFromGluon(tp)) {
			vec_q.push_back(tp);
			if(std::abs(tp->pdgId())==5)		vec_bquark.push_back(tp);
		}

	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//sort all vectors which are related to on-time and delayed chi1 according to their time, first entry is largest time (i.e. most delayed particle)
	sort(vec_llp.begin(),vec_llp.end(),compareTdecay); //based on their decay vertex
	sort(vec_chi0_llp.begin(),vec_chi0_llp.end(),compareTprod); //based on their production vertex

	// Clear scalar particle vectors
	m_scalar_pt->clear();
	m_scalar_eta->clear();
	m_scalar_phi->clear();
	m_scalar_m->clear();
	m_scalar_pdgid->clear();
	m_scalar_status->clear();

	// Store number of scalar particles
	m_scalar_N = vec_scalar.size();

	// Loop over scalar particles
	for(auto scalar : vec_scalar) {
		// Fill vectors for branches
		m_scalar_pt->push_back(scalar->pt());
		m_scalar_eta->push_back(scalar->eta());
		m_scalar_phi->push_back(scalar->phi());
		m_scalar_m->push_back(scalar->m());
		m_scalar_pdgid->push_back(scalar->pdgId());		
		m_scalar_status->push_back(scalar->status());
	}

	// Clear BSM vectors (chi1)
	m_bsm_pt->clear();
	m_bsm_eta->clear();
	m_bsm_phi->clear();
	m_bsm_m->clear();
	m_bsm_pdgid->clear();
	m_bsm_status->clear();

	// Store number of BSM particles
	m_bsm_N = vec_llp.size();

	// Loop over BSM particles
	for(auto bsm : vec_llp) {
		// Fill vectors for branches
		m_bsm_pt->push_back(bsm->pt());
		m_bsm_eta->push_back(bsm->eta());
		m_bsm_phi->push_back(bsm->phi());
		m_bsm_m->push_back(bsm->m());
		m_bsm_pdgid->push_back(bsm->pdgId());		
		m_bsm_status->push_back(bsm->status());
	}

	//Clear llp vectors
	m_llp1_pt=-999;
	m_llp1_phi=-999;
	m_llp1_e=-999;
	m_llp1_px=-999;
	m_llp1_py=-999;
	m_llp1_pz=-999;
	m_llp1_Lxy=-999;
	m_llp1_Lz=-999999;
	m_llp1_t=-999;
	m_llp1_t_eff=-999;
	m_llp1_boost=-999;
	m_llp1_ctau=-999;

	m_llp2_pt=-999;
	m_llp2_e=-999;
	m_llp2_phi=-999;
	m_llp2_px=-999;
	m_llp2_py=-999;
	m_llp2_Lxy=-999;
	m_llp2_Lz=-999999;
	m_llp2_t=-999;
	m_llp2_t_eff=-999;
	m_llp2_ctau=-999;
	m_llp2_boost=-999;

	m_llp_del_phi=-999;
	m_llp_del_phi_eff=-999;
	m_llp_del_R=-999;

	const xAOD::TruthParticle* nonpromptLLP = nullptr;
	const xAOD::TruthParticle* promptLLP = nullptr;
	int i = 0;
	double nonpromptLLPt = -1;

	if (vec_llp.size() != 2)		ATH_MSG_WARNING("This event contains " << vec_llp.size() << " BSM particles. This is NOT expected");
	for(auto llp : vec_llp){
		if(i==0){
			nonpromptLLP = llp;
			m_llp1_pt = llp->pt()/Truth::GeV;
			m_llp1_phi = llp->phi();
			m_llp1_e = llp->e();
			m_llp1_px = llp->px();
			m_llp1_py = llp->py();
			m_llp1_pz = llp->pz();
			//FIXME change unit of m_llp1_pt and m_llp1_t
			if(llp->hasDecayVtx())
			{
				const xAOD::TruthVertex_v1* decVtx = llp->decayVtx();
				nonpromptLLPt = decVtx->t();
				m_llp1_Lxy = decVtx->perp();
				m_llp1_Lz = decVtx->z();
				m_llp1_t = decVtx->t()/3e2;// Convert to ns
				double gamma = llp->e()/llp->m();
				m_llp1_boost = pow(1-pow(gamma,-2),0.5);
				//FIXME add member variable m_llp1_gamma
				double r = sqrt(pow(decVtx->perp(),2)+pow(decVtx->z(),2));
				m_llp1_ctau = r/gamma;
			}
			else		ATH_MSG_WARNING("This event contains a BSM particle without decay vertex. This is NOT expected");
		}
		if(i==1){
			promptLLP = llp;
			m_llp2_pt = llp->pt()/Truth::GeV;
			m_llp2_e = llp->e();
			m_llp2_phi = llp->phi();
			m_llp2_px = llp->px();
			m_llp2_py = llp->py();
			//FIXME change unit of m_llp2_pt and m_llp2_t
			if(llp->hasDecayVtx())
			{
				const xAOD::TruthVertex_v1* decVtx = llp->decayVtx();
				m_llp2_Lxy = decVtx->perp();
				m_llp2_Lz = decVtx->z();
				m_llp2_t = decVtx->t()/3e2;// Convert to ns
				double gamma = llp->e()/llp->m();
				m_llp2_boost = pow(1-pow(gamma,-2),0.5);
				//FIXME add member variable m_llp2_gamma
				//m_pho_llp_del_t = std::abs((decVtx->t()-nonpromptLLPt)/3e2); //Convert from mm/c to ns
				m_chi0_llp_del_t = std::abs((decVtx->t()-nonpromptLLPt)/3e2); //Convert from mm/c to ns
				double r = sqrt(pow(decVtx->perp(),2)+pow(decVtx->z(),2));
				m_llp2_ctau = r/gamma;
			}
			else            ATH_MSG_WARNING("This event contains a BSM particle without decay vertex. This is NOT expected");
			//hist("h_llp_decayDiffRatio")->Fill((m_llp2_t-m_llp1_t)/(m_llp2_t+m_llp1_t));
			//FIXME add member variable m_llp_decayDiffRatio
		}
		i++;
		
	} 

	if(promptLLP != 0 && nonpromptLLP != 0 && nonpromptLLP->hasDecayVtx() && promptLLP->hasDecayVtx()){

                auto llpv1 = nonpromptLLP->decayVtx()->v4();
                auto llpv2 = promptLLP->decayVtx()->v4();
                m_llp_del_phi = std::abs(llpv1.DeltaPhi(llpv2));
                m_llp_del_R = std::abs(llpv1.DeltaR(llpv2));

        }

	TLorentzVector scalarReco;
	m_scalarReco_pt = -999.;
	m_scalarReco_pz = -999.;
	m_scalarReco_m = -999.;
	if (promptLLP != 0 and nonpromptLLP != 0){
		scalarReco = promptLLP->p4()+nonpromptLLP->p4();
		m_scalarReco_pt = scalarReco.Pt();
		m_scalarReco_pz = scalarReco.Pz();
		m_scalarReco_m = scalarReco.M();
	}
	else			ATH_MSG_WARNING("This event contains a nullptr as BSM particle.");


	// Clear LLP chi0 vector
	m_chi0_llp_N=-999;
	m_chi0_llp1_pt =-999;
	m_chi0_llp1_e =-999;
	m_chi0_llp1_px =-999;
	m_chi0_llp1_py =-999;
	m_chi0_llp1_eta =-999;
	m_chi0_llp1_phi =-999;
	m_chi0_llp2_pt =-999;
	m_chi0_llp2_e =-999;
	m_chi0_llp2_px =-999;
	m_chi0_llp2_py =-999;
	m_chi0_llp2_pz =-999;
	m_chi0_llp2_eta =-999;
	m_chi0_llp2_phi =-999;
	m_chi0_llp_del_phi =-999;
	m_chi0_llp_del_R =-999;
	m_chi0_status->clear();
	m_chi0_llp_status->clear();

	// Loop over chi0s
	for(auto chi0 : vec_chi0) {
		m_chi0_status->push_back(chi0->status());
	}

	// Loop over truth matched chi0s
	for(auto chi0_llp : vec_chi0_llp) {
		m_chi0_llp_status->push_back(chi0_llp->status());
	}

	m_chi0_llp_N = vec_chi0_llp.size();
	if(vec_chi0_llp.size() == 2){ 					//change from "==2" to ">=2" makes this more inclusive
		auto promptChi0 =  vec_chi0_llp.at(0);
		auto notpromptChi0 = vec_chi0_llp.at(1);
		m_chi0_llp1_pt = promptChi0->pt();
		m_chi0_llp1_e = promptChi0->e();
		m_chi0_llp1_px = promptChi0->px();
		m_chi0_llp1_py = promptChi0->py();
		m_chi0_llp1_eta = promptChi0->eta();
		m_chi0_llp1_phi = promptChi0->phi();
		m_chi0_llp2_pt = notpromptChi0->pt();
		m_chi0_llp2_e = notpromptChi0->e();
		m_chi0_llp2_px = notpromptChi0->px();
		m_chi0_llp2_py = notpromptChi0->py();
		m_chi0_llp2_pz = notpromptChi0->pz();
		m_chi0_llp2_eta = notpromptChi0->eta();
		m_chi0_llp2_phi = notpromptChi0->phi();
		m_chi0_llp_del_phi = std::abs(promptChi0->p4().DeltaPhi(notpromptChi0->p4()));
                m_chi0_llp_del_phi = std::abs(promptChi0->p4().DrEtaPhi(notpromptChi0->p4()));
	}
	else			ATH_MSG_WARNING("This event contains " << vec_chi0_llp.size() << " BSM Chi0 particles. This is NOT expected.");


	// Clear neutrino vectors
        m_nu_pt->clear();
        m_nu_e->clear();
        m_nu_px->clear();
        m_nu_py->clear();
        m_nu_eta->clear();
        m_nu_phi->clear();
        m_nu_m->clear();
        m_nu_pdgid->clear();
        m_nu_status->clear();

        // Store number of neutrinos
        m_nu_N = vec_nu.size();

        // Loop over neutrinos
        for(auto nu : vec_nu) {
                // Fill vectors for branches
                m_nu_pt->push_back(nu->pt());
                m_nu_px->push_back(nu->px());
                m_nu_py->push_back(nu->py());
                m_nu_e->push_back(nu->e());
                m_nu_eta->push_back(nu->eta());
                m_nu_phi->push_back(nu->phi());
                m_nu_m->push_back(nu->m());
                m_nu_pdgid->push_back(nu->pdgId());
                m_nu_status->push_back(nu->status());
        }

	// Clear lepton vectors
	m_lep_pt->clear();
	m_lep_eta->clear();
	m_lep_phi->clear();
	m_lep_m->clear();
	m_lep_pdgid->clear();
	m_lep_status->clear();

	// Store number of leptons
	m_lep_N = vec_lep.size();

	// Loop over leptons
	for(auto lep : vec_lep) {
		// Fill vectors for branches
		m_lep_pt->push_back(lep->pt());
		m_lep_eta->push_back(lep->eta());
		m_lep_phi->push_back(lep->phi());
		m_lep_m->push_back(lep->m());
		m_lep_pdgid->push_back(lep->pdgId());
		m_lep_status->push_back(lep->status());
	}


	// Clear electron vectors
	m_el_pt->clear();
	m_el_eta->clear();
	m_el_phi->clear();
	m_el_m->clear();
	m_el_status->clear();

	// Store number of electrons
	m_el_N = vec_el.size();

	// Loop over electrons
	for(auto el : vec_el) {
		// Fill vectors for branches
		m_el_pt->push_back(el->pt());
		m_el_eta->push_back(el->eta());
		m_el_phi->push_back(el->phi());
		m_el_m->push_back(el->m());		
		m_el_status->push_back(el->status());
	}
	
	// Clear muon vectors
	m_mu_pt->clear();
	m_mu_eta->clear();
	m_mu_phi->clear();
	m_mu_m->clear();
	m_mu_status->clear();

	// Store number of muons
	m_mu_N = vec_mu.size();

	// Loop over muons
	for(auto mu : vec_mu) {
		// Fill vectors for branches
		m_mu_pt->push_back(mu->pt());
		m_mu_eta->push_back(mu->eta());
		m_mu_phi->push_back(mu->phi());
		m_mu_m->push_back(mu->m());
		m_mu_status->push_back(mu->status());
	}


	// Clear quark vectors
	m_q_pt->clear();
	m_q_eta->clear();
	m_q_phi->clear();
	m_q_m->clear();
	m_q_pdgid->clear();
	m_q_status->clear();

	// Store number of quarks
	m_q_N = vec_q.size();
	// Loop over quarks
	for(auto q : vec_q) {
		// Fill vectors for branches
		m_q_pt->push_back(q->pt());
		m_q_eta->push_back(q->eta());
		m_q_phi->push_back(q->phi());
		m_q_m->push_back(q->m());
		m_q_pdgid->push_back(q->pdgId());
		m_q_status->push_back(q->status());
	}


	// Clear b quark vectors
	m_bq_charge->clear();
	m_bq_pt->clear();
	m_bq_eta->clear();
	m_bq_phi->clear();
	m_bq_m->clear();
	m_bq_status->clear();
	m_bq_N = vec_bquark.size();

	m_bq_prodVtx_perp->clear();
	m_bq_prodVtx_z->clear();
	m_bq_prodVtx_t->clear();
	m_bq_prodVtx_eta->clear();
	m_bq_prodVtx_phi->clear();
	m_bq_prodVtx_nIn->clear();
	m_bq_prodVtx_nOut->clear();
	m_bq_prodVtx_ID->clear();
	m_bq_prodVtx_barcode->clear();

	//loop over bquarks
	for(auto quark : vec_bquark){

		m_bq_charge->push_back(quark->charge());
		m_bq_pt->push_back(quark->pt());
		m_bq_eta->push_back(quark->eta());
		m_bq_phi->push_back(quark->phi());
		m_bq_m->push_back(quark->m());
		m_bq_status->push_back(quark->status());
		
		m_bq_prodVtx_perp -> push_back(quark->prodVtx()->perp());
		m_bq_prodVtx_z -> push_back(quark->prodVtx()->z());
		m_bq_prodVtx_t -> push_back(quark->prodVtx()->t());
		m_bq_prodVtx_eta -> push_back(quark->prodVtx()->eta());
		m_bq_prodVtx_phi -> push_back(quark->prodVtx()->phi());
		m_bq_prodVtx_nIn -> push_back(quark->prodVtx()->nIncomingParticles());
		m_bq_prodVtx_nOut -> push_back(quark->prodVtx()->nOutgoingParticles());
		m_bq_prodVtx_ID -> push_back(quark->prodVtx()->id());
		m_bq_prodVtx_barcode -> push_back(quark->prodVtx()->barcode());
	}


	//Clear jet vectors (all flavours)
	m_j_pt->clear();
        m_j_eta->clear();
        m_j_phi->clear();
        m_j_m->clear();
        m_j_status->clear();
        m_j_truthcode->clear();
        m_j_N = vec_jet.size();

        for(auto jet : vec_jet){
                int truthcode;
                jet->getAttribute("HadronConeExclTruthLabelID",truthcode);
                m_j_truthcode->push_back(truthcode);
                m_j_pt->push_back(jet->pt());
                m_j_eta->push_back(jet->eta());
                m_j_phi->push_back(jet->phi());
                m_j_m->push_back(jet->m());
                //m_j_status->push_back(-1);
        }


	// Clear b jet vectors
	m_bj_pt->clear();
	m_bj_eta->clear();
	m_bj_phi->clear();
	m_bj_m->clear();
	m_bj_status->clear();
	m_bj_truthcode->clear();
	m_bj_N = vec_bjet.size();

	//for(auto jet : vec_bjet){
	for(auto jet : vec_bjet){
		int truthcode;
		jet->getAttribute("HadronConeExclTruthLabelID",truthcode); 
		m_bj_truthcode->push_back(truthcode);
		m_bj_pt->push_back(jet->pt());
		m_bj_eta->push_back(jet->eta());
		m_bj_phi->push_back(jet->phi());
		m_bj_m->push_back(jet->m());
        	//m_bj_status->push_back(-1);
	}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//HIGGS
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////Method 1. on shell Higgs

	// Clear Higgs particle vectors
	m_higgs_pt->clear();
	m_higgs_eta->clear();
	m_higgs_phi->clear();
	m_higgs_m->clear();
	m_higgs_pdgid->clear();
	m_higgs_status->clear();
	m_higgs_N = vec_higgs.size();

	// Loop over Higgs particles
	for(auto higgs : vec_higgs) {
		m_higgs_pt->push_back(higgs->pt());
		m_higgs_eta->push_back(higgs->eta());
		m_higgs_phi->push_back(higgs->phi());
		m_higgs_m->push_back(higgs->m());
		m_higgs_pdgid->push_back(higgs->pdgId());		
		m_higgs_status->push_back(higgs->status());
	}


//Method 2. Reconstruct off-shell higgs based on their decay products (i.e. 2.1 b QUARKS and 2.2 b JETS) 

	//2.0 initialise some variables
	m_onTimeHiggsReco_et = -999999;
        m_onTimeHiggsReco_pt = -999999;
        m_onTimeHiggsReco_pz = -999999;
        m_onTimeHiggsReco_eta = -999999;
        m_onTimeHiggsReco_phi = -999999;
        m_onTimeHiggsReco_m = -999999;

        m_outOfTimeHiggsReco_et = -999999;
        m_outOfTimeHiggsReco_pt = -999999;
        m_outOfTimeHiggsReco_pz = -999999;
        m_outOfTimeHiggsReco_eta = -999999;
        m_outOfTimeHiggsReco_phi = -999999;
        m_outOfTimeHiggsReco_m = -999999;

	m_onTimeHiggsPseudoApproach2_et = -999999;
        m_onTimeHiggsPseudoApproach2_pt = -999999;
        m_onTimeHiggsPseudoApproach2_pz = -999999;
        m_onTimeHiggsPseudoApproach2_eta = -999999;
        m_onTimeHiggsPseudoApproach2_phi = -999999;
        m_onTimeHiggsPseudoApproach2_m = -999999;
        m_onTimeHiggsPseudoApproach2_ProdVtx_t = -999999;
        m_onTimeHiggsPseudoApproach2_ProdVtx_Lxy = -999999;
        m_onTimeHiggsPseudoApproach2_ProdVtx_Lz = -999999;

        m_outOfTimeHiggsPseudoApproach2_et = -999999;
        m_outOfTimeHiggsPseudoApproach2_pt = -999999;
        m_outOfTimeHiggsPseudoApproach2_pz = -999999;
        m_outOfTimeHiggsPseudoApproach2_eta = -999999;
        m_outOfTimeHiggsPseudoApproach2_phi = -999999;
        m_outOfTimeHiggsPseudoApproach2_m = -999999;
        m_outOfTimeHiggsPseudoApproach2_ProdVtx_t = -999999;
        m_outOfTimeHiggsPseudoApproach2_ProdVtx_Lxy = -999999;
        m_outOfTimeHiggsPseudoApproach2_ProdVtx_Lz = -999999;

	//2.1 loop over all b quarks, form pairs	
	std::vector<std::vector<const xAOD::TruthParticle*>> vec_bbbarPairsCLUSTER;
        std::vector<float> vec_timingCLUSTER;

        int b1_idxCLUSTER = 0;
        for(auto b1 : vec_bquark){
                int b2_idxCLUSTER = 0;
                for(auto b2 : vec_bquark){
                        if(     b1->charge() * b2->charge() < 0                         &&
                                b1->prodVtx()->t() == b2->prodVtx()->t()        &&
                                b1->prodVtx()->x() == b2->prodVtx()->x()        &&
                                b1->prodVtx()->y() == b2->prodVtx()->y()        &&
                                b1->prodVtx()->z() == b2->prodVtx()->z()
                        ){
                                std::vector<const xAOD::TruthParticle*> pair;
                                pair.push_back(b1);
                                pair.push_back(b2);
                                if(b1_idxCLUSTER < b2_idxCLUSTER){
                                        vec_bbbarPairsCLUSTER.push_back(pair);
                                        vec_timingCLUSTER.push_back(b1->prodVtx()->t());
                                }
                        }
                b2_idxCLUSTER += 1;
                }
        b1_idxCLUSTER += 1;
        }

	int idxOnTime = -1;
        int idxOutOfTime = -1;
	std::vector<TLorentzVector*> vec_HiggsRECO(2);
	std::vector<TLorentzVector*> vec_HiggsCLUSTER(2);

	if (vec_bbbarPairsCLUSTER.size() != 2){
		ATH_MSG_WARNING("Cluster Algorithm found " << vec_bbbarPairsCLUSTER.size() << " pairs of opposite sign, same production vertex b quarks. We expect 2 pairs.");
	}
	else{
		//reconstruct Higgs candidates
		vec_HiggsRECO.at(0) = new TLorentzVector(vec_bbbarPairsCLUSTER.at(0).at(0)->p4() + vec_bbbarPairsCLUSTER.at(0).at(1)->p4());	
		vec_HiggsRECO.at(1) = new TLorentzVector(vec_bbbarPairsCLUSTER.at(1).at(0)->p4() + vec_bbbarPairsCLUSTER.at(1).at(1)->p4());
		
		//2.2 loop over jets, assign them to their closest b quark, and the corresponding bbbar pair	
		std::vector<std::vector<const xAOD::Jet_v1*>> vec_jetPairsCLUSTER(2);
		for(auto bjet : vec_bjet){
                	//search for closest b quarks in dR < 0.4:
                	float dR_Bjet_closestBquark = 999; //huge default value for delta R between this jet and its closest b quark
                	const xAOD::TruthParticle* closestBquark = nullptr; //set nullptr as initial value for closest b quark
                	float bj_eta = bjet->eta(); //eta of this jet
                	float bj_phi = bjet->phi(); //phi of this jet
                	for(auto bquark : vec_bquark){
                        	float bq_eta = bquark->eta(); //eta of quark
                        	float bq_phi = bquark->phi(); //phi of quark
                        	float dEta = bq_eta - bj_eta; //delta eta between jet and quark
                        	float dPhi = bq_phi - bj_phi; //delta phi between jet and quark
                        	float dR = sqrt(dEta * dEta + dPhi * dPhi); //delta R between jet and quark
                        	//"update" delta R between jet and closest quark, if the delta R between current jet and quark is smaller than the previous value
                        	if (dR < dR_Bjet_closestBquark){
                                	dR_Bjet_closestBquark = dR;
                                	closestBquark = bquark;
                        	}
                	}
                	//if there is a quark within dR < 0.4, we assign the jet to the corresponding bbbar pair
                	if (dR_Bjet_closestBquark < .4){
                		int idx = 0;
				for(auto bbbarpair : vec_bbbarPairsCLUSTER){
					for(auto bquark : bbbarpair){
						//check of bquark is the same as closestBquark		
						if(
							bquark->pt()  == closestBquark->pt()            &&
                                        		bquark->eta() == closestBquark->eta()           &&
                                        		bquark->phi() == closestBquark->phi()           
														)
						{
							vec_jetPairsCLUSTER.at(idx).push_back(bjet);

						}
					}
					idx += 1;
				}
			}
        	}
		for (int i = 0; i < 2; ++i) {
                        if (vec_jetPairsCLUSTER.at(i).size() > 0){
                                auto higgs = vec_jetPairsCLUSTER.at(i).at(0)->p4();
                                for (size_t j = 1; j < vec_jetPairsCLUSTER.at(i).size(); ++j) {
                                        higgs += vec_jetPairsCLUSTER.at(i).at(j)->p4();
                                }
                                vec_HiggsCLUSTER.at(i) = new TLorentzVector(higgs);
                        }
			else	ATH_MSG_WARNING("Cluster Algorithm has found a b quark pair without assigned jets.");
		}

		//Fill some variables
		
		if (vec_timingCLUSTER.at(0) < vec_timingCLUSTER.at(1)){
			idxOnTime = 1;
			idxOutOfTime = 0;
		}
		else if (vec_timingCLUSTER.at(0) > vec_timingCLUSTER.at(1)){
			idxOnTime = 0;
			idxOutOfTime = 1;
		}

		if (idxOnTime >= 0 && vec_HiggsRECO.at(idxOnTime) != nullptr){
                        m_onTimeHiggsReco_et = vec_HiggsRECO.at(idxOnTime)->Et();
                        m_onTimeHiggsReco_pt = vec_HiggsRECO.at(idxOnTime)->Pt();
                        m_onTimeHiggsReco_pz = vec_HiggsRECO.at(idxOnTime)->Pz();
                        m_onTimeHiggsReco_eta = vec_HiggsRECO.at(idxOnTime)->Eta();
                        m_onTimeHiggsReco_phi = vec_HiggsRECO.at(idxOnTime)->Phi();
                        m_onTimeHiggsReco_m = vec_HiggsRECO.at(idxOnTime)->M();
                        //m_onTimeHiggsReco_ProdVtx_t = vec_timingCLUSTER.at(idxOnTime); //vec_RecoHiggs.at(onTimeRECO)->prodVtx()->t(); //FIXME: Check if this works
                        //m_onTimeHiggsReco_ProdVtx_Lxy = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->perp(); //FIXME: Check if this works
                        //m_onTimeHiggsReco_ProdVtx_Lz = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->z(); //FIXME: Check if this works
                }
		if (idxOutOfTime >= 0 && vec_HiggsRECO.at(idxOutOfTime) != nullptr){
                        m_outOfTimeHiggsReco_et = vec_HiggsRECO.at(idxOutOfTime)->Et();
                        m_outOfTimeHiggsReco_pt = vec_HiggsRECO.at(idxOutOfTime)->Pt();
                        m_outOfTimeHiggsReco_pz = vec_HiggsRECO.at(idxOutOfTime)->Pz();
                        m_outOfTimeHiggsReco_eta = vec_HiggsRECO.at(idxOutOfTime)->Eta();
                        m_outOfTimeHiggsReco_phi = vec_HiggsRECO.at(idxOutOfTime)->Phi();
                        m_outOfTimeHiggsReco_m = vec_HiggsRECO.at(idxOutOfTime)->M();
                        //m_outOfTimeHiggsReco_ProdVtx_t = vec_timingCLUSTER.at(idxOutOfTime); //vec_RecoHiggs.at(onTimeRECO)->prodVtx()->t(); //FIXME: Check if this works
                        //m_outOfTimeHiggsReco_ProdVtx_Lxy = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->perp(); //FIXME: Check if this works
                        //m_outOfTimeHiggsReco_ProdVtx_Lz = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->z(); //FIXME: Check if this works
                }
        	if (	idxOnTime >= 0 && vec_HiggsRECO.at(idxOnTime) != nullptr &&
			idxOutOfTime >= 0 && vec_HiggsRECO.at(idxOutOfTime) != nullptr
		){
                	m_HiggsReco_del_phi = std::abs(vec_HiggsRECO.at(idxOnTime)->DeltaPhi(*vec_HiggsRECO.at(idxOutOfTime)));
        	}
		if (idxOnTime >= 0 && vec_HiggsCLUSTER.at(idxOnTime) != nullptr){
			m_onTimeHiggsPseudoApproach2_et = vec_HiggsCLUSTER.at(idxOnTime)->Et();
                        m_onTimeHiggsPseudoApproach2_pt = vec_HiggsCLUSTER.at(idxOnTime)->Pt();
                        m_onTimeHiggsPseudoApproach2_pz = vec_HiggsCLUSTER.at(idxOnTime)->Pz();
                        m_onTimeHiggsPseudoApproach2_eta = vec_HiggsCLUSTER.at(idxOnTime)->Eta();
                        m_onTimeHiggsPseudoApproach2_phi = vec_HiggsCLUSTER.at(idxOnTime)->Phi();
                        m_onTimeHiggsPseudoApproach2_m = vec_HiggsCLUSTER.at(idxOnTime)->M();
                        //m_onTimeHiggsPseudoApproach2_ProdVtx_t = vec_bbbarPairs_Approach2.at(onTime_Approach2).at(0)->prodVtx()->t();
                        //m_onTimeHiggsPseudoApproach2_ProdVtx_Lxy = vec_bbbarPairs_Approach2.at(onTime_Approach2).at(0)->prodVtx()->perp();
                        //m_onTimeHiggsPseudoApproach2_ProdVtx_Lz = vec_bbbarPairs_Approach2.at(onTime_Approach2).at(0)->prodVtx()->z();
		}
		if (idxOutOfTime >= 0 && vec_HiggsCLUSTER.at(idxOutOfTime) != nullptr){
			m_outOfTimeHiggsPseudoApproach2_et = vec_HiggsCLUSTER.at(idxOutOfTime)->Et();
                        m_outOfTimeHiggsPseudoApproach2_pt = vec_HiggsCLUSTER.at(idxOutOfTime)->Pt();
                        m_outOfTimeHiggsPseudoApproach2_pz = vec_HiggsCLUSTER.at(idxOutOfTime)->Pz();
                        m_outOfTimeHiggsPseudoApproach2_eta = vec_HiggsCLUSTER.at(idxOutOfTime)->Eta();
                        m_outOfTimeHiggsPseudoApproach2_phi = vec_HiggsCLUSTER.at(idxOutOfTime)->Phi();
                        m_outOfTimeHiggsPseudoApproach2_m = vec_HiggsCLUSTER.at(idxOutOfTime)->M();
                        //m_onTimeHiggsPseudoApproach2_ProdVtx_t = vec_bbbarPairs_Approach2.at(onTime_Approach2).at(0)->prodVtx()->t();
                        //m_onTimeHiggsPseudoApproach2_ProdVtx_Lxy = vec_bbbarPairs_Approach2.at(onTime_Approach2).at(0)->prodVtx()->perp();
                        //m_onTimeHiggsPseudoApproach2_ProdVtx_Lz = vec_bbbarPairs_Approach2.at(onTime_Approach2).at(0)->prodVtx()->z();
                }
		if (    idxOnTime >= 0 && vec_HiggsCLUSTER.at(idxOnTime) != nullptr &&
			idxOutOfTime >= 0 && vec_HiggsCLUSTER.at(idxOutOfTime) != nullptr
		){
                	m_HiggsPseudoApproach2_del_phi = std::abs(vec_HiggsCLUSTER.at(idxOnTime)->DeltaPhi(*vec_HiggsCLUSTER.at(idxOutOfTime)));;
        	}
	}//end else

	
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Hardware Trigger Decision (based on MET in N-1 and jet in N)
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

	//We mimic the TT strucutre of the hardware trigger
	//details on noise cuts, and TT strucutre https://gitlab.cern.ch/atlas-l1calo/firfilteranalysisrun3/-/blob/tobias-working-branch/source/AutoCorAnalysis/AutoCorAnalysis/mapping.h?ref_type=heads
	
	//TTs are typ. 0.1x0.1 in eta x phi. In the forward region (eta > 2.5) the Trigger Tower size in eta becomes a bit larger: 
	const unsigned short ETA_BINS = 66;
	const float etaBinCuts[ETA_BINS+1] = {
		-4.9,-4.475,-4.050,-3.625,-3.2,
		-3.1,-2.9,-2.7,-2.5,
		-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,
		-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
		0.0,
		0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,
		1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,
		2.7,2.9,3.1,
  		3.2,3.625,4.050,4.475,4.9
	};
	//phi is uniformely binned /from -3.2 to 3.2 in 0.1 steps
	const unsigned short PHI_BINS = 64;
	float phiBinCuts[PHI_BINS+1] ;
	for (int i = 0 ; i < PHI_BINS+1 ; ++i)		phiBinCuts[i] = -3.2 + i * 0.1;
	
	//Define vector to store all truth particles within the different TT
	std::vector<std::vector<std::vector<std::vector<const xAOD::TruthParticle*>>>> vec_TriggerTowerParticles(
		2, std::vector<std::vector<std::vector<const xAOD::TruthParticle*>>>( 					//ECAL and HCAL
    			ETA_BINS, std::vector<std::vector<const xAOD::TruthParticle*>>(					//Eta bins
        			PHI_BINS, std::vector<const xAOD::TruthParticle*>(					//Phi bins
		))));

	//The total energy per Trigger Tower must be higher than a threshold (eta dependent as defined below) in order to not account for noise
	double Etreshold[ETA_BINS] = {
	    11500, 10000, 7500, 4000, 2500, 3000, 3200, 3500,
	    1850, 1850, 1850, 1850, 1800, 1750, 1750, 1750,
	    1800, 1800, 1900, 1900, 2000, 2000, 2000, 2050,
	    2100, 2150, 2200, 2200, 2250, 2300, 2300, 2300,
	    2300, 2300,
	    2300, 2300, 2300, 2250, 2200, 2200, 2150, 2100,
	    2050, 2000, 2000, 2000, 1900, 1900, 1800, 1800,
	    1750, 1750, 1750, 1800, 1850, 1850, 1850, 1850,
	    3500, 3200, 3000, 2500, 4000, 7500, 10000, 11500
	}; //MeV
	//FIXME: Currently the PPr noisecuts of 2024 are used. Also, the same noise cuts are used for ECAL and HCAL layer.
	
	//List of pdgids that are invisible in the Calorimeters:
	std::vector<int> detectorEscape =      {12, //electron neutrino
						14, //muon neutrino
						16, //tau neutrino
						13, //muon
						m_pdgIdBSM, //BSM Higgs/Scalar
						m_pdgIdHS, //Chi1 
						m_pdgIdchi0 //Chi0
	};

	bool digiTime = true;  //if set to true, times (to differentiate between N-1 and N) are corrected (synchronised) for readout: t -> t - r/c with r being the travel distance

	//1. We assign all detector stable truth particles that are visible in the Calorimeter, and that have a production vertex in the tracker, ECAL, or HCAL, to the corresponding TTs 
	for (auto part : vec_finalStateParticles){ 												//vec_finalStateParticles includes all particles with status = 1
    		if (std::find(detectorEscape.begin(), detectorEscape.end(), std::abs(part->pdgId())) != detectorEscape.end()) continue;		//skip particles that are invisible to the CALO, e.g. Chi0, muons, neutrinos
		if (!part->hasProdVtx())									continue; 			//this should be unnecessary, as final state particles always should have a production vertex. 
		if (std::abs(part->prodVtx()->perp()) > 3500 or std::abs(part->prodVtx()->z()) > 5500)      	continue;			//skip particles that are produced outside the CALOs

		//In order to assign these truth particles to the corresponding we need eta and phi in the CALO
		//For particles that are produced in the CALO (not the tracker!) we could simply take the eta and phi of the production vertex
		//For particles produced in the tracker, we have to extrapolate the eta, phi to their entry into the CALO. Ehis is implemented by (simply geometry) a linear extrapolation [assuming that the tracks have no curvature, which should be fine at the GeV level] 	
		auto CaloEntry = extrapolate(part); //this returns the (x,y,z,t,t_readout) of the particle either at the production vertex (when inside the CALO) or at the extrapolated CALO entry
		if (CaloEntry.at(3) == 1e10){
			ATH_MSG_WARNING("Extrapolation failed, e.g. due to vx = vy = vz = 0");
			continue;
		}
		float phi = std::atan2(CaloEntry.at(1) , CaloEntry.at(0)); //based on x,y in the CALO
		float r = std::sqrt(CaloEntry.at(0) * CaloEntry.at(0) + CaloEntry.at(1) * CaloEntry.at(1) + CaloEntry.at(2) * CaloEntry.at(2));
		float eta = 0.5 * std::log((r + CaloEntry.at(2)) / (r - CaloEntry.at(2))); //based on x,y,z in the CALO

		//now we assign it to the corresponding Trigger Tower in eta/phi/ ECAL and HCAL respectively
		//1.i) assign the corresponding eta bin
		int etaBin = -1;
		for (int i = 0; i < ETA_BINS; ++i) {
    			if (eta >= etaBinCuts[i] && eta < etaBinCuts[i+1]) {
			        etaBin = i;
        			break;
    			}
		}
		if (etaBin == -1) 		continue; //the particle is outside the L1CALO coverage. Skip the particle

		//1.ii) wrapp phi into -pi,pi interval. Probably not needed, as per definition of phi=std::atan2(...) phi should already be in the "correct" interval -pi,pi 
		int phiBin = -1;
		if (phi > M_PI) phi -= 2 * M_PI; 
		if (phi <= -M_PI) phi += 2 * M_PI; 
		//1.iii) assign the correspinding phi bin
		for (int i = 0; i < PHI_BINS; ++i) {
                	if (phi >= phiBinCuts[i] && phi < phiBinCuts[i+1]) {       
	       			phiBin = i;
                                break;
                        }
		}
	
		//1.iv) assing particles with production vertex in tracker or the ECAL to the corresponding ECAL TT, and particles produced in the HCAL to the HCAL TT
		float rho = part->prodVtx()->perp();
		float z = part->prodVtx()->z();
		bool isEM = rho < 2000 and z < 4300; //rough values of the outer boundary of the ECAL. I.e. all particles with prod Vtx in tracker or ECAL are added to the TT of the ECAL layer
		vec_TriggerTowerParticles.at(isEM).at(etaBin).at(phiBin).push_back(part);	
	}


	///////////////////////////
	///////////////////////////
	// MET in N-1
	///////////////////////////
	///////////////////////////
	
	//2. We mimic the MET algorithm of the Hardware Trigger
	//Therefore approximate the MET in N-1 as the negative ET in N-1. 
	//Hence, we calculate the momentum balance of all visible particles within the time window N-1	
	float momentumBalance_x = 0 , momentumBalance_y = 0 , momentumBalance_z = 0;
	//for comparison, we also calculate the energy balance by weighting the total energy with sin / cos phi
        float energyBalance_x = 0 , energyBalance_y = 0;
	
	//2.i) loop over all TTs, and all final state particles in this TT, add their energies, and store the energy/ momenta sums if the energy exceeds the noise cut
	for (int calo = 0 ; calo <= 1 ; ++calo){ //ECAL and HCAL layer
		for (int etaBin = 0 ; etaBin < ETA_BINS ; ++etaBin){
			float threshold = float(Etreshold[etaBin]); //FIXME: currently same noise cuts for HCAL and ECAL
			for (int phiBin = 0 ; phiBin < PHI_BINS ; ++phiBin){
				//We estimate phi by the center of the TT
				float phiLow = phiBinCuts[phiBin];	
				float phiHigh = phiBinCuts[phiBin+1];
				float phi = 0.5 * (phiLow + phiHigh);
				
				//we sum momenta and energies off all particles of the "current" TT:
				float px = 0;
				float py = 0;
				float pz = 0;
				float e = 0;
				for (auto part : vec_TriggerTowerParticles.at(calo).at(etaBin).at(phiBin)){
					//In order to assign the particles to N-1, we use the readout time of the (extrapolatet) CaloEntry
					auto CaloEntry = extrapolate(part);
					if (CaloEntry.at(3) == 1e10)			continue; //the particle has a production vertex outside the CALO
                                        float time = CaloEntry.at(3);
                                        if (digiTime){					//digiTime is just a bool to switch off the flight-time correction (digiTime=False only recommended for comparison/debugging)
                                                time = CaloEntry.at(4);
                                        }
					time = time / 3e2; //conversion from mm to ns
                                        if (time < 0 or time > 10)              continue; //CALO is sensitive to 0 < t < 10 in BC N-1
					//sum the momenta of all particles; eventually used to estimate momentum balance
					px += part->px();
					py += part->py();
					pz += part->pz();
					e += part->e(); //sum energies of all particles; eventually conmpared to noise cut
				}
				if (e >= threshold){ 
					//if TT exceeds the threshold, we store the summed momenta
					momentumBalance_x += px;
					momentumBalance_y += py;
					momentumBalance_z += pz;	
					//for comparison, we store also the summed energy weighted by sin phi and cos phi
					energyBalance_x += (e * std::sin(phi));
					energyBalance_y += (e * std::cos(phi));
				}
			}
		}
	}

	//2.ii) Calcualte MET 
	//momentum balance
	m_L1_MET_Nminus1_px = momentumBalance_x;
        m_L1_MET_Nminus1_py = momentumBalance_y;
        m_L1_MET_Nminus1_pz = momentumBalance_z;
	//MET
	m_L1_MET_Nminus1_mpt = std::sqrt(momentumBalance_x * momentumBalance_x + momentumBalance_y * momentumBalance_y); //MET calculatet on momentum balance
	m_L1_MET_Nminus1_met = std::sqrt(energyBalance_x * energyBalance_x + energyBalance_y * energyBalance_y); //MET based on energy balance (for comparison only)
        //Eta and Phi
	m_L1_MET_Nminus1_phi = std::atan2(-momentumBalance_y , -momentumBalance_x); //phi based on momentum balance
	float momentum = std::sqrt(momentumBalance_x * momentumBalance_x + momentumBalance_y * momentumBalance_y + momentumBalance_z * momentumBalance_z);
        m_L1_MET_Nminus1_eta = 0.5 * std::log((momentum + momentumBalance_z) / (momentum - momentumBalance_z));; //eta based on momentum balance

	//2.iii) For comparison, we calculate the momentum/energy balance of all invisible particles ("Truth MET") in N-1
	float truthMET_Nminus1_PX = 0;
	float truthMET_Nminus1_PY = 0;
        float truthMET_Nminus1_EX = 0; 
        float truthMET_Nminus1_EY = 0; 

	//a. neutrinos 
	for(auto nu : vec_nu) {
		if (! nu->hasProdVtx())							continue;
        	//again, we use the readout time of the CALO entry to divide the event into N-1 and N       	
		auto CaloEntry = extrapolate(nu);
		if (CaloEntry.at(3) == 1e10)			continue; //skip the particle if it is produced outside the CALO
                float time = CaloEntry.at(3);
                if (digiTime){ 
                	time = CaloEntry.at(4);
                }
                time = time / 3e2; //conversion mm to ns
                if (time < 0 or time > 10)              continue; //only take neutrinos in N-1, i.e. in 0-10 ns
		
		//add momenta to the overall momentum balance
		truthMET_Nminus1_PX += nu->px();
		truthMET_Nminus1_PY += nu->py();
		//also add energy balance for comparison/debugging reasons
		float phi = nu->phi(); //FIXME we should use the CaloEntry (x,y) coordinates to calculate phi 
		float e = nu->e();
		truthMET_Nminus1_EX += e * std::sin(phi);
		truthMET_Nminus1_EY += e * std::cos(phi);
	}

	//b. muons
        for(auto mu : vec_mu) {
                if (! mu->hasProdVtx())                                                 continue;
                //again, we use the readout time of the CALO entry to divide the event into N-1 and N
                auto CaloEntry = extrapolate(mu);
                if (CaloEntry.at(3) == 1e10)                    continue; //skip the particle if it is produced outside the CALO
                float time = CaloEntry.at(3);
                if (digiTime){
                        time = CaloEntry.at(4);
                }
                time = time / 3e2; //conversion mm to ns
                if (time < 0 or time > 10)              continue; //only take muons in N-1, i.e. in 0-10 ns

                //add momenta to the overall momentum balance
                truthMET_Nminus1_PX += mu->px();
                truthMET_Nminus1_PY += mu->py();
                //also add energy balance for comparison/debugging reasons
                float phi = mu->phi(); //FIXME  we should use the CaloEntry (x,y) coordinates to calculate phi
                float e = mu->e();
                truthMET_Nminus1_EX += e * std::sin(phi);
                truthMET_Nminus1_EY += e * std::cos(phi);
        }

	//c. BSM particles
	//we loop over the leading (i=0) and subleading Chi1s and Chi0s
	for (int i = 0 ; i < 2 ; ++i) {
		
		//leading/subleading chi1 and chi0
		auto chi1 = vec_llp.at(i);
		auto chi0 = vec_chi0_llp.at(i);

		//when the chi0 is produced within the tracker/CALO, and the readout time is within N-1, we consider the Chi0 to contribute to the MET. Otherwise the Chi1 is contributing instead (because the chi1->chi0+h decay has not happend in N-1)
		auto CaloEntry = extrapolate(chi0);
                float time = CaloEntry.at(3); 
                if (digiTime){
                        time = CaloEntry.at(4);
                }
                time = time / 3e2; //conversion from mm to ns
                //if (time < 0 or time > 10){
                if (time > 0 and time < 10){ 
			truthMET_Nminus1_PX += chi0->px();
                	truthMET_Nminus1_PY += chi0->py();
			//for comparison / debugging reasons we also compute the energy balance
                	float phi = chi0->phi(); //FIXME  we should use the CaloEntry (x,y) coordinates to calculate phi
                	float e = chi0->e();
                	truthMET_Nminus1_EX += e * std::sin(phi);
                	truthMET_Nminus1_EY += e * std::cos(phi);
		}
		else { //use chi1 instead of chi0
		        truthMET_Nminus1_PX += chi1->px();
                        truthMET_Nminus1_PY += chi1->py();
                        float phi = chi1->phi(); //FIXME  we should use the CaloEntry (x,y) coordinates to calculate phi
                        float e = chi1->e();
                        truthMET_Nminus1_EX += e * std::sin(phi);
                        truthMET_Nminus1_EY += e * std::cos(phi);
		}
	}

	//calculate MET
	m_truth_MET_Nminus1_mpt = std::sqrt(truthMET_Nminus1_PX * truthMET_Nminus1_PX + truthMET_Nminus1_PY * truthMET_Nminus1_PY ); //MET calculated as negative momentum balance (which is typically called MET)
	m_truth_MET_Nminus1_met = std::sqrt(truthMET_Nminus1_EX * truthMET_Nminus1_EX + truthMET_Nminus1_EY * truthMET_Nminus1_EY ); //MET calculated bases on the energy balance (untypical, will not be used in the further analysis, was just calculated for comparison / debugging)

	
	///////////////////////////
	///////////////////////////
	// jet in N
	///////////////////////////
	///////////////////////////
	
	//3. We mimic the small R jet algorithm of the Hardware Trigger
	//We search for local energy maxima, and add the TTs in an R=0.4 "cone"	around the maximum
	
	//Define vector to store the energy of each TT in the time window N
        std::vector<std::vector<std::vector<float>>> vec_TriggerTowerEnergies_N(
                2, std::vector<std::vector<float>>(                                     //ECAL and HCAL
                        ETA_BINS, std::vector<float>(                                   //Eta bins
                                PHI_BINS                                      		//Phi bins
                )));

	//3.i) We only consider TTs above the noise cut
	//Loop over all particles in the individual TTs in event N
	for (int calo = 0 ; calo <= 1 ; ++calo){ //ECAL and HCAL layer
                for (int etaBin = 0 ; etaBin < ETA_BINS ; ++etaBin){
                        float threshold = float(Etreshold[etaBin]); //FIXME: currently the same noisecuts are used for the ECAL and HCAL
                        for (int phiBin = 0 ; phiBin < PHI_BINS ; ++phiBin){
                                //We sum the energies/momenta of all particles in this TT
				float px = 0;
                                float py = 0;
                                float pz = 0;
                                float e = 0; 
                                for (auto part : vec_TriggerTowerParticles.at(calo).at(etaBin).at(phiBin)){
					//similar to the MET algorithm, we extrapolate the calo entry time, in order to only consider particles in the time window N
					auto CaloEntry = extrapolate(part);
					if (CaloEntry.at(3) == 1e10)			continue; //skip particles outside the tracker,ECAL,HCAL
                                        float time = CaloEntry.at(3);
                                        if (digiTime){
                                                time = CaloEntry.at(4); //correct entry time by the time of flight
                                        }
                                        time = time / 3e2; //conversion from mm to ns
                                        if (time < 25 or time > 35)              continue; //Sensitive time window of BC N is between 25 and 35 ns
					//consider the energy/momentum of this particle for the TT energy/momentum
					px += part->px();
                                        py += part->py();
                                        pz += part->pz();
                                        e += part->e();
                                }
				//compare TT energy to noise cut
				if (e >= threshold){ 
					//trigger tower exceeds threshold, we store the summed energy of all particles
                                        vec_TriggerTowerEnergies_N.at(calo).at(etaBin).at(phiBin) = e;
				}
                                else { 
					//trigger tower below threshold. Energy is set to 0
                                        vec_TriggerTowerEnergies_N.at(calo).at(etaBin).at(phiBin) = 0;
                                }
                        }
                }
        }

	//3.ii) loop over all TT, and for each eta x phi TT sum the ECAL and HCAL layer, and search for local maxima
	//local maxima are searched in a window of 5x5 Trigger Tower window (~0.5x0.5 in eta x phi) around the central TT
	//the central TT is classified as maximum, if its energie is larger than the 24 neighbouring TTs in the 5x5 search window. 
	//In this case the eta and phi bins are stored in the following list:
	std::vector<std::vector<int>> localMaxima;

	for (int etaBin = 0 ; etaBin < ETA_BINS ; ++etaBin){
		for (int phiBin = 0 ; phiBin < PHI_BINS ; ++phiBin){
			// the TTs with etaBin and phiBin are defined as central TT
			// sum energies of ECAL and HCAL of the central TT
			float e_centralTT = vec_TriggerTowerEnergies_N.at(0).at(etaBin).at(phiBin) + vec_TriggerTowerEnergies_N.at(1).at(etaBin).at(phiBin);
			//loop over all neighbours in 5x5 window around central TT
			//check whether the central TT has the highest energy in this window
			bool max = true; //set `max` initially to true 
			if (e_centralTT == 0)		max = false; //skip this central TT since it has not exceeded the noise cut
			//compare energy of central TT to all neigbouring TTs
			for (int slidingIndexETA = -2; slidingIndexETA <= 2 ; ++slidingIndexETA){ //5 eta bins from -2,...,0,...,+2
				int etaBinNEIGHBOUR = etaBin + slidingIndexETA;
				//at the edge of the calo there are less neighbours we can compare to:
				if (etaBinNEIGHBOUR < 0 or etaBinNEIGHBOUR >= ETA_BINS)			continue;	
				for (int slidingIndexPHI = -2; slidingIndexPHI <= 2 ; ++slidingIndexPHI){	 //5 phi bins from -2,...,0,...,+2
					int phiBinNEIGHBOUR = phiBin + slidingIndexPHI;
	                                //use the Cyclic of phi at the edge of the CALO
					if (phiBinNEIGHBOUR >= PHI_BINS) phiBinNEIGHBOUR -= PHI_BINS;                   
	                                if (phiBinNEIGHBOUR < 0) phiBinNEIGHBOUR += PHI_BINS;                           
					//we don't want to compare the central TT with itself:
					if (etaBinNEIGHBOUR == etaBin and phiBinNEIGHBOUR == phiBin)	continue;
					//sum the ECAL and HCAL layers of the "current" neighbouring TT					
					float e_neighbourTT = vec_TriggerTowerEnergies_N.at(0).at(etaBinNEIGHBOUR).at(phiBinNEIGHBOUR) + vec_TriggerTowerEnergies_N.at(1).at(etaBinNEIGHBOUR).at(phiBinNEIGHBOUR);
					//compare the central and the neighbouring TTs
					if (e_neighbourTT > e_centralTT)		max = false;	
				} 
			}
			//in case max = true the central TT has the highest energy in this window, and we store the central TT as a local maximum
			if (max){
				localMaxima.push_back({etaBin,phiBin});
				//ATH_MSG_ALWAYS("Found a local maximum (" << etaBin << "," << phiBin << ")!");
			}
		}
	}

	//3.iii) now we form smallR jets by summing all TT in a 9x9 window (actually not all, but those in a specific pattern to mimic a circular shape) around the local maxima
	//we store the (energy,eta,phi) values in the following list
	std::vector<std::tuple<float, float, float>> jets;
	for (auto max : localMaxima){
		int etaBin = max.at(0);
		int phiBin = max.at(1);
		//within a 9x9 window a small jet is defined by the following TT pattern
		// 0	0	0	0	0	0	0	0	0
		// 0	0	1	1	1	1	1	0	0
		// 0	1	1	1	1	1	1	1	0
		// 0	1	1	1	1	1	1	1	0
		// 0	1	1	1	1	1	1	1	0
		// 0	1	1	1	1	1	1	1	0
		// 0	1	1	1	1	1	1	1	0
		// 0	1	1	1	1	1	1	1	0
		// 0	0	1	1	1	1	1	0	0
		// 0	0	0	0	0	0	0	0	0
		//see e.g. fig. 3.3 in https://cds.cern.ch/record/2688511/files/ATL-COM-DAQ-2019-146.pdf
		//compute the energy of the jet by summing the corresponding TTs
		float e_smallRjet = 0;
		for (int slidingIndexETA = -4; slidingIndexETA <= 4 ; ++slidingIndexETA){ //9 bins in eta
                        int etaBinNEIGHBOUR = etaBin + slidingIndexETA;
			if (etaBinNEIGHBOUR < 0 or etaBinNEIGHBOUR >= ETA_BINS)         		continue; //on the edge of the calo there are less TTs to be considered
                        for (int slidingIndexPHI = -4; slidingIndexPHI <= 4 ; ++slidingIndexPHI){ //9 bins in phi
                                int phiBinNEIGHBOUR = phiBin + slidingIndexPHI;
				//Use cyclic of phi
				if (phiBinNEIGHBOUR >= PHI_BINS) phiBinNEIGHBOUR -= PHI_BINS; 			
                		if (phiBinNEIGHBOUR < 0) phiBinNEIGHBOUR += PHI_BINS;				
				//only use TT in the pattern as defined above
                                if (std::abs(slidingIndexETA) == 4 or  std::abs(slidingIndexPHI) == 4)		continue; //skip TT at the edge of the 9x9 window
                                if (std::abs(slidingIndexETA) == 3 and std::abs(slidingIndexPHI) == 3)		continue; //skip TT at the "corner" of the remaining 8x8 window
				e_smallRjet += vec_TriggerTowerEnergies_N.at(0).at(etaBinNEIGHBOUR).at(phiBinNEIGHBOUR) + vec_TriggerTowerEnergies_N.at(1).at(etaBinNEIGHBOUR).at(phiBinNEIGHBOUR);
                        }
                }

		//take eta and phi in the center of the central TT
		float etaLow = etaBinCuts[etaBin];
                float etaHigh = etaBinCuts[etaBin+1];
                float eta = 0.5*(etaLow+etaHigh);

		float phiLow = phiBinCuts[phiBin];
                float phiHigh = phiBinCuts[phiBin+1];
                float phi = 0.5*(phiLow+phiHigh);

		ATH_MSG_ALWAYS("L1 JET in N with energy " << e_smallRjet << " MeV, eta " << eta << ", and phi " << phi << " (N.B. phi of MET in N-1 was " << m_L1_MET_Nminus1_phi << ")");

		//jet acceptance (at the level of our hardware trigger) is |eta| < 3.2
		if (etaLow > -3.2 and etaHigh < 3.2) { 
                        jets.emplace_back(e_smallRjet, eta, phi);
                }
	}

	//3.iv) store some variables
	//sort the jets by descending energy and keep only the top 6 (our trigger was implemented in such a way to store the six most energetic jets)
	std::sort(jets.begin(), jets.end(),
 		[](const auto& a, const auto& b) {	
		return std::get<0>(a) > std::get<0>(b); 
    	});
        if (jets.size() > 6)            jets.resize(6);
	
	//We need at least on "leading" jet with ET > 40 GeV and delta Phi(jet,MET) < 1.0
        bool leading = false; //keep track if leading jet was found

        // Now fill the output vectors
	m_L1_jet_N_e->clear();
	m_L1_jet_N_et->clear();
        m_L1_jet_N_eta->clear();
        m_L1_jet_N_phi->clear();
        m_L1_delPhi->clear();

	m_L1_leadingJet_N_e = -1;
        m_L1_leadingJet_N_et = -1;
        m_L1_leadingJet_N_eta = -100;
        m_L1_leadingJet_N_phi = -10;
        m_L1_leading_delPhi = -1;

        for (const auto& jet : jets) {
                float e = std::get<0>(jet);
                float eta = std::get<1>(jet);
                float et = e / std::cosh(eta);
                float phi = std::get<2>(jet);

		//calculate Delta phi between the jet (in N) and the MET in N-1
                TLorentzVector L1_jet;
                L1_jet.SetPtEtaPhiM(-1., eta, phi, -1.);
                TLorentzVector L1_MET;
                L1_MET.SetPtEtaPhiM(-1., 0 , m_L1_MET_Nminus1_phi, -1.);
                float delPhi = std::abs(L1_MET.DeltaPhi(L1_jet));

		//store the most energetic jet with ET > 40 and delta Phi < 1.0
		if ((not leading) and (et/1000 > 40) and (delPhi < 1.0)){
                        m_L1_leadingJet_N_e = e;
                        m_L1_leadingJet_N_et = et;
                        m_L1_leadingJet_N_eta = eta;
                        m_L1_leadingJet_N_phi = phi;
                        m_L1_leading_delPhi = std::abs(delPhi);
                        leading = true;
                }

		//store all (up to six) jets
                m_L1_jet_N_e->push_back(e);
                m_L1_jet_N_et->push_back(e/std::cosh(eta));
                m_L1_jet_N_eta->push_back(eta);
                m_L1_jet_N_phi->push_back(phi);
                m_L1_delPhi->push_back(std::abs(delPhi));
        }

        m_L1_jet_N_N = m_L1_jet_N_e->size();

	///////////////////////////
	///////////////////////////
	// High Level Trigger in N
	///////////////////////////
	///////////////////////////

	//4. mimic the HLT decision and offline RECO of Event N
	//in N we need an 0.4 Jet within eta < 2.5 with pT > 20
	//furthermore, a significant fraction of the jet must be in the HCAL (E_EM < 6.4% E_HAD)
	//finally, there mustn't be tracks with pT > 2 GeV in a 0.2 cone around the jet axis

	//Tracks in N
	//first, we define tracks as charged TruthParticles, with production vertex in BC N (readout Time between 25 and 35 ns), and inside the inner detector
	std::vector<const xAOD::TruthParticle*> particlesWithTracks_N;
	for (auto part : vec_finalStateParticles){ // vec_finalStateParticles loops over all detector stable particles

		if (part->charge() == 0)		continue; //skip neutral particles
		if (! part->hasProdVtx())		continue; //this condition should be trivial for detector stable particles
		
		//we approximate that the particle is readout at its production vertex (no extrapolation as above)
		auto vtx = part->prodVtx();
		float t0 = vtx->t();
                float rho0 = vtx->perp();
                float z0 = vtx->z();
                float r0 = std::sqrt(rho0 * rho0 + z0 * z0);
		//we correct the time t0 by the distance to the interaction vertex r0:
                float readoutTime = (t0 - r0); 
                readoutTime = readoutTime / 3e2; //conversion from mm to ns
		
		float time = t0;
		if (digiTime)		time = readoutTime;
                if (time < 25 or time > 35)       continue; //only consider tracks in N //FIXME does tracker have the same time resolution as CALO, or is the window 25-x larger/smaller than in the CALO?	
		if (std::abs(z0) > 2500 or rho0 > 1000)		continue; //only consider particles in the inner detector
		
		//the surviving particles are defined as tracks
		particlesWithTracks_N.push_back(part);
	}
	if (particlesWithTracks_N.size() != 0){
		ATH_MSG_ALWAYS("There are " << particlesWithTracks_N.size() << " tracks in the considered time and space window");
	}

	//mimic the HLT decision and offline RECO of Event N
        m_offline_jet_N_e->clear();
        m_offline_jet_N_et->clear();
        m_offline_jet_N_eta->clear();
        m_offline_jet_N_phi->clear();

	//Now we use antiKt4 truth jets to check whether at least one jet fulfils all criteria, i.e.
	//i) pT > 20, |eta| < 2.5
	//ii) E_EM < 6.4% E_HAD
	//iii) no tracks with pT > 2 GeV in 0.2 cone around jet axis
	
	for(auto jet : vec_jet){ //vec_jet contains all AntiKt4TruthDressedWZJets, see e.g. https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TruthDAOD

		float jet_pt = jet->pt(); 
		if (jet_pt < 20000)				continue; //i) pT threshold of 20 GeV		
		float jet_e = jet->e();

		//We use all Jet Constitutents and their production vertices to define the jet origin
		//We access the JC as xAOD::IParticle and downcast it to access the xAOD::TruthParticle information
		auto constitutents = jet->getConstituents().asIParticleVector(); //this has the type std::vector<const xAOD::IParticle *>
		//initialise an empty list for the downcasted truth particles and their production vertices
		std::vector<const xAOD::TruthParticle*> JC_truthParticle_list; 
		std::vector<const xAOD::TruthVertex_v1*> prodVtx_list; 
		for (const xAOD::IParticle* constituent : constitutents){
			const xAOD::TruthParticle* part = dynamic_cast<const xAOD::TruthParticle*>(constituent);
			if (! part->hasProdVtx())		continue;	
			JC_truthParticle_list.push_back(part);
			//add the production vertex to the list if it is not on the list so far
			auto vtx0 =  part->prodVtx();	
			bool found = false; //keep track if this vtx is already in the prodVtx_list
			for (const auto* existing_vtx : prodVtx_list) {
			    if (existing_vtx == vtx0) {
			        found = true;
			        break;
			    }
			}
			if (!found) {
			    prodVtx_list.push_back(vtx0);
			}
		}

		if (prodVtx_list.empty())			continue;
		//We define the jet origin as the JC production vertex with the smallest production time and a momentum fraction (i.e. pT_vertex/pT_jet) of at least 30% //FIXME how sensible is this value, and this approach?
		//Therefore, we sort the production vertices based on their readout time 
		//FIXME should I sort them based on momentum fraction instead (i.e. p_vertex/p_jet)? I see arguments for both sortings.
		std::sort(prodVtx_list.begin() , prodVtx_list.end() ,
				[] (const xAOD::TruthVertex_v1* a , 
				    const xAOD::TruthVertex_v1* b) {
				return (a->t() < b->t()); //sort vertices by time
				//return (a->t() - std::sqrt(a->perp() * a->perp() + a->z() * a->z())) < (b->t() - std::sqrt(b->perp() * b->perp() + b->z() * b->z())); //sort vertices by readout time
				});
		//We search for the "first" prodVtx with momentum fraction > 0.3
		int idx = 0;
		for (auto vtx : prodVtx_list){
			idx += 1;
			//We caclulate the momentum fractions for i) pT sum of all outgoing particles of this vertex in comparison to pT of the jet and ii) pT sum of all outgoing particles of this vertex, which additionally are jet constitutents of the jet in comparison to the jet pT
                	float PT_all = 0 , PT = 0; //PT_all is used for i) and PT is used for ii)
			//loop over all outgoing particles of this vertex
                        for (int i = 0 ; i < (int) vtx->nOutgoingParticles() ; ++i){
                        	auto part = vtx->outgoingParticle(i);
                                PT_all += part->pt(); //add the pT regardless if this outgoing particle is part of the jet or not
                                if (std::find(JC_truthParticle_list.begin(), JC_truthParticle_list.end(), part) != JC_truthParticle_list.end()) {
                                        PT += part->pt(); //add the pT only if this outgoing particle is part of the jet
                                }
                        }
			//take the first vertex with pT(vertex)/pT(jet) > 0.3
			if (PT_all / jet_pt > .3){
				break; 
			}
		}	

		//grab this vertex. If none of the vertices has a momentum fraction > 30%, we currently just take the last vertex. //FIXME what to do instead? Take the first vertex or maybe the vertex with highest momentum fraction?
		auto vtx0 = prodVtx_list.at(idx - 1);	
		float t0 = vtx0->t();
		float rho0 = vtx0->perp();
		float z0 = vtx0->z();
		float r0 = std::sqrt(rho0 * rho0 + z0 * z0);
		float readoutTime = (t0 - r0);
		//ATH_MSG_ALWAYS("Originated at " << t0 << " " << rho0 << " " << z0 << " " << readoutTime);
		readoutTime = readoutTime / 3e2; //no extrapolation needed, as EMF cut requires prod Vtx within CALOs
			
		//calculate eta and phi of jet axis based on (x,y,z) of the jet origin
		float phi0 = std::atan2(vtx0->y() , vtx0->x());
		float eta0 = 0.5 * std::log((r0 + z0) / (r0 - z0));
		
		float jet_et = jet_e / std::cosh(eta0);
	
		//only consider jets in the time Window of BC N	
		float time = t0;
		if (digiTime)		time = readoutTime;	
		if (time < 25 or time > 35)	continue;
		ATH_MSG_ALWAYS("Survived timing cut");

		//i) only consider jets in eta < 2.5 (this is the acceptance of the high level trigger)
		if (std::abs(eta0) > 2.5)			continue;

                ATH_MSG_ALWAYS("HLT JET in N with energy " << jet_e << " MeV, eta " << eta0 << ", and phi " << phi0);

		//ii) as a guesstimate, we only consider jets that originate at least one hadronic interaction length away from the HCAL entry (which should approximate the E_EM < 6.4% E_HAD condition)
		if (rho0 > 3500 or std::abs(z0) > 5500)		continue; //skip jets originating outside the HCAL
		if (rho0 < 1400 and std::abs(z0) < 3700)	continue; //skip jets originating in the tracker
		//For jets originating in the ECAL, we calculate the remaining interaction length within the ECAL
		float lambda = 0;
		if ((std::abs(z0) > 3700 and std::abs(z0) < 4300 and rho0 < 2000) or (rho0 > 1400 and rho0 < 2000 and std::abs(z0) < 3700)){
			lambda = remainingInteractionLength(eta0 , rho0 , z0);
		}
		ATH_MSG_ALWAYS("  >> Remaining hadronic interaction length in ECAL: " << lambda << " (production vertex at rho0 " << rho0 << " mm, and z0: " << z0 << " mm)");
		if (std::abs(lambda) >= 1.)			continue; //skip jets with more than 1 interaction length in the ECAL
		//jets originating in the HCAL automatically pass the EMF cut
		ATH_MSG_ALWAYS("Survived EMF cut");

		//iii) lastly, we apply the isolation cut	
		bool foundTrack = false;
		for (auto part : particlesWithTracks_N){
			//loop over charged particles with production vertex in BC N, inside the tracker
			//we are only worried about tracks with pT > 2GeV		
			if (part->pt() < 2000) 			continue;
			auto trackVtx = part->prodVtx();
			float etaTrack = trackVtx->v4().Eta(); //We define eta based on the production vertex
			float phiTrack = trackVtx->v4().Phi(); //We define phi based on the production vertex

			//Define jet axis based on eta and phi of the jet origin (see above)
			TLorentzVector jetAxis;
               		jetAxis.SetPtEtaPhiM(-1., eta0, phi0 , -1.);
			//Define Track based on eta and phi of particle's production vertex
               		TLorentzVector track;
               		track.SetPtEtaPhiM(-1., etaTrack , phiTrack , -1.);
			//Calculate Delta R
			float delR = std::abs(jetAxis.DrEtaPhi(track));
			//A critical track is found if Delta R < 0.2
			if (delR > 0.2)				continue;
			//if the loop survives until here, we have found a track with pT > 2 GeV in a cone of dR = 0.2 around the jet axis
			foundTrack = true;
		}

		//if a critical track was found, we skip this jet
		if (foundTrack)					continue;
		ATH_MSG_ALWAYS("Survived isolation cut");

		//all cuts are survived
		//if the loop survives until here, we declare the jet as a "good" offline jet
	        m_offline_jet_N_e->push_back(jet_e);
	        m_offline_jet_N_et->push_back(jet_et); 
	        m_offline_jet_N_eta->push_back(eta0);
       		m_offline_jet_N_phi->push_back(phi0);
	} //end of jet-loop
        m_offline_jet_N_N = m_offline_jet_N_e->size();

	tree ("analysis")->Fill ();

	return StatusCode::SUCCESS;

}




StatusCode MCValAlg :: finalize ()
{
	// This method is the mirror image of initialize(), meaning it gets
	// called after the last event has been processed on the worker node
	// and allows you to finish up any objects you created in
	// initialize() before they are written to disk.  This is actually
	// fairly rare, since this happens separately for each worker node.
	// Most of the time you want to do your post-processing on the
	// submission node after all your histogram outputs have been
	// merged.
	return StatusCode::SUCCESS;
}

MCValAlg :: ~MCValAlg () {

	// Delete the allocated vectors to avoid memory leaks
	delete m_scalar_pt;
	delete m_scalar_eta;
	delete m_scalar_phi;
	delete m_scalar_m;
	delete m_scalar_pdgid;
	delete m_scalar_status;

	delete m_higgs_pt;
	delete m_higgs_eta;
	delete m_higgs_phi;
	delete m_higgs_m;
	delete m_higgs_pdgid;
	delete m_higgs_status;

	delete m_bsm_pt;
	delete m_bsm_eta;
	delete m_bsm_phi;
	delete m_bsm_m;
	delete m_bsm_pdgid;
	delete m_bsm_status;

        delete m_nu_pt;
        delete m_nu_e;
        delete m_nu_px;
        delete m_nu_py;
        delete m_nu_eta;
        delete m_nu_phi;
        delete m_nu_m;
        delete m_nu_pdgid;
        delete m_nu_status;

	delete m_lep_pt;
	delete m_lep_eta;
	delete m_lep_phi;
	delete m_lep_m;
	delete m_lep_pdgid;
	delete m_lep_status;

	delete m_chi0_status;
	delete m_chi0_llp_status;

	delete m_el_pt;
	delete m_el_eta;
	delete m_el_phi;
	delete m_el_m;
	delete m_el_status;

	delete m_mu_pt;
	delete m_mu_eta;
	delete m_mu_phi;
	delete m_mu_m;
	delete m_mu_status;

	delete m_q_pt;
	delete m_q_eta;
	delete m_q_phi;
	delete m_q_m;
	delete m_q_pdgid;
	delete m_q_status;

	delete m_bq_charge;
	delete m_bq_pt;
	delete m_bq_eta;
	delete m_bq_phi;
	delete m_bq_m;
	delete m_bq_status;
	delete m_bq_prodVtx_perp;
	delete m_bq_prodVtx_z;
	delete m_bq_prodVtx_t;
	delete m_bq_prodVtx_eta;
	delete m_bq_prodVtx_phi;
	delete m_bq_prodVtx_nIn;
	delete m_bq_prodVtx_nOut;
	delete m_bq_prodVtx_ID;
	delete m_bq_prodVtx_barcode;

	delete m_j_truthcode;
        delete m_j_pt;
        delete m_j_eta;
        delete m_j_phi;
        delete m_j_status;
        delete m_j_m;

	delete m_bj_truthcode;
	delete m_bj_pt;
	delete m_bj_eta;
	delete m_bj_phi;
    	delete m_bj_status;
	delete m_bj_m;

       	delete m_L1_jet_N_e;
       	delete m_L1_jet_N_et;
        delete m_L1_jet_N_eta;
        delete m_L1_jet_N_phi;
	delete m_L1_delPhi;

        delete m_offline_jet_N_e;
        delete m_offline_jet_N_et;
        delete m_offline_jet_N_eta;
        delete m_offline_jet_N_phi;

}


