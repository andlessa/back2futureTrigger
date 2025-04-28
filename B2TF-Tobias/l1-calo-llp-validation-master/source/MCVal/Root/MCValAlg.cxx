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



StatusCode MCValAlg :: initialize ()
{

	ATH_MSG_ALWAYS("THE PDGID OF THE SCALAR IS " << m_pdgIdBSM);
	ATH_MSG_ALWAYS("THE PDGID OF THE CHI1 IS " << m_pdgIdHS);
	ATH_MSG_ALWAYS("THE PDGID OF THE CHI0 IS " << m_pdgIdchi0);

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

	//Higgs Boson (most likely empty, as we're expecting the higgs to be off-shell)
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
	m_onTimeHiggsReco_et = -999.;
	m_onTimeHiggsReco_pt = -999.;
	m_onTimeHiggsReco_pz = -999.;
	m_onTimeHiggsReco_eta = -999.;
	m_onTimeHiggsReco_phi = -999.;
	m_onTimeHiggsReco_m = -999.;
	m_onTimeHiggsReco_ProdVtx_t = -999.;
	m_onTimeHiggsReco_ProdVtx_Lxy = -999.;	
	m_onTimeHiggsReco_ProdVtx_Lz = -999.;

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
	m_outOfTimeHiggsReco_et = -999.;
	m_outOfTimeHiggsReco_pt = -999.;
	m_outOfTimeHiggsReco_pz = -999.;
	m_outOfTimeHiggsReco_eta = -999.;
	m_outOfTimeHiggsReco_phi = -999.;
	m_outOfTimeHiggsReco_m = -999.;
	m_outOfTimeHiggsReco_ProdVtx_t = -999.;
	m_outOfTimeHiggsReco_ProdVtx_Lxy = -999.;
        m_outOfTimeHiggsReco_ProdVtx_Lz = -999.;

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
	m_HiggsReco_del_phi = -999.;
	mytree->Branch ("HiggsReco_del_phi" , &m_HiggsReco_del_phi);	
	
	//clustered Higgs bosons from b jets
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
	//FIXME declaration of member variables to -999.?
	mytree->Branch("nonpromptLLPpt",&m_llp1_pt);
	mytree->Branch("nonpromptLLPLxy",&m_llp1_Lxy);
	mytree->Branch("nonpromptLLPLz",&m_llp1_Lz);
	mytree->Branch("nonpromptLLPt",&m_llp1_t);
	mytree->Branch("nonpromptLLPt_eff",&m_llp1_t_eff);
	mytree->Branch("nonpromptLLPctau",&m_llp1_ctau);
	mytree->Branch("nonpromptLLPboost",&m_llp1_boost);

	//on-time Chi1
	mytree->Branch("promptLLPpt",&m_llp2_pt);
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
	mytree->Branch("chi0_llp1_eta",&m_chi0_llp1_eta);
	mytree->Branch("chi0_llp1_phi", &m_chi0_llp1_phi);

	//on-time-chi1-matched chi0
	mytree->Branch("chi0_llp2_pt", &m_chi0_llp2_pt);
	mytree->Branch("chi0_llp2_eta", &m_chi0_llp2_eta);
	mytree->Branch("chi0_llp2_phi", &m_chi0_llp2_phi);


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
bool compareT(const xAOD::TruthParticle* part1, const xAOD::TruthParticle* part2)
{
	if(part1->hasDecayVtx() && part2-> hasDecayVtx())		return (part1->decayVtx()->t()>part2->decayVtx()->t());
	else return true;
}

//same fct as above, one of them should be obsolete
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
	const xAOD::JetContainer* truth_largeR_jets = nullptr;
	ANA_CHECK ( evtStore()->retrieve(truth_largeR_jets, "AntiKt10TruthTrimmedPtFrac5SmallR20Jets" ) );
	
	//truth vertices >> currently not used, though might be useful
	const xAOD::TruthVertexContainer *truth_vertices = nullptr;
	ANA_CHECK ( evtStore()->retrieve(truth_vertices, "TruthVertices") );

	//create vectors to store the different particle types
	std::vector<const xAOD::TruthParticle*> vec_scalar;
	std::vector<const xAOD::TruthParticle*> vec_higgs;
	std::vector<const xAOD::TruthParticle*> vec_llp;
	std::vector<const xAOD::TruthParticle*> vec_chi0;
	std::vector<const xAOD::TruthParticle*> vec_chi0_llp;
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

		//Electrons
		if(Truth::isFinalElectron(tp)) {
			vec_el.push_back(tp);
			vec_lep.push_back(tp);
		}

		//Muons
		//if(Truth::isFinalMuon(tp) && !Truth::isFromPhoton(tp)) {
		//if(Truth::isFinalMuon(tp) && !Truth::isFromPhoton(tp) && Truth::isFromParticle(tp,m_pdgIdBSM)) {
		if(Truth::isFinalMuon(tp)) {
			vec_mu.push_back(tp);
			vec_lep.push_back(tp);
		}

		//Quarks
		//if(Truth::isFinalQuark(tp) && !Truth::isFromGluon(tp)) {
		//if(Truth::isFinalQuark(tp) && !Truth::isFromGluon(tp) && Truth::isFromParticle(tp,m_pdgIdBSM)) {
		if(Truth::isFinalQuark(tp) && !Truth::isFromGluon(tp)) {
			vec_q.push_back(tp);
			if(std::abs(tp->pdgId())==5)		vec_bquark.push_back(tp);
		}

	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//sort all vectors which are related to on-time and delayed chi1 according to their time, first entry is largest time (i.e. most delayed particle)
	sort(vec_llp.begin(),vec_llp.end(),compareTdecay); //based on their decay vertex
	//sort(vec_llp.begin(),vec_llp.end(),comparePt); //based on their decay vertex
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
	m_llp1_Lxy=-999;
	m_llp1_Lz=-999999;
	m_llp1_t=-999;
	m_llp1_t_eff=-999;
	m_llp1_boost=-999;
	m_llp1_ctau=-999;

	m_llp2_pt=-999;
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
	m_chi0_llp1_eta =-999;
	m_chi0_llp1_phi =-999;
	m_chi0_llp2_pt =-999;
	m_chi0_llp2_eta =-999;
	m_chi0_llp2_phi =-999;
	m_chi0_llp_del_phi =-999;
	m_chi0_llp_del_R =-999;
	m_chi0_status->clear();
	m_chi0_llp_status->clear();

	//FIXME Maybe store chi0 variables in nTuples?	
	// Store number of chi0
	//hist("h_chi0_N")->Fill(vec_chi0.size());
	//hist("h_chi0_llp_N")->Fill(vec_chi0_llp.size());

	// Loop over chi0s
	for(auto chi0 : vec_chi0) {
		// Fill histograms
		//hist("h_chi0_pt")->Fill(chi0->pt()/Truth::GeV);
		//hist("h_chi0_eta")->Fill(chi0->eta());
		//hist("h_chi0_phi")->Fill(chi0->phi());
		m_chi0_status->push_back(chi0->status());
		//hist("h_chi0_status")->Fill(chi0->status());
	}

	// Loop over truth matched chi0s
	for(auto chi0_llp : vec_chi0_llp) {
		// Fill histograms and tuples
		//hist("h_chi0_llp_pt")->Fill(chi0_llp->pt()/Truth::GeV);
		//hist("h_chi0_llp_eta")->Fill(chi0_llp->eta());
		//hist("h_chi0_llp_phi")->Fill(chi0_llp->phi());
		m_chi0_llp_status->push_back(chi0_llp->status());
		//hist("h_chi0_llp_status")->Fill(chi0_llp->status());
		//const xAOD::TruthVertex_v1* decVtx = chi0_llp->prodVtx();
		//hist("chi0Lxy")->Fill(decVtx->perp());
		//hist("chi0Lz")->Fill(decVtx->z());
	}

	m_chi0_llp_N = vec_chi0_llp.size();
	if(vec_chi0_llp.size() == 2){ 					//change from "==2" to ">=2" makes this more inclusive
		auto v1 = vec_chi0_llp.at(0)->p4();
		auto v2 = vec_chi0_llp.at(1)->p4();
		m_chi0_llp1_pt = v1.Pt()/Truth::GeV;
		m_chi0_llp1_eta = v1.Eta();
		m_chi0_llp1_phi = v1.Phi();
		m_chi0_llp2_pt = v2.Pt()/Truth::GeV;
		m_chi0_llp2_eta = v2.Eta();
		m_chi0_llp2_phi = v2.Phi();
		m_chi0_llp_del_phi = std::abs(v1.DeltaPhi(v2));
		m_chi0_llp_del_R = v1.DrEtaPhi(v2);
	}
	else			ATH_MSG_WARNING("This event contains " << vec_chi0_llp.size() << " BSM Chi0 particles. This is NOT expected.");


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


//HIGGS
//Method 1. on shell Higgs

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
                        //m_onTimeHiggsReco_ProdVtx_t = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->t(); //FIXME: Check if this works
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
                        //m_outOfTimeHiggsReco_ProdVtx_t = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->t(); //FIXME: Check if this works
                        //m_outOfTimeHiggsReco_ProdVtx_Lxy = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->perp(); //FIXME: Check if this works
                        //m_outOfTimeHiggsReco_ProdVtx_Lz = vec_RecoHiggs.at(onTimeRECO)->prodVtx()->z(); //FIXME: Check if this works
                }
        	if (	idxOnTime >= 0 && vec_HiggsRECO.at(idxOnTime) != nullptr &&
			idxOutOfTime >= 0 && vec_HiggsRECO.at(idxOutOfTime) != nullptr
		){
                	m_HiggsReco_del_phi = std::abs(vec_HiggsRECO.at(idxOnTime)->DeltaPhi(*vec_HiggsRECO.at(idxOutOfTime)));
        	}
		if (idxOnTime >= 0 && vec_HiggsCLUSTER.at(idxOnTime) != nullptr){
			ATH_MSG_ALWAYS("WITHIN OT");
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
                        ATH_MSG_ALWAYS("WITHIN OOT");
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

	
//Delta Phi correction for out-of-calorimeter decays (previously applied)

	//make corrections for out-of calorimeter decays, for now test ignorring size requirements
	if(promptLLP != 0 && nonpromptLLP != 0 && nonpromptLLP->hasDecayVtx() && promptLLP->hasDecayVtx()){
		
		auto llpv1 = nonpromptLLP->decayVtx()->v4();
		auto llpv2 = promptLLP->decayVtx()->v4();
		m_llp_del_phi = std::abs(llpv1.DeltaPhi(llpv2));
		m_llp_del_R = std::abs(llpv1.DeltaR(llpv2));

		//We calculate "effective" phi and timing values of the Calorimeter Cluster for the case that the chi1 particle decays before entering the Calorimeter
		//To do so, we propagate the photon phi and timing information to "effective" phi and timing values for the CALO cluster
		//The corrections are calculated in the transverse plane, solving the geometric situation 
		//We make this correction separately for both chi1's
		//FIXME Change names of pho1 and pho2 to higgs1 and higgs2 respectively	
		//Correction for non-prompt chi1:
		TLorentzVector *pho1 = nullptr;
		TLorentzVector *pho2 = nullptr;
		bool skip = true;
                if (    idxOnTime >= 0 && vec_HiggsCLUSTER.at(idxOnTime) != nullptr &&
                        idxOutOfTime >= 0 && vec_HiggsCLUSTER.at(idxOutOfTime) != nullptr
                ){
			pho1 = vec_HiggsCLUSTER.at(idxOnTime);
			pho2 = vec_HiggsCLUSTER.at(idxOutOfTime);
			skip = false;
		}			
		else		ATH_MSG_WARNING("At least one of the Pseudo Higgs Candidates is not accessible.");

	
		//Calculate (and store) the decay position
		float llpv1_Lxy = llpv1.Perp();
		//Initial value for the effective time is the decay time of the chi1
		m_llp1_t_eff = m_llp1_t*3e2; //in mm/c
		//In the case the decay happens before entering the CALO, we make some corrections:
	
		if (llpv1_Lxy < R_CALO && skip == false){
			//We need the phi of the chi1:
			double phi = nonpromptLLP->p4().Phi();
			//The photon might have a different angular orientation, hence we calculate Delta phi between chi1 and the photon:
			double delta = -nonpromptLLP->p4().DeltaPhi(*pho1); //negative sign due to definition of TLorentzVector::DeltaPhi
			m_pho_llp1_phi_corr = delta;
			//The effective angle depends on the decay position, the inner radius of the calorimeter, the phi of the chi1 and the delta phi between chi1 and the photon
			//For the calculation we use the sine theorem
			double abs_correction = abs(delta) - asin(sin(abs(delta)) * llpv1_Lxy / R_CALO); 
			double sgn_correction = delta/abs(delta);
			double phi_eff = phi + delta/abs(delta) * abs_correction; //NOT in range (-pi,pi)... need to be converted:
			while (phi_eff > 3.142) {
				phi_eff -= 2 * 3.142;
			}
			while (phi_eff <= -3.142) {
				phi_eff += 2 * 3.142;
			}

			//Modify the phi value of the decay (as we need this value to calculate the effective delta phi between both chi1's)
			llpv1.SetPhi(phi_eff);

			//Timing correction: FIXME implement timing correction
			//We need the time of the chi1 decay:
			double t = nonpromptLLP->decayVtx()->t();
			//we correct for the photon travel and the calibration time
			//The photon travel distance s can be calculated by the cosine theorem, i.e. s depends on the decay Vaetex position, the calo radius and the correction angle
			double s = sqrt(llpv1_Lxy * llpv1_Lxy + R_CALO * R_CALO - 2 * llpv1_Lxy * R_CALO * cos(abs_correction)); //using cosine theorem
			m_pho_llp1_travelDistance = s; //in mm
			//The effective time is the decay time of the chi1, summed with the photon travel time. The calibration time is corrected for by subtracting the CALO inner radius.
			double t_eff = t + s - R_CALO; //in mm/c
			m_llp1_t_eff = t_eff; //conversion to ns
		}

		//Analouge corrections are performed for the prompt chi1
		float llpv2_Lxy = llpv2.Perp();
		m_llp2_t_eff = m_llp2_t*3e2; //convert to mm/c

		if (llpv2_Lxy < R_CALO && skip == false){
			double phi = promptLLP->p4().Phi(); 
			double delta = -promptLLP->p4().DeltaPhi(*pho2); //negative sign due to definition of TLorentzVector::DeltaPhi
			m_pho_llp2_phi_corr = delta;
			double abs_correction = abs(delta) - asin(sin(abs(delta)) * llpv2_Lxy / R_CALO); 
			double sgn_correction = delta/abs(delta);
			double phi_eff = phi + delta/abs(delta) * abs_correction; //NOT in range (-pi,pi)... need to be converted:
			while (phi_eff > 3.142) {
				phi_eff -= 2 * 3.142;
			}
			while (phi_eff <= -3.142) {
				phi_eff += 2 * 3.142;
			}

			llpv2.SetPhi(phi_eff);

			double t = promptLLP->decayVtx()->t();
			double s = sqrt(llpv2_Lxy * llpv2_Lxy + R_CALO * R_CALO - 2 * llpv2_Lxy * R_CALO * cos(abs_correction)); //using cosine theorem
			m_pho_llp2_travelDistance = s;
			
			//hist("h_pho_llp2_travelDistance")->Fill(s);
			double t_eff = t + s - R_CALO;
			m_llp2_t_eff = t_eff;
		}
		
		//The effective delta phi is calculated based on the effective phi values of both chi1:
		m_llp_del_phi_eff = std::abs(llpv1.DeltaPhi(llpv2));
	
	}

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
}


