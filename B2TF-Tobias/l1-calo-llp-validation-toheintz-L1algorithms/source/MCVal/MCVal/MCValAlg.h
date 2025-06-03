#ifndef MCVal_MCValAlg_H
#define MCVal_MCValAlg_H

#include <AnaAlgorithm/AnaAlgorithm.h>

#include <xAODTruth/TruthParticleContainer.h>
#include <xAODMissingET/MissingETContainer.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include "MCVal/TruthUtils.h"

class MCValAlg : public EL::AnaAlgorithm
{
	public:
		// this is a standard algorithm constructor
		MCValAlg (const std::string& name, ISvcLocator* pSvcLocator);

		~MCValAlg () override;

		// these are the functions inherited from Algorithm
		virtual StatusCode initialize () override;
		virtual StatusCode execute () override;
		virtual StatusCode finalize () override;

	private:
		// Algorithm properties
		int m_pdgIdBSM;
		int m_pdgIdHS;
		int m_pdgIdchi0;

		// Event properties to save
		unsigned int m_runNumber = 0; ///< Run number
		unsigned long long m_eventNumber = 0; ///< Event number

		float m_met;
		float m_met_phi;

		std::vector<float> *m_scalar_pt;
		std::vector<float> *m_scalar_eta;
		std::vector<float> *m_scalar_phi;
		std::vector<float> *m_scalar_m;
		std::vector<int> *m_scalar_pdgid;
		std::vector<int> *m_scalar_status;
		int m_scalar_N;

		float m_L1_MET_Nminus1_met;
		float m_L1_MET_Nminus1_mpt;
		float m_L1_MET_Nminus1_px;
		float m_L1_MET_Nminus1_py;
		float m_L1_MET_Nminus1_pz;
		float m_L1_MET_Nminus1_phi;
		float m_L1_MET_Nminus1_eta;
	
		float m_truth_MET_Nminus1_mpt;
		float m_truth_MET_Nminus1_met;

		std::vector<float> *m_L1_jet_N_e;
		std::vector<float> *m_L1_jet_N_et;
                std::vector<float> *m_L1_jet_N_eta;
                std::vector<float> *m_L1_jet_N_phi;
                int m_L1_jet_N_N;

                std::vector<float> *m_offline_jet_N_e;
                std::vector<float> *m_offline_jet_N_et;
                std::vector<float> *m_offline_jet_N_eta;
                std::vector<float> *m_offline_jet_N_phi;
                int m_offline_jet_N_N;

		float m_L1_leadingJet_N_e; 
                float m_L1_leadingJet_N_et; 
                float m_L1_leadingJet_N_eta; 
                float m_L1_leadingJet_N_phi; 
                float m_L1_leading_delPhi; 
		
                std::vector<float> *m_L1_delPhi;
		
		std::vector<float> *m_higgs_pt;
		std::vector<float> *m_higgs_eta;
		std::vector<float> *m_higgs_phi;
		std::vector<float> *m_higgs_m;
		std::vector<int> *m_higgs_pdgid;
		std::vector<int> *m_higgs_status;
		int m_higgs_N;

		std::vector<float> *m_bsm_pt;
		std::vector<float> *m_bsm_eta;
		std::vector<float> *m_bsm_phi;
		std::vector<float> *m_bsm_m;
		std::vector<int> *m_bsm_pdgid;
		std::vector<int> *m_bsm_status;
		int m_bsm_N;

		Double_t m_scalarReco_pt;
		Double_t m_scalarReco_pz;
		Double_t m_scalarReco_m;

		Double_t m_onTimeHiggsReco_et;
		Double_t m_onTimeHiggsReco_pt;
		Double_t m_onTimeHiggsReco_pz;
		Double_t m_onTimeHiggsReco_eta;
		Double_t m_onTimeHiggsReco_phi;
		Double_t m_onTimeHiggsReco_m;
		Double_t m_onTimeHiggsReco_ProdVtx_t;
		Double_t m_onTimeHiggsReco_ProdVtx_Lxy;
		Double_t m_onTimeHiggsReco_ProdVtx_Lz;

		Double_t m_outOfTimeHiggsReco_et;
		Double_t m_outOfTimeHiggsReco_pt;
		Double_t m_outOfTimeHiggsReco_pz;
		Double_t m_outOfTimeHiggsReco_eta;
		Double_t m_outOfTimeHiggsReco_phi;
		Double_t m_outOfTimeHiggsReco_m;
		Double_t m_outOfTimeHiggsReco_ProdVtx_t;
		Double_t m_outOfTimeHiggsReco_ProdVtx_Lxy;
		Double_t m_outOfTimeHiggsReco_ProdVtx_Lz;

		float m_HiggsReco_del_phi;
		float m_HiggsPseudo_del_phi;
		float m_HiggsPseudoApproach2_del_phi;

		std::vector<float> *m_dR_Bquark_closestBjet;
		std::vector<float> *m_dR_Bjet_closestBquark;
		int m_N_Bquark_closestBjet;
		int m_N_Bjet_closestBquark;
		int m_N_bbbarPairs;

		float m_onTimeHiggsPseudo_et;
		float m_onTimeHiggsPseudo_pt;
		float m_onTimeHiggsPseudo_pz;	
		float m_onTimeHiggsPseudo_eta;
		float m_onTimeHiggsPseudo_phi;
		float m_onTimeHiggsPseudo_m;
		float m_onTimeHiggsPseudo_ProdVtx_t;
		float m_onTimeHiggsPseudo_ProdVtx_Lxy;
		float m_onTimeHiggsPseudo_ProdVtx_Lz;

		float m_outOfTimeHiggsPseudo_et;
		float m_outOfTimeHiggsPseudo_pt;
		float m_outOfTimeHiggsPseudo_pz;	
		float m_outOfTimeHiggsPseudo_eta;
		float m_outOfTimeHiggsPseudo_phi;
		float m_outOfTimeHiggsPseudo_m;
		float m_outOfTimeHiggsPseudo_ProdVtx_t;
		float m_outOfTimeHiggsPseudo_ProdVtx_Lxy;
		float m_outOfTimeHiggsPseudo_ProdVtx_Lz;

		float m_onTimeHiggsPseudoApproach2_et;
		float m_onTimeHiggsPseudoApproach2_pt;
		float m_onTimeHiggsPseudoApproach2_pz;	
		float m_onTimeHiggsPseudoApproach2_eta;
		float m_onTimeHiggsPseudoApproach2_phi;
		float m_onTimeHiggsPseudoApproach2_m;
		float m_onTimeHiggsPseudoApproach2_ProdVtx_t;
		float m_onTimeHiggsPseudoApproach2_ProdVtx_Lxy;
		float m_onTimeHiggsPseudoApproach2_ProdVtx_Lz;

		float m_outOfTimeHiggsPseudoApproach2_et;
		float m_outOfTimeHiggsPseudoApproach2_pt;
		float m_outOfTimeHiggsPseudoApproach2_pz;	
		float m_outOfTimeHiggsPseudoApproach2_eta;
		float m_outOfTimeHiggsPseudoApproach2_phi;
		float m_outOfTimeHiggsPseudoApproach2_m;
		float m_outOfTimeHiggsPseudoApproach2_ProdVtx_t;
		float m_outOfTimeHiggsPseudoApproach2_ProdVtx_Lxy;
		float m_outOfTimeHiggsPseudoApproach2_ProdVtx_Lz;


		std::vector<int> *m_truthParticles_pdgid;

		std::vector<float> *m_nu_pt;
		std::vector<float> *m_nu_px;
		std::vector<float> *m_nu_py;
		std::vector<float> *m_nu_e;
                std::vector<float> *m_nu_eta;
                std::vector<float> *m_nu_phi;
                std::vector<float> *m_nu_m;
                std::vector<int> *m_nu_pdgid;
                std::vector<int> *m_nu_status;
                int m_nu_N;

		std::vector<float> *m_lep_pt;
		std::vector<float> *m_lep_eta;
		std::vector<float> *m_lep_phi;
		std::vector<float> *m_lep_m;
		std::vector<int> *m_lep_pdgid;
		std::vector<int> *m_lep_status;
		int m_lep_N;

		std::vector<float> *m_el_pt;
		std::vector<float> *m_el_eta;
		std::vector<float> *m_el_phi;
		std::vector<float> *m_el_m;
		std::vector<int> *m_el_status;
		int m_el_N;

		std::vector<float> *m_mu_pt;
		std::vector<float> *m_mu_eta;
		std::vector<float> *m_mu_phi;
		std::vector<float> *m_mu_m;
		std::vector<int> *m_mu_status;
		int m_mu_N;

		std::vector<float> *m_q_pt;
		std::vector<float> *m_q_eta;
		std::vector<float> *m_q_phi;
		std::vector<float> *m_q_m;
		std::vector<int> *m_q_pdgid;
		std::vector<int> *m_q_status;
		int m_q_N;

		std::vector<double> *m_bq_charge;
		std::vector<float> *m_bq_pt;
		std::vector<float> *m_bq_eta;
		std::vector<float> *m_bq_phi;
		std::vector<float> *m_bq_m;
		std::vector<int> *m_bq_status;
		int m_bq_N;
		std::vector<float> *m_bq_prodVtx_perp;
		std::vector<float> *m_bq_prodVtx_z;
		std::vector<float> *m_bq_prodVtx_t;
		std::vector<float> *m_bq_prodVtx_eta;
		std::vector<float> *m_bq_prodVtx_phi;
		std::vector<int> *m_bq_prodVtx_nIn;
		std::vector<int> *m_bq_prodVtx_nOut;
		std::vector<int> *m_bq_prodVtx_ID;
		std::vector<int> *m_bq_prodVtx_barcode;

                std::vector<int> *m_j_truthcode;
                std::vector<float> *m_j_pt;
                std::vector<float> *m_j_eta;
                std::vector<float> *m_j_phi;
                std::vector<float> *m_j_m;
                std::vector<int> *m_j_status;

                int m_j_N;

		std::vector<int> *m_bj_truthcode;
		std::vector<float> *m_bj_pt;
		std::vector<float> *m_bj_eta;
		std::vector<float> *m_bj_phi;
		std::vector<float> *m_bj_m;
		std::vector<int> *m_bj_status;

		int m_bj_N;

		std::vector<float> *m_jet_pt;
		std::vector<float> *m_jet_eta;
		std::vector<float> *m_jet_phi;
		std::vector<float> *m_jet_m;
		int m_jet_N;

		std::vector<int> *m_pho_status;
		std::vector<int> *m_pho_llp_status;
		std::vector<int> *m_chi0_status;
		std::vector<int> *m_chi0_llp_status;



		float m_llp1_pt;
		float m_llp1_phi;
		float m_llp1_e;
		float m_llp1_px;
		float m_llp1_py;
		float m_llp1_pz;
		float m_llp1_Lxy;
		float m_llp1_Lz;
		float m_llp1_t;
		float m_llp1_t_eff;
		float m_llp1_ctau;
		float m_llp1_boost;
		float m_llp2_pt;
		float m_llp2_e;
		float m_llp2_phi;
		float m_llp2_px;
		float m_llp2_py;
		float m_llp2_Lxy;
		float m_llp2_Lz;
		float m_llp2_t;
		float m_llp2_t_eff;
		float m_llp2_ctau;
		float m_llp2_boost;
		float m_llp_del_phi;
		float m_llp_del_phi_eff;
		float m_llp_del_R;
		float m_pho_llp1_pt;
		float m_pho_llp1_eta;
		float m_pho_llp1_phi;
		float m_pho_llp1_travelDistance;

		int m_pho_llp_N;
		int m_chi0_llp_N;

		float m_pho_llp1_r;
		float m_pho_llp1_phi_initial;
		float m_pho_llp1_phi_corr;
		float m_pho_llp1_phi_eff;

		float m_pho_llp2_pt;
		float m_pho_llp2_eta;
		float m_pho_llp2_phi;

		float m_pho_llp2_r;
		float m_pho_llp2_phi_initial;
		float m_pho_llp2_phi_corr;
		float m_pho_llp2_phi_eff;
		float m_pho_llp2_travelDistance;


		float m_pho_llp_del_phi;
		float m_pho_llp_del_R;
		float m_pho_llp_del_t;
		float m_chi0_llp1_pt;
		float m_chi0_llp1_e;
		float m_chi0_llp1_px;
		float m_chi0_llp1_py;
		float m_chi0_llp1_eta;
		float m_chi0_llp1_phi;
		float m_chi0_llp2_pt;
		float m_chi0_llp2_e;
		float m_chi0_llp2_px;
		float m_chi0_llp2_py;
		float m_chi0_llp2_pz;
		float m_chi0_llp2_eta;
		float m_chi0_llp2_phi;
		float m_chi0_llp_del_phi;
		float m_chi0_llp_del_R;
		float m_chi0_llp_del_t;

		std::vector<int> *m_onTimeHiggs_bbbar;
		std::vector<int> *m_onTimeHiggs_matchedJets;
		std::vector<int> *m_onTimeHiggs_matchedQuarks;
		std::vector<int> *m_outOfTimeHiggs_bbbar;
		std::vector<int> *m_outOfTimeHiggs_matchedJets;
		std::vector<int> *m_outOfTimeHiggs_matchedQuarks;


		int m_vertex_N;
		std::vector<int> *m_vertex_Nin;
		int m_vertex_N_sIn_chi1Out;
		std::vector<float> *m_vertex_t_sIn_chi1Out;
		int m_vertex_N_chi1In_chi0Out;
		std::vector<float> *m_vertex_t_chi1In_chi0Out;
		std::vector<int> *m_vertex_pdgid_in;
		std::vector<int> *m_vertex_pdgid_out;
		std::vector<int> *m_vertex_Nout;

		const double R_CALO = 2000; //inner radius of Calorimeter

};

#endif
