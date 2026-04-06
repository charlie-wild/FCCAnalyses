#ifndef MYUTILS_ANALYZERS_H
#define MYUTILS_ANALYZERS_H
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/VertexData.h"

#include "TLorentzVector.h"
#include "VertexingUtils.h"
#include <vector>

namespace FCCAnalyses{

namespace myUtils{

  struct FCCAnalysesComposite{
    TLorentzVector particle;
    ROOT::VecOps::RVec<int> index;//index in the RP
    edm4hep::VertexData vertex;
    int charge;
    int mc_index;
  };


  struct FCCAnalysesComposite2{
    TLorentzVector particle;
    int vertex; //index in the vertex collection
    int charge;
    int mc_index;
  };

  struct filter_PV{
    filter_PV(bool arg_pv);
    bool m_pv=true;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
								      ROOT::VecOps::RVec<int> index);
  };

  struct build_D0 {
    build_D0(float arg_mass, float arg_p, bool arg_filterPV);
    float m_mass=0.05;
    float m_p=1.;
    bool m_filterPV=true;
    ROOT::VecOps::RVec<FCCAnalysesComposite> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							 ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							 ROOT::VecOps::RVec<int> pions,
							 ROOT::VecOps::RVec<int> kaons,
							 ROOT::VecOps::RVec<int> pvindex);
  };



  struct build_tau23pi {
    build_tau23pi( float arg_masslow, float arg_masshigh, float arg_p, float arg_angle, bool arg_rho);
    float m_masslow=0.05;
    float m_masshigh=3.0;
    float m_p=1.;
    float m_angle=1.;
    bool m_rho = true;
    ROOT::VecOps::RVec<FCCAnalysesComposite2> operator() (ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);
  };



  struct sel_PV {
    sel_PV(bool arg_closest);
    bool m_closest;
    VertexingUtils::FCCAnalysesVertex operator()(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> pv);
  };


  struct sel_PID {
    sel_PID(int arg_PDG);
    int m_PDG=211;
    ROOT::VecOps::RVec<int> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);
  };


  ROOT::VecOps::RVec<edm4hep::TrackState> get_pseudotrack(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);


  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_tau23pi(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_B2Kstee(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_B2Kstmumu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_Bd2KstNuNu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							     ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_Bs2PhiNuNu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							     ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_Bd2MuMu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> build_D02KK(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
								   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop, bool mass_cut=true);

  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Pi02photonphoton(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop, bool mass_cut=true, bool momentum_cut=true);

  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Pi02photonphoton_wider(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Pi02photonphoton_cut_3(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Pi02photonphoton_cut_5(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Pi02photonphoton_cut_7(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<FCCAnalysesComposite> build_B02D0Pi0(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                  ROOT::VecOps::RVec<FCCAnalysesComposite2> D0,
                  ROOT::VecOps::RVec<FCCAnalysesComposite> Pi0,
                  ROOT::VecOps::RVec<float>& thrust,
                  bool mass_cut=false,
                  bool momentum_cut=true);

  ROOT::VecOps::RVec<int> get_daughter_index_from_B0Candidate(ROOT::VecOps::RVec<FCCAnalysesComposite>& B0Candidate, int daughter_index);

  ROOT::VecOps::RVec<int> variable_access_B0_Daughters_int(ROOT::VecOps::RVec<int>& index_of_daughter, ROOT::VecOps::RVec<int>& daughter_integer_variable);

  ROOT::VecOps::RVec<float> variable_access_B0_Daughters_flt(ROOT::VecOps::RVec<int>& index_of_daughter, ROOT::VecOps::RVec<float>& daughter_integer_variable);

  ROOT::VecOps::RVec<int> indRP_matched_to_list( ROOT::VecOps::RVec<int>&  mcParticles_indices,
						 ROOT::VecOps::RVec<int>& recind,
						 ROOT::VecOps::RVec<int>& mcind,
						 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
						 ROOT::VecOps::RVec<edm4hep::MCParticleData>& mc);

  bool array_searcher(int index, ROOT::VecOps::RVec<int> index_array);

  ROOT::VecOps::RVec<int> truthmatch_D0(ROOT::VecOps::RVec<FCCAnalysesComposite2>& D0, ROOT::VecOps::RVec<int>& RP_truthmatched_indices, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& vertex);

  ROOT::VecOps::RVec<int> truthmatch_Pi0(ROOT::VecOps::RVec<FCCAnalysesComposite>& Pi0, ROOT::VecOps::RVec<int>& RP_truthmatched_indices);

  ROOT::VecOps::RVec<int> truthmatch_B0(ROOT::VecOps::RVec<FCCAnalysesComposite>& B0, ROOT::VecOps::RVec<int>& RP_truthmatched_indices); 

  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> get_MCVertexObject(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
									     ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> merge_VertexObjet(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> in);

  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> get_VertexObject(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver,
									 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
									 ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
									 ROOT::VecOps::RVec<int> recin,
									 ROOT::VecOps::RVec<int> mcin);

  ROOT::VecOps::RVec<float> get_Vertex_mass(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
					    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco);

  ROOT::VecOps::RVec<float> get_Vertex_x(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_y(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_z(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_xErr(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_yErr(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_zErr(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_chi2(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<int> get_Vertex_isPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<int> get_Vertex_ntracks(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_d2PV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
					    int comp);

  ROOT::VecOps::RVec<float> get_Vertex_d2PVError(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						 int comp);

  ROOT::VecOps::RVec<int> get_Vertex_indMC(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
					   ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver);

  ROOT::VecOps::RVec<int> get_Vertex_indMC(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_Vertex_d2MC(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
					    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver,
					    ROOT::VecOps::RVec<int> mcind,
					    int comp);

  std::vector<std::vector<int>> get_Vertex_ind(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  float get_d0(TVector3 x, TVector3 p);
  float get_z0(TVector3 x, TVector3 p);


  ROOT::VecOps::RVec<TVector3> get_MCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex);

  ROOT::VecOps::RVec<float> get_MCVertex_x(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex);

  ROOT::VecOps::RVec<float> get_MCVertex_y(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex);

  ROOT::VecOps::RVec<float> get_MCVertex_z(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex);

  ROOT::VecOps::RVec<int> get_NTracksMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex);

  std::vector< std::vector<int> > get_MCindMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex);

  std::vector<std::vector<int>> get_MCpdgMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
						  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

  std::vector<std::vector<int>> get_MCpdgMotherMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
							ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

  std::vector<std::vector<int>> get_MCpdgGMotherMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
							 ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

  ROOT::VecOps::RVec<int> get_MCMother1(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<int> get_MCMother2(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<int> get_MCDaughter1(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					  ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<int> get_MCDaughter2(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					  ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<int> get_MCMother2(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					ROOT::VecOps::RVec<int> ind);


  ROOT::VecOps::RVec<FCCAnalysesComposite> add_truthmatched(ROOT::VecOps::RVec<FCCAnalysesComposite> comp,
							    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
							    ROOT::VecOps::RVec<int> rp2mc,
							    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							    ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<FCCAnalysesComposite2> add_truthmatched2(ROOT::VecOps::RVec<FCCAnalysesComposite2> comp,
							      ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
							      ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							      ROOT::VecOps::RVec<int> rp2mc,
							      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							      ROOT::VecOps::RVec<int> ind);


  ROOT::VecOps::RVec<int> get_trueVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
					 ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					 ROOT::VecOps::RVec<int> ind,
					 int mother,
					 int grandmother);

  
  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_vtxind(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_vtxFDinfo(ROOT::VecOps::RVec<FCCAnalysesComposite2> in, ROOT::VecOps::RVec<float> FD_info);


  bool isPV(edm4hep::ReconstructedParticleData recop,
	    ROOT::VecOps::RVec<int> pvindex);

  ROOT::VecOps::RVec<int> getMC_daughter(int daughterindex,
					 ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
					 ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<int> getMC_parent(int parentindex,
				       ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
				       ROOT::VecOps::RVec<int> ind);

  int getMC_parent(int parentindex,
		   edm4hep::MCParticleData in,
		   ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<float> get_flightDistanceVertex(ROOT::VecOps::RVec<FCCAnalysesComposite> in,
						     VertexingUtils::FCCAnalysesVertex pv);

  ROOT::VecOps::RVec<float> get_flightDistanceVertex(ROOT::VecOps::RVec<FCCAnalysesComposite> in,
						     edm4hep::VertexData pv);

  float get_distanceVertex(edm4hep::VertexData v1,
			   edm4hep::VertexData v2,
			   int comp);

  float get_distanceErrorVertex(edm4hep::VertexData v1,
				edm4hep::VertexData v2,
				int comp);

  float get_distance(TVector3 v1,
		     TVector3 v2);

  float get_distance(edm4hep::Vector3f v1,
		     TVector3 v2,
		     int comp);



  ROOT::VecOps::RVec<int> get_compmc(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<TLorentzVector> getFCCAnalysesComposite_particle(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> getFCCAnalysesComposite_index(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<edm4hep::VertexData> getFCCAnalysesComposite_vertex(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  int getFCCAnalysesComposite_N(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite> in);


  int getFCCAnalysesComposite_N(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
							 ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
						    int index);

  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_type(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
						    int index);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_p(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						      int type);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_p(ROOT::VecOps::RVec<FCCAnalysesComposite> in,
						      int type);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_B(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						      ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_p(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						      ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
						      int index,
						      int type);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_E(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_E(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);  

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Pt(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Pt(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_charge(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_vertex(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_vtx_position(ROOT::VecOps::RVec<FCCAnalysesComposite2> in, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> verteces, int i);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_vtx_position_err(ROOT::VecOps::RVec<FCCAnalysesComposite2> in, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> verteces, int i);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_vtx_chi2(ROOT::VecOps::RVec<FCCAnalysesComposite2> in, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> verteces);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Phi(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Phi(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Eta(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Eta(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Theta(ROOT::VecOps::RVec<FCCAnalysesComposite> in);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_Theta(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);

  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_mcvertex(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
							   ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_q(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
						    int index);


  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_d0(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						       ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						       int index);

  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_z0(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
						       ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
						       int index);

  ROOT::VecOps::RVec<edm4hep::TrackState> getFCCAnalysesComposite_track(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
									ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_trackd0(ROOT::VecOps::RVec<edm4hep::TrackState> in);
  ROOT::VecOps::RVec<float> get_trackz0(ROOT::VecOps::RVec<edm4hep::TrackState> in);

  ROOT::VecOps::RVec<int> getFCCAnalysesComposite_truthMatch(ROOT::VecOps::RVec<FCCAnalysesComposite2> in);


  ROOT::VecOps::RVec<FCCAnalysesComposite> build_Bu2D0Pi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							 ROOT::VecOps::RVec<FCCAnalysesComposite> d0,
							 ROOT::VecOps::RVec<int> pions);


  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> PID(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							     ROOT::VecOps::RVec<int> recind,
							     ROOT::VecOps::RVec<int> mcind,
							     ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> get_RP_atVertex(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
									 ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);



  float build_invmass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
		      ROOT::VecOps::RVec<int> index);

  TLorentzVector build_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
			   ROOT::VecOps::RVec<int> index);

  ROOT::VecOps::RVec<float> get_Vertex_thrusthemis_angle(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							 ROOT::VecOps::RVec<float> thrust);

  ROOT::VecOps::RVec<float> get_DVertex_thrusthemis_angle(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
							  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
							  ROOT::VecOps::RVec<float> thrust);

  ROOT::VecOps::RVec<int> get_Vertex_thrusthemis_emin(ROOT::VecOps::RVec<float> angle,
						      float eneg,
						      float epos);

  ///index ==1 -> positive angle == minimum energy
  ///index ==0 -> negative angle == maximum energy
  ROOT::VecOps::RVec<int> get_Vertex_thrusthemis(ROOT::VecOps::RVec<float> angle,
						 int index);


  ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>> build_rho(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
										       ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
										       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop);

  ROOT::VecOps::RVec<float> get_mass(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>> in,
				     int index);

  ROOT::VecOps::RVec<float> get_px(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>> in,
				   int index);

  ROOT::VecOps::RVec<float> get_py(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>> in,
				   int index);

  ROOT::VecOps::RVec<float> get_pz(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>> in,
				   int index);

  /// return the momenta of the input ReconstructedParticles
  float get_p(edm4hep::ReconstructedParticleData in);
  float get_E(edm4hep::ReconstructedParticleData in);

  ROOT::VecOps::RVec<edm4hep::TrackState> get_truetrack(ROOT::VecOps::RVec<int> in,
							ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
							ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

  int get_PV_ntracks(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);
  int hasPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  float get_dPV2DV_min(ROOT::VecOps::RVec<float> in);
  float get_dPV2DV_max(ROOT::VecOps::RVec<float> in);
  float get_dPV2DV_ave(ROOT::VecOps::RVec<float> in);

  int get_Npos(ROOT::VecOps::RVec<float> in);
  int get_Nneg(ROOT::VecOps::RVec<float> in);

  ROOT::VecOps::RVec<edm4hep::MCParticleData> build_truerho(ROOT::VecOps::RVec<int> vertexind,
							    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
							    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);


  ROOT::VecOps::RVec<float> getFCCAnalysesComposite_anglethrust(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
								ROOT::VecOps::RVec<float> thrust);

  int has_anglethrust_emin(ROOT::VecOps::RVec<float> angle);

  ROOT::VecOps::RVec<int> D0_and_Pi0_same_hemisphere(ROOT::VecOps::RVec<FCCAnalysesComposite>& B0_candidates, ROOT::VecOps::RVec<FCCAnalysesComposite2>& D0_candidates,
								ROOT::VecOps::RVec<FCCAnalysesComposite>& Pi0_candidates, ROOT::VecOps::RVec<float>& thrust);

  int D0_and_Pi0_same_hemisphere(FCCAnalysesComposite& B0_candidates, FCCAnalysesComposite2& D0,
								FCCAnalysesComposite& Pi0, ROOT::VecOps::RVec<float>& thrust);

  TVector3 get_MC_PV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver,
                ROOT::VecOps::RVec<int> mcind);

  float get_B0D0_FD_x(TVector3 PV, edm4hep::MCParticleData K);

  float get_B0D0_FD_y(TVector3 PV, edm4hep::MCParticleData K);

  float get_B0D0_FD_z(TVector3 PV, edm4hep::MCParticleData K);

  float get_B0D0_FD(TVector3 PV, edm4hep::MCParticleData K);

  VertexingUtils::FCCAnalysesVertex get_PV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex);

  ROOT::VecOps::RVec<float> get_IP(VertexingUtils::FCCAnalysesVertex PV, ROOT::VecOps::RVec<FCCAnalysesComposite2> D0, ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertices);

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> get_D0Tracks(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& vertex, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recop, 
    ROOT::VecOps::RVec<FCCAnalysesComposite2>& D0Candidates, ROOT::VecOps::RVec<edm4hep::TrackState>& tracks);

  ROOT::VecOps::RVec<edm4hep::TrackState> get_D0Tracks(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& vertex, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recop,
    FCCAnalysesComposite2& D0, ROOT::VecOps::RVec<edm4hep::TrackState>& tracks);

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>  ReconstructedD0( ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>&  tracks_vectors, bool AddMassConstraints);

  ROOT::VecOps::RVec<edm4hep::TrackState>  ReconstructedD0( ROOT::VecOps::RVec<edm4hep::TrackState>&  tracks_vectors, bool AddMassConstraints);

  ROOT::VecOps::RVec<TVector3> Momentum_ReconstructedD0( ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& D0_pseudoTracks, bool is_neutral=true);

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> get_photon_fake_tracks_from_Pi0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recop,ROOT::VecOps::RVec<FCCAnalysesComposite>& Pi0Candidates);

  ROOT::VecOps::RVec<float> component_access( ROOT::VecOps::RVec<TVector3> TVectors , int component );

  ROOT::VecOps::RVec<edm4hep::TrackState> get_photon_fake_tracks_from_Pi0(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recop,FCCAnalysesComposite& Pi0);

  ROOT::VecOps::RVec<edm4hep::TrackState>  tracks_for_fitting_the_B0_vertex(ROOT::VecOps::RVec<edm4hep::TrackState>  ReconstructedD0,ROOT::VecOps::RVec<edm4hep::TrackState> BachelorPi0);

  ROOT::VecOps::RVec<edm4hep::TrackState>  tracks_for_fitting_the_B0_vertex_pi0(ROOT::VecOps::RVec<edm4hep::TrackState>  ReconstructedD0,ROOT::VecOps::RVec<edm4hep::TrackState> BachelorPi0);

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> B0_vertex_fit(ROOT::VecOps::RVec<FCCAnalysesComposite>& B0_candidates, ROOT::VecOps::RVec<FCCAnalysesComposite>& Pi0_candidates, ROOT::VecOps::RVec<FCCAnalysesComposite2>& D0_candidates,
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& vertex, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recopatvtx, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recop, ROOT::VecOps::RVec<edm4hep::TrackState>& tracks);

  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> B0_Vertex_fitter(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& B0_tracks, int num_tracks);

  ROOT::VecOps::RVec<int> B0_Vertex_n_tracks(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& B0_tracks);

  ROOT::VecOps::RVec<edm4hep::VertexData> B0_Vertex_getter(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& B0_vertices);

  ROOT::VecOps::RVec<float> B0_Vertex_position(ROOT::VecOps::RVec<edm4hep::VertexData>& B0_vertices, int componant);

  ROOT::VecOps::RVec<float> B0_Vertex_chi2(ROOT::VecOps::RVec<edm4hep::VertexData>& B0_vertices);

  ROOT::VecOps::RVec<edm4hep::TrackState>  ReconstructedPi0( ROOT::VecOps::RVec<edm4hep::TrackState>&  tracks, bool AddMassConstraints);

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> B0_vertex_fit_pi0(ROOT::VecOps::RVec<FCCAnalysesComposite>& B0_candidates, ROOT::VecOps::RVec<FCCAnalysesComposite>& Pi0_candidates, ROOT::VecOps::RVec<FCCAnalysesComposite2>& D0_candidates,
        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>& vertex, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recopatvtx, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recop, ROOT::VecOps::RVec<edm4hep::TrackState>& tracks);

  float get_B0_MCVertex_x(edm4hep::MCParticleData D0);

  float get_B0_MCVertex_y(edm4hep::MCParticleData D0);

  float get_B0_MCVertex_z(edm4hep::MCParticleData D0);

  float get_B0_MCVertex(edm4hep::MCParticleData D0);

  ROOT::VecOps::RVec<float> get_opening_angle(ROOT::VecOps::RVec<float> Px1, ROOT::VecOps::RVec<float> Px2, ROOT::VecOps::RVec<float>Py1, ROOT::VecOps::RVec<float> Py2, ROOT::VecOps::RVec<float>Pz1, ROOT::VecOps::RVec<float> Pz2, ROOT::VecOps::RVec<float> P1, ROOT::VecOps::RVec<float> P2);

}//end NS myUtils

}//end NS FCCAnalyses
#endif