#include "TreeReader.h"  // call libraries from ROOT and C++

using namespace fastjet;

void JetTreesRecluster(TString InputFileList, TString OutputFile, bool removeelectrons, bool donhitcut){

	typedef ROOT::Math::PxPyPzEVector LorentzVector;
	
    // Define R values
    std::vector<float> R_values;
    for (int i = 1; i <= 10; i++) R_values.push_back(i * 0.1);

    double minCstPt            = 0.2 ;				 // minimum pT of objects
    double maxCstPt            = 100.;  			 // maximum pT of objects
    double minJetPt            = 1.0 ;  			 // minimum jet pT
    double ghostMaxRap         = 3.5;   			 // maximum rapidity of ghosts
    double ghostArea           = 0.01; 				 // area per ghost
    int numGhostRepeat         = 1;                  // reuse count of ghosts


	// Read the list of input file(s)
	fstream FileList;
	FileList.open(Form("%s",InputFileList.Data()), ios::in);
	if(!FileList.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << InputFileList.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> FileListVector;
	string FileChain;
	while(getline(FileList, FileChain)){FileListVector.push_back(FileChain.c_str());}
	FileList.close();	
	TChain *mychain = new TChain("events");
	for (std::vector<TString>::iterator listIterator = FileListVector.begin(); listIterator != FileListVector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			mychain->Add(*listIterator);
		}else{cout << "File: " << *listIterator << " failed!" << endl;}
	}

	// Reading trees
	std::unique_ptr<TTreeReader> tree_reader;
	TreeReader(mychain, tree_reader); 
	
	// Make Output
	TFile *OutFile = TFile::Open(Form("%s",OutputFile.Data()),"RECREATE");	

	int NEVENTS = 0;
	// Reco Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> RecoJet_pt;
	std::vector<float> RecoJet_eta;
	std::vector<float> RecoJet_phi;
	std::vector<float> RecoJet_E;
	std::vector<float> RecoJet_M;
	std::vector<bool> RecoJet_hasElectron;
	std::vector<float> RecoJet_maxPtPart_pt; 
	// Reco jet constituents
	std::vector<std::vector<float>> RecoJet_constituent_pt;
	std::vector<std::vector<float>> RecoJet_constituent_eta;
	std::vector<std::vector<float>> RecoJet_constituent_phi;
	std::vector<std::vector<int>> RecoJet_constituent_nhits;

	// Gen Jets (Variable-length vectors for multiple jets per event)
	std::vector<float> GenJet_pt;
	std::vector<float> GenJet_eta;
	std::vector<float> GenJet_phi;
	std::vector<float> GenJet_E;
	std::vector<float> GenJet_M;
	std::vector<bool> GenJet_hasElectron;
	std::vector<bool> GenJet_hasNeutral;
	std::vector<float> GenJet_maxPtPart_pt;
	// Gen jet constituents
	std::vector<std::vector<float>> GenJet_constituent_pt;
	std::vector<std::vector<float>> GenJet_constituent_eta;
	std::vector<std::vector<float>> GenJet_constituent_phi;

    // === Create trees ===
    std::vector<TTree*> trees;
    for (auto R : R_values) {

        TString rStr = Form("JetTree_R0p%d", int(R * 10)); // e.g. R=0.1 → R0p1
        TTree *tree = new TTree(rStr, rStr);

        tree->Branch("NEVENTS", &NEVENTS, "NEVENTS/I");

        // Reco
        tree->Branch("RecoJet_pt", &RecoJet_pt);
        tree->Branch("RecoJet_eta", &RecoJet_eta);
        tree->Branch("RecoJet_phi", &RecoJet_phi);
        tree->Branch("RecoJet_E", &RecoJet_E);
        tree->Branch("RecoJet_M", &RecoJet_M);
        tree->Branch("RecoJet_hasElectron", &RecoJet_hasElectron);
        tree->Branch("RecoJet_maxPtPart_pt", &RecoJet_maxPtPart_pt);
        tree->Branch("RecoJet_constituent_pt", &RecoJet_constituent_pt);
        tree->Branch("RecoJet_constituent_eta", &RecoJet_constituent_eta);
        tree->Branch("RecoJet_constituent_phi", &RecoJet_constituent_phi);
		tree->Branch("RecoJet_constituent_nhits", &RecoJet_constituent_nhits);
        // Gen
        tree->Branch("GenJet_pt", &GenJet_pt);
        tree->Branch("GenJet_eta", &GenJet_eta);
        tree->Branch("GenJet_phi", &GenJet_phi);
        tree->Branch("GenJet_E", &GenJet_E);
        tree->Branch("GenJet_M", &GenJet_M);
        tree->Branch("GenJet_hasElectron", &GenJet_hasElectron);
        tree->Branch("GenJet_hasNeutral", &GenJet_hasNeutral);
        tree->Branch("GenJet_maxPtPart_pt", &GenJet_maxPtPart_pt);
        tree->Branch("GenJet_constituent_pt", &GenJet_constituent_pt);
        tree->Branch("GenJet_constituent_eta", &GenJet_constituent_eta);
        tree->Branch("GenJet_constituent_phi", &GenJet_constituent_phi);

        trees.push_back(tree);
    }
	
	// Loop over events	
	int globalEvent = 0;
	while(tree_reader->Next()) {	
	
	    if(globalEvent%50000 == 0) cout << "Events Processed: " << globalEvent << endl;
	    globalEvent++;
	    NEVENTS = globalEvent;

        // Build pseudojets
        std::vector<PseudoJet> particles_reco;
        for (unsigned int i = 0; i < TrkRecoPx->GetSize(); ++i) {
            TVector3 mom((*TrkRecoPx)[i], (*TrkRecoPy)[i], (*TrkRecoPz)[i]);
            if (mom.Pt() < minCstPt || mom.Pt() > maxCstPt) continue;
            if (donhitcut) { if ( (*TrkRecoNhits)[i] < 4 ) continue; }
            if (removeelectrons){
				// Find electron
			    int chargePartIndex = i; 
			    int elecIndex = -1;
			    float elecIndexWeight = -1.0;
		    	for(unsigned int itrkass = 0; itrkass < TrkPartAssocRec->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
					if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
					    if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
							elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
							elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
			      		}
			  		}
		      	}
				if((*TrkMCGenPDG)[elecIndex] == 11) continue;
            }
            PseudoJet p((*TrkRecoPx)[i], (*TrkRecoPy)[i], (*TrkRecoPz)[i], (*TrkRecoE)[i]);
            p.set_user_index(i);
            particles_reco.push_back(p);
        }

        std::vector<PseudoJet> particles_gen;
        for (unsigned int i = 0; i < TrkGenPx->GetSize(); ++i) {
            TVector3 mom((*TrkGenPx)[i], (*TrkGenPy)[i], (*TrkGenPz)[i]);
            if (mom.Pt() < minCstPt || mom.Pt() > maxCstPt) continue;
            if ((*TrkGenCharge)[i] == 0) continue;
            if (removeelectrons) { if( (*TrkGenPDG)[i] == 11 ) continue; }
            PseudoJet p((*TrkGenPx)[i], (*TrkGenPy)[i], (*TrkGenPz)[i], (*TrkGenE)[i]);
            p.set_user_index(i);
            particles_gen.push_back(p);
        }

        // --- Loop over R values --- making jets
        for (size_t iR = 0; iR < R_values.size(); ++iR) { 
        
            double R = R_values[iR];
            TTree *tree = trees[iR];
            // Define algorithm
            JetAlgorithm algo = antikt_algorithm;
            RecombinationScheme scheme = E_scheme;
			// Jet definition
            JetDefinition jet_def(algo, R, scheme);
            GhostedAreaSpec ghost_spec(ghostMaxRap, numGhostRepeat, ghostArea);
            AreaType atype = active_area;
            AreaDefinition area_def(atype, ghost_spec);
                        
			// Clear all vectors for the new event
			RecoJet_pt.clear();
			RecoJet_eta.clear();
			RecoJet_phi.clear();
			RecoJet_E.clear();
			RecoJet_M.clear();
			RecoJet_hasElectron.clear();
			RecoJet_maxPtPart_pt.clear();
	        RecoJet_constituent_pt.clear(); 
	        RecoJet_constituent_eta.clear();
	        RecoJet_constituent_phi.clear(); 
	        RecoJet_constituent_nhits.clear();
		
			GenJet_pt.clear();
			GenJet_eta.clear();
			GenJet_phi.clear();
			GenJet_E.clear();
			GenJet_M.clear();
			GenJet_hasElectron.clear();
			GenJet_hasNeutral.clear();
			GenJet_maxPtPart_pt.clear();
			GenJet_constituent_pt.clear(); 
	        GenJet_constituent_eta.clear();
	        GenJet_constituent_phi.clear(); 

            // --- Reco clustering ---
            ClusterSequenceArea cs_reco(particles_reco, jet_def, area_def);
            std::vector<PseudoJet> jets_reco = sorted_by_pt(cs_reco.inclusive_jets(minJetPt));
            for (auto &jet : jets_reco) {
                RecoJet_pt.push_back(jet.pt());
                RecoJet_eta.push_back(jet.eta());
                RecoJet_phi.push_back(jet.phi_std());
                RecoJet_E.push_back(jet.e());
                RecoJet_M.push_back(jet.m());

                bool hasElectron = false;
                float maxPtReco = -1.0;
                std::vector<float> cpt, ceta, cphi;
                std::vector<int> chits;
                cpt.clear(); ceta.clear(); cphi.clear(); chits.clear();

                for (auto &c : jet.constituents()) {
                    int idx = c.user_index();
                    TVector3 v3((*TrkRecoPx)[idx], (*TrkRecoPy)[idx], (*TrkRecoPz)[idx]);
                    cpt.push_back(v3.Pt());
                    ceta.push_back(v3.Eta());
                    cphi.push_back(v3.Phi());
                    chits.push_back((*TrkRecoNhits)[idx]);
                    if (v3.Pt() > maxPtReco) maxPtReco = v3.Pt();
					// Find electron
			    	int chargePartIndex = idx; 
			    	int elecIndex = -1;
			    	float elecIndexWeight = -1.0;
		    		for(unsigned int itrkass = 0; itrkass < TrkPartAssocRec->GetSize(); itrkass++){ // Loop Over All ReconstructedChargedParticleAssociations
						if((*TrkPartAssocRec)[itrkass] == chargePartIndex){ // Select Entry Matching the ReconstructedChargedParticle Index
						    if((*TrkPartAssocWeight)[itrkass] > elecIndexWeight){ // Find Particle with Greatest Weight = Contributed Most Hits to Track
								elecIndex = (*TrkPartAssocSim)[itrkass]; // Get Index of MCParticle Associated with ReconstructedChargedParticle
								elecIndexWeight = (*TrkPartAssocWeight)[itrkass];
			     	 		}
			  			}
		      		}
					if((*TrkMCGenPDG)[elecIndex] == 11) hasElectron = true;
                }

                RecoJet_constituent_pt.push_back(cpt);
                RecoJet_constituent_eta.push_back(ceta);
                RecoJet_constituent_phi.push_back(cphi);
                RecoJet_constituent_nhits.push_back(chits);
                RecoJet_hasElectron.push_back(hasElectron);
                RecoJet_maxPtPart_pt.push_back(maxPtReco);
            }
            
        	// --- Gen clustering ---
            ClusterSequenceArea cs_gen(particles_gen, jet_def, area_def);
            std::vector<PseudoJet> jets_gen = sorted_by_pt(cs_gen.inclusive_jets(minJetPt));
            for (auto &jet : jets_gen) {
                GenJet_pt.push_back(jet.pt());
                GenJet_eta.push_back(jet.eta());
                GenJet_phi.push_back(jet.phi_std());
                GenJet_E.push_back(jet.e());
                GenJet_M.push_back(jet.m());

                bool hasGenElectron = false;
                bool hasGenNeutral = false;
                float maxPtGen = -1.0;
                std::vector<float> gpt, geta, gphi;
				gpt.clear(); geta.clear(); gphi.clear(); 

                for (auto &c : jet.constituents()) {
                    int idx = c.user_index();
                    TVector3 gv((*TrkGenPx)[idx], (*TrkGenPy)[idx], (*TrkGenPz)[idx]);
                    gpt.push_back(gv.Pt());
                    geta.push_back(gv.Eta());
                    gphi.push_back(gv.Phi());
                    if (gv.Pt() > maxPtGen) maxPtGen = gv.Pt();
                    if ((*TrkGenPDG)[idx] == 11) hasGenElectron = true;
                    if ((*TrkGenCharge)[idx] == 0) hasGenNeutral = true;
                }

                GenJet_constituent_pt.push_back(gpt);
                GenJet_constituent_eta.push_back(geta);
                GenJet_constituent_phi.push_back(gphi);
                GenJet_hasElectron.push_back(hasGenElectron);
                GenJet_hasNeutral.push_back(hasGenNeutral);
                GenJet_maxPtPart_pt.push_back(maxPtGen);
            }
			 
			tree->Fill();      
            
		}
		
	}
	
    // Write and close
    for (auto *tree : trees) tree->Write();
    OutFile->Close();

	cout << "Total number of events: " << globalEvent << endl;    
    

}

