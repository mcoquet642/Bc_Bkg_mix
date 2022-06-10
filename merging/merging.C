void merging(){
//    std::vector<string> filelist = {"../trimu_bkg_Jpsi_bkg_mu_M0-100_PT_0.root", "../trimu_prompt_Jpsi_bkg_mu_M0-100_PT_0.root", "../trimu_nonprompt_Jpsi_bkg_mu_M0-100_PT_0.root", "../sig_sig/trimuon_signal_signal.root"};
    //std::vector<string> filelist = {"../test_bkg.root", "../test_prompt.root", "../test_nonprompt.root", "../sig_sig/trimuon_signal_signal.root"};
//    std::vector<string> filelist = {"../test_bkg.root", "../test_prompt_01.root", "../test_nonprompt_01.root", "../sig_sig/trimuon_signal_signal_01.root"};
    std::vector<string> filelist = {"../test_bkg.root", "../test_prompt.root", "../test_nonprompt_02.root", "../sig_sig/trimuon_signal_signal_02.root"};
    std::vector<TTree*> ttreelist;
    std::vector<TFile*> tfilelist;
    std::vector<TFile*> outtfilelist;
    std::vector<TTree*> clonedTrees;
    for (int i = 0; i < (int)filelist.size(); ++i) {
        tfilelist.push_back(TFile::Open(filelist[i].c_str()));
        ttreelist.push_back(dynamic_cast<TTree*>(tfilelist.back()->Get("BcDecay")));
        outtfilelist.push_back(TFile::Open(("out" + std::to_string(i) + ".root").c_str(), "RECREATE"));
        clonedTrees.push_back(dynamic_cast<TTree*>(ttreelist[0]->CloneTree(0)));
//	ttreelist[i]->Print();
    }
    cout << (int)ttreelist.size() << endl;
    for (int i = 0; i < (int)ttreelist.size(); ++i) {
        for (auto&& t : clonedTrees) ttreelist[i]->CopyAddresses(t);
        std::mt19937 mt(std::random_device{}());
        std::uniform_int_distribution<size_t> dist(0, ttreelist.size()-1);
        auto localEntries = ttreelist[i]->GetEntries();
	cout << localEntries << endl;
        for (Long64_t j = 0; j < localEntries; ++j) {
            ttreelist[i]->GetEntry(j);
            clonedTrees[dist(mt)]->Fill();
        }
        for (auto&& t : clonedTrees) ttreelist[i]->CopyAddresses(t, true);
        tfilelist[i]->Close();
    }

    TTree* shTree = clonedTrees[0]->CloneTree(0);
    TList *list = new TList;
    for (int i = 0; i < (int)clonedTrees.size(); ++i) {
        auto localEntries = clonedTrees[i]->GetEntries();

        std::vector<Long64_t> ev(localEntries);
        std::iota(ev.begin(), ev.end(), 0);
        std::shuffle(ev.begin(), ev.end(), std::mt19937{std::random_device{}()});

        clonedTrees[i]->CopyAddresses(shTree);
        clonedTrees[i]->LoadBaskets();
        for (Long64_t j = 0; j < localEntries; ++j) {
            clonedTrees[i]->GetEntry(ev[j]);
            shTree->Fill();
        }
        clonedTrees[i]->CopyAddresses(shTree, true);
        list->Add(clonedTrees[i]);
//	clonedTrees[i]->Print();
    }
//    TFile f("out.root", "recreate");
//    TFile f("out_01.root", "recreate");
    TFile f("out_02.root", "recreate");
    TTree *newtree = TTree::MergeTrees(list);
    newtree->SetName("BcDecay");
    f.cd();
    newtree->Write();
    f.Close();
    for (int i = 0; i < (int)clonedTrees.size(); ++i) {
    	clonedTrees[i]->Delete();
    }
}
