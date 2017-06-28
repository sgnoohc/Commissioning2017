// TasUtil tool
#include "tasutil.cc"

// CORE tools
#include "CORE/CMS3.cc"
#include "CORE/Base.cc"
#include "CORE/TriggerSelections.cc"
#include "CORE/ElectronSelections.cc"
#include "CORE/MuonSelections.cc"
#include "CORE/IsolationTools.cc"
#include "CORE/Tools/goodrun.cc"
#include "CORE/Tools/JetCorrector.cc"

// typedef for simple naming
typedef std::map<TString, std::map<TString, TH1*>> Hist_DB;
void makeHists(TString prefix, Hist_DB& hist_DB);
void fillSingleLepHists(TString prefix, int lepid, int ilep, Hist_DB& hist_DB);
void fillNlepHists(TString prefix, std::vector<int>& good_muon_idx, std::vector<int>& good_elec_idx, Hist_DB& hist_DB);
void fillDilepHists(TString prefix, int lepid0, int lepid1, int ilep0, int ilep1, Hist_DB& hist_DB);
void saveHists(Hist_DB& hist_DB, TString output_prefix);
void fill(TString prefix, std::vector<int>& good_muon_idx, std::vector<int>& good_elec_idx, Hist_DB& hist_DB);
float Mll(std::vector<int>& good_muon_idx, std::vector<int>& good_elec_idx);
bool passHLTs(std::vector<TString>);
void ScanChain(TChain* chain, TString output_prefix, int nevents=-1);

//_________________________________________________________________________________________________
void ScanChain(TChain* chain, TString output_name, int nevents)
{

    // -~-~-~-~-~-~-~-~-~-~-
    // Set triggers to study
    // -~-~-~-~-~-~-~-~-~-~-

    // List of triggers intrested
    std::vector<TString> trigger_names;
    trigger_names.push_back("HLT_IsoTkMu27_v");
    trigger_names.push_back("HLT_IsoTkMu24_v");
    trigger_names.push_back("HLT_IsoMu27_v");
    trigger_names.push_back("HLT_IsoMu24_v");
    trigger_names.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    trigger_names.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    trigger_names.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    trigger_names.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    trigger_names.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    trigger_names.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");

    std::map<TString, std::vector<TString>> combined_trigger_names;
    combined_trigger_names["SingleEl"] = {"HLT_Ele27_WPTight_Gsf_v", "HLT_Ele35_WPTight_Gsf_v", "HLT_Ele38_WPTight_Gsf_v", "HLT_Ele40_WPTight_Gsf_v"};
    combined_trigger_names["SingleMu"] = {"HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50"};
    combined_trigger_names["DoubleEl"] = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
    combined_trigger_names["DoubleMu"] = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"};
    combined_trigger_names["DilepEMu"] = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};

    // -~-~-~-~-~-~-~-~-~
    // Set Good Runs List
    // -~-~-~-~-~-~-~-~-~

    const char* json_file = "lumi/jsons/json_DCSONLY_170622_snt.txt";
    set_goodrun_file(json_file);

    // -~-~-~-~-~-~
    // Create Hists
    // -~-~-~-~-~-~

    // All our TH1s saved in a map
    Hist_DB hist_DB;

    // Create list of histograms for a given triggers
    for (auto& trigger_name : trigger_names)
    {
        makeHists(trigger_name, hist_DB);
        makeHists(trigger_name+"_onZ", hist_DB);
    }

    // Create list of histograms for a given triggers
    for (auto& map_combined_trigger_name : combined_trigger_names)
    {
        TString combined_trigger_name = map_combined_trigger_name.first;
        makeHists(combined_trigger_name, hist_DB);
        makeHists(combined_trigger_name+"_onZ", hist_DB);
    }

    // -~-~-~-~-~
    // Event Loop
    // -~-~-~-~-~

    TasUtil::Looper<CMS3> looper(chain, &cms3, nevents);
    while (looper.nextEvent())
    {

        // -~-~-~-~-
        // Check GRL
        // -~-~-~-~-

        // If not a good run event skip
        if (!goodrun(cms3.evt_run(), cms3.evt_lumiBlock()))
            continue;

        // -~-~-~-~-~-~-~-~-~-
        // Select Good Leptons
        // -~-~-~-~-~-~-~-~-~-

        // Good muon idx
        std::vector<int> good_muon_idx;

        // Good elec idx
        std::vector<int> good_elec_idx;

        // Loop over Muons
        for (unsigned int imu = 0; imu < cms3.mus_p4().size(); ++imu)
        {
            // If ID fail pass
            if (!muonID(imu, id_level_t::HAD_loose_v4))
                continue;
            good_muon_idx.push_back(imu);
        }

        // Loop over electrons
        for (unsigned int iel = 0; iel < cms3.els_p4().size(); ++iel)
        {
            // If ID fail pass
            if (!electronID(iel, id_level_t::HAD_loose_v4))
                continue;
            good_elec_idx.push_back(iel);
        }

        // -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
        // Loop over individual triggers to study
        // -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

        // Loop over triggers individually
        for (auto& trigger_name : trigger_names)
        {
            // If the trigger doesn't fire then skip
            if (passHLTTriggerPattern(trigger_name))
                continue;

            // -~-~-~-~-~-~-~-~-~-
            // Fill the histograms
            // -~-~-~-~-~-~-~-~-~-

            fill(trigger_name, good_muon_idx, good_elec_idx, hist_DB);
        }

        // -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
        // Loop over combined triggers to study
        // -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

        for (auto& map_combined_trigger_name : combined_trigger_names)
        {
            TString combined_trigger_name = map_combined_trigger_name.first;
            std::vector<TString> trigger_names = map_combined_trigger_name.second;
            if (!passHLTs(trigger_names))
                continue;

            // -~-~-~-~-~-~-~-~-~-
            // Fill the histograms
            // -~-~-~-~-~-~-~-~-~-

            fill(combined_trigger_name, good_muon_idx, good_elec_idx, hist_DB);
        }

    }

    // Save the histograms to output.root
    saveHists(hist_DB, output_name);

}

//_________________________________________________________________________________________________
void makeHists(TString prefix, Hist_DB& hist_DB)
{
    // Make histograms
    // Lepton multiplicity related histograms
    hist_DB[prefix]["h_nlep"]                   = new TH1D(prefix + "_" + "h_nlep"                   , "" , 50    , 0.    , 50.    );
    hist_DB[prefix]["h_1lep_flav"]              = new TH1D(prefix + "_" + "h_1lep_flav"              , "" , 30    , -15.  , 15.    );
    hist_DB[prefix]["h_2lep_flav"]              = new TH1D(prefix + "_" + "h_2lep_flav"              , "" , 6     , 0.    , 6.     );
    // Single lepton related histograms
    hist_DB[prefix]["h_lep_pt"]                 = new TH1D(prefix + "_" + "h_lep_pt"                 , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_el_pt"]                  = new TH1D(prefix + "_" + "h_el_pt"                  , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_mu_pt"]                  = new TH1D(prefix + "_" + "h_mu_pt"                  , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_lep_eta"]                = new TH1D(prefix + "_" + "h_lep_eta"                , "" , 10000 , -2.5  , 2.5    );
    hist_DB[prefix]["h_el_eta"]                 = new TH1D(prefix + "_" + "h_el_eta"                 , "" , 10000 , -2.5  , 2.5    );
    hist_DB[prefix]["h_mu_eta"]                 = new TH1D(prefix + "_" + "h_mu_eta"                 , "" , 10000 , -2.5  , 2.5    );
    hist_DB[prefix]["h_lep_phi"]                = new TH1D(prefix + "_" + "h_lep_phi"                , "" , 10000 , -3.15 , 3.15   );
    hist_DB[prefix]["h_el_phi"]                 = new TH1D(prefix + "_" + "h_el_phi"                 , "" , 10000 , -3.15 , 3.15   );
    hist_DB[prefix]["h_mu_phi"]                 = new TH1D(prefix + "_" + "h_mu_phi"                 , "" , 10000 , -3.15 , 3.15   );
    // Dilepton related histograms
    hist_DB[prefix]["h_mll_ll"]                 = new TH1D(prefix + "_" + "h_mll_ll"                 , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_mll_ee"]                 = new TH1D(prefix + "_" + "h_mll_ee"                 , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_mll_em"]                 = new TH1D(prefix + "_" + "h_mll_em"                 , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_mll_mm"]                 = new TH1D(prefix + "_" + "h_mll_mm"                 , "" , 10000 , 0.    , 10000. );
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_ll"] = new TH1D(prefix + "_" + "h_mll_zoomed_on_Zpeak_ll" , "" , 10000 , 50.   , 150.   );
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_ee"] = new TH1D(prefix + "_" + "h_mll_zoomed_on_Zpeak_ee" , "" , 10000 , 50.   , 150.   );
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_em"] = new TH1D(prefix + "_" + "h_mll_zoomed_on_Zpeak_em" , "" , 10000 , 50.   , 150.   );
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_mm"] = new TH1D(prefix + "_" + "h_mll_zoomed_on_Zpeak_mm" , "" , 10000 , 50.   , 150.   );
    // Lepton multiplicity related histograms
    hist_DB[prefix]["h_nlep"]                   -> Sumw2();
    hist_DB[prefix]["h_1lep_flav"]              -> Sumw2();
    hist_DB[prefix]["h_2lep_flav"]              -> Sumw2();
    // Single lepton related histograms
    hist_DB[prefix]["h_lep_pt"]                 -> Sumw2();
    hist_DB[prefix]["h_el_pt"]                  -> Sumw2();
    hist_DB[prefix]["h_mu_pt"]                  -> Sumw2();
    hist_DB[prefix]["h_lep_eta"]                -> Sumw2();
    hist_DB[prefix]["h_el_eta"]                 -> Sumw2();
    hist_DB[prefix]["h_mu_eta"]                 -> Sumw2();
    hist_DB[prefix]["h_lep_phi"]                -> Sumw2();
    hist_DB[prefix]["h_el_phi"]                 -> Sumw2();
    hist_DB[prefix]["h_mu_phi"]                 -> Sumw2();
    // Dilepton related histograms
    hist_DB[prefix]["h_mll_ll"]                 -> Sumw2();
    hist_DB[prefix]["h_mll_ee"]                 -> Sumw2();
    hist_DB[prefix]["h_mll_em"]                 -> Sumw2();
    hist_DB[prefix]["h_mll_mm"]                 -> Sumw2();
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_ll"] -> Sumw2();
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_ee"] -> Sumw2();
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_em"] -> Sumw2();
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_mm"] -> Sumw2();
}

//_________________________________________________________________________________________________
void fillSingleLepHists(
        TString prefix,
        int lepid,
        int ilep,
        Hist_DB& hist_DB)
{
    double wgt = 1.;
    if (!cms3.evt_isRealData())
        wgt = cms3.evt_scale1fb();
    const LorentzVector& p4 = abs(lepid) == 11 ? cms3.els_p4()[ilep] : cms3.mus_p4()[ilep];
    hist_DB[prefix]["h_lep_pt"]  ->Fill(p4.pt(), wgt);
    hist_DB[prefix]["h_lep_eta"] ->Fill(p4.eta(), wgt);
    hist_DB[prefix]["h_lep_phi"] ->Fill(p4.phi(), wgt);
    switch (abs(lepid))
    {
        case 11:
            hist_DB[prefix]["h_el_pt"]   ->Fill(p4.pt(), wgt);
            hist_DB[prefix]["h_el_eta"]  ->Fill(p4.eta(), wgt);
            hist_DB[prefix]["h_el_phi"]  ->Fill(p4.phi(), wgt);
            break;
        case 13:
            hist_DB[prefix]["h_mu_pt"]   ->Fill(p4.pt(), wgt);
            hist_DB[prefix]["h_mu_eta"]  ->Fill(p4.eta(), wgt);
            hist_DB[prefix]["h_mu_phi"]  ->Fill(p4.phi(), wgt);
            break;
    }
}

//_________________________________________________________________________________________________
void fillNlepHists(
        TString prefix,
        std::vector<int>& good_muon_idx,
        std::vector<int>& good_elec_idx,
        Hist_DB& hist_DB)
{
    double wgt = 1.;
    if (!cms3.evt_isRealData())
        wgt = cms3.evt_scale1fb();

    // Fill Nlep
    hist_DB[prefix]["h_nlep"]->Fill(good_muon_idx.size() + good_elec_idx.size(), wgt);

    if (good_muon_idx.size() + good_elec_idx.size() == 1)
    {
        if (good_muon_idx.size() == 1)
        {
            hist_DB[prefix]["h_1lep_flav"]->Fill(cms3.mus_charge()[good_muon_idx[0]] * 13, wgt);
        }
        else if (good_elec_idx.size() == 1)
        {
            hist_DB[prefix]["h_1lep_flav"]->Fill(cms3.els_charge()[good_elec_idx[0]] * 11, wgt);
        }
    }
    else if (good_muon_idx.size() + good_elec_idx.size() == 2)
    {
        // 0 = SSee
        // 1 = OSee
        // 2 = SSem
        // 3 = OSem
        // 4 = SSmm
        // 5 = OSmm
        if (good_muon_idx.size() == 2)
        {
            if (cms3.mus_charge()[good_muon_idx[0]] != -cms3.mus_charge()[good_muon_idx[1]])
                hist_DB[prefix]["h_2lep_flav"]->Fill(4., wgt);
            else
                hist_DB[prefix]["h_2lep_flav"]->Fill(5., wgt);
        }
        else if (good_elec_idx.size() == 2)
        {
            if (cms3.els_charge()[good_elec_idx[0]] != -cms3.els_charge()[good_elec_idx[1]])
                hist_DB[prefix]["h_2lep_flav"]->Fill(0., wgt);
            else
                hist_DB[prefix]["h_2lep_flav"]->Fill(1., wgt);
        }
        else
        {
            if (cms3.mus_charge()[good_muon_idx[0]] != -cms3.els_charge()[good_elec_idx[0]])
                hist_DB[prefix]["h_2lep_flav"]->Fill(2., wgt);
            else
                hist_DB[prefix]["h_2lep_flav"]->Fill(3., wgt);
        }
    }

}

//_______________________________________________________________________________________________
void fillDilepHists(
        TString prefix,
        int lepid0,
        int lepid1,
        int ilep0,
        int ilep1,
        Hist_DB& hist_DB)
{
    double wgt = 1.;
    if (!cms3.evt_isRealData())
        wgt = cms3.evt_scale1fb();
    const LorentzVector& p4_0 = abs(lepid0) == 11 ? cms3.els_p4()[ilep0] : cms3.mus_p4()[ilep0];
    const LorentzVector& p4_1 = abs(lepid1) == 11 ? cms3.els_p4()[ilep1] : cms3.mus_p4()[ilep1];
    float mll = (p4_0 + p4_1).M();
    hist_DB[prefix]["h_mll_ll"]      ->Fill(mll, wgt);
    hist_DB[prefix]["h_mll_zoomed_on_Zpeak_ll"]  ->Fill(mll, wgt);
    switch (lepid0 * lepid1)
    {
        case -121:
            hist_DB[prefix]["h_mll_ee"]      ->Fill(mll, wgt);
            hist_DB[prefix]["h_mll_zoomed_on_Zpeak_ee"]  ->Fill(mll, wgt);
            break;
        case -143:
            hist_DB[prefix]["h_mll_em"]      ->Fill(mll, wgt);
            hist_DB[prefix]["h_mll_zoomed_on_Zpeak_em"]  ->Fill(mll, wgt);
            break;
        case -169:
            hist_DB[prefix]["h_mll_mm"]      ->Fill(mll, wgt);
            hist_DB[prefix]["h_mll_zoomed_on_Zpeak_mm"]  ->Fill(mll, wgt);
            break;
    }
}

//_________________________________________________________________________________________________
void saveHists(Hist_DB& hist_DB, TString output_name)
{
    TFile* ofile = new TFile(output_name, "recreate");
    TasUtil::print("Writing outputs to: " + output_name);

    for (auto& map_map_hist : hist_DB)
    {
        TString trigger_name = map_map_hist.first;
        for (auto& map_hist : map_map_hist.second)
        {
            TString hist_name = map_hist.first;
            TH1* hist = map_hist.second;
            ofile->cd();
            hist->Write();
        }
    }

    ofile->Close();
}

//_________________________________________________________________________________________________
float Mll(std::vector<int>& good_muon_idx, std::vector<int>& good_elec_idx)
{
    if (good_muon_idx.size() + good_elec_idx.size() != 2)
        return -999;
    float mll = -999;
    if (good_muon_idx.size() == 2)
    {
        mll = (cms3.mus_p4()[good_muon_idx[0]] + cms3.mus_p4()[good_muon_idx[1]]).M();
        if (cms3.mus_charge()[good_muon_idx[0]] != -cms3.mus_charge()[good_muon_idx[1]])
            return -999;
    }
    else if (good_elec_idx.size() == 2)
    {
        mll = (cms3.els_p4()[good_elec_idx[0]] + cms3.els_p4()[good_elec_idx[1]]).M();
        if (cms3.els_charge()[good_elec_idx[0]] != -cms3.els_charge()[good_elec_idx[1]])
            return -999;
    }
    else
    {
        mll = (cms3.mus_p4()[good_muon_idx[0]] + cms3.els_p4()[good_elec_idx[0]]).M();
        if (cms3.mus_charge()[good_muon_idx[0]] != -cms3.els_charge()[good_elec_idx[0]])
            return -999;
    }
    return mll;
}

//_________________________________________________________________________________________________
bool passHLTs(std::vector<TString> trigger_names)
{
    for (auto& trigger_name : trigger_names)
    {
        if (passHLTTriggerPattern(trigger_name))
            return true;
    }
    return false;
}

//_________________________________________________________________________________________________
void fill(TString trigger_name, std::vector<int>& good_muon_idx, std::vector<int>& good_elec_idx, Hist_DB& hist_DB)
{

    // Fill lepton multiplicity related histograms
    fillNlepHists(trigger_name, good_muon_idx, good_elec_idx, hist_DB);

    // Fill all leptons
    for (auto& idx : good_muon_idx) fillSingleLepHists(trigger_name, 13, idx, hist_DB);
    for (auto& idx : good_elec_idx) fillSingleLepHists(trigger_name, 11, idx, hist_DB);

    // Select, opposite-sign and on Z. (Mll returns -999 for same sign pair)
    if ((good_elec_idx.size() + good_muon_idx.size() == 2) && abs(Mll(good_muon_idx, good_elec_idx) - 90.) < 10.)
    {
        TString prefix = trigger_name + "_onZ";
        for (auto& idx : good_muon_idx) fillSingleLepHists(prefix, 13, idx, hist_DB);
        for (auto& idx : good_elec_idx) fillSingleLepHists(prefix, 11, idx, hist_DB);
        int lepid0 = -999;
        int lepid1 = -999;
        int ilep0 = -999;
        int ilep1 = -999;
        if (good_muon_idx.size() == 2)
        {
            lepid0 =  13;
            lepid1 = -13;
            ilep0 = good_muon_idx[0];
            ilep1 = good_muon_idx[1];
        }
        else if (good_elec_idx.size() == 2)
        {
            lepid0 =  11;
            lepid1 = -11;
            ilep0 = good_elec_idx[0];
            ilep1 = good_elec_idx[1];
        }
        else
        {
            lepid0 =  11;
            lepid1 = -13;
            ilep0 = good_elec_idx[0];
            ilep1 = good_muon_idx[0];
        }
        fillDilepHists(prefix, lepid0, lepid1, ilep0, ilep1, hist_DB);

    }
}
