
# Commissioning work

## The datasets

The new 2017 data sits here:
    /hadoop/cms/store/user/namin/ProjectMetis/

    #(ls'ed the directory June 26 10:54 am)

    $ ls /hadoop/cms/store/user/namin/ProjectMetis/
    DoubleEG_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03    MET_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03
    DoubleEG_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03    MET_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    DoubleEG_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03    MuonEG_Run2017A-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    DoubleMuon_Run2017A-PromptReco-v1_MINIAOD_CMS4_V00-00-03  MuonEG_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03
    DoubleMuon_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03  MuonEG_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03
    DoubleMuon_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03  MuonEG_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    DoubleMuon_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03  SingleElectron_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03
    HTMHT_Run2017A-PromptReco-v1_MINIAOD_CMS4_V00-00-03       SingleElectron_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03
    HTMHT_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03       SingleElectron_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    HTMHT_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03       SingleMuon_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03
    HTMHT_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03       SingleMuon_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03
    JetHT_Run2017A-PromptReco-v1_MINIAOD_CMS4_V00-00-03       SingleMuon_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    JetHT_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03       SinglePhoton_Run2017A-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    JetHT_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03       SinglePhoton_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03
    JetHT_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03       SinglePhoton_Run2017A-PromptReco-v3_MINIAOD_CMS4_V00-00-03
    MET_Run2017A-PromptReco-v1_MINIAOD_CMS4_V00-00-03         SinglePhoton_Run2017B-PromptReco-v1_MINIAOD_CMS4_V00-00-03
    MET_Run2017A-PromptReco-v2_MINIAOD_CMS4_V00-00-03         TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8

Some test ntuple results from Sicheng
    /hadoop/cms/store/group/snt/run2_data_test/DoubleMuon_Run2017A-PromptReco-v2/V00-00-03_workaround/merged_ntuple_*.root

## Setting up a simple CMS4 looper that runs over the 2017 dataset ntuples.

First, setup some tools to the current working directory.

    cp ~phchang/public_html/tasutil/tasutil.h .
    cp ~phchang/public_html/tasutil/tasutil.cc .

Check out cmstas/CORE

    git clone git@github.com:cmstas/CORE
    cd CORE
    git checkout cms4
    cd ../

Then in file ScanChan.C copy paste

    #include "CORE/CMS3.cc"
    #include "CORE/Base.cc"
    #include "tasutil.cc"

    void ScanChain(TChain* chain)
    {
        TasUtil::Looper<CMS3> looper(chain, &cms3, 1000);
        while (looper.nextEvent())
        {
            // Do stuff here for your commissioning work
        }
    }

In file doAll.C copy paste

    {
        // Compile your ScanChain.C with "+".
        gROOT->LoadMacro("ScanChain.C+");

        // Setup the TChain that gets fed to ScanChain.C
        TChain *ch = new TChain("Events");

        // Then add a list of samples you want to loop over
        ch->Add("/hadoop/cms/store/group/snt/run2_data_test/DoubleMuon_Run2017A-PromptReco-v2/V00-00-03_workaround/merged_ntuple_*.root");
        ch->Add("/hadoop/cms/store/group/snt/run2_data_test/DoubleEG_Run2017A-PromptReco-v2/V00-00-03_workaround/merged_ntuple_*.root");

        // LOOP!
        ScanChain(ch);
    }

Happy coding!
