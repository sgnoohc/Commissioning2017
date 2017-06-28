//  .
// ..: P. Chang, philip@physics.ucsd.edu

// N.B.: I kept a very compact style in listing the functions in order to easily move around with {,} in vim.
//       Therefore, if one wants to add new feature please make the function extremely short as possible.
//       If the function is too long, maybe it doesn't belong here!

#ifndef tasutil_h
#define tasutil_h

// C
#include <map>
#include <vector>
#include <stdarg.h>
#include <functional>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TChainElement.h"
#include "TTreeCache.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"

namespace TasUtil
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Printing functions
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // No namespace given in order to minimize typing
    // (e.g. TasUtil::print v. TasUtil::NAMESPACE::print)
    void print   (TString msg="", const char* fname="", int flush_before=0, int flush_after=0);
    void error   (TString msg, const char* fname, int is_error=1);
    void warning (TString msg, const char* fname);
    void announce(TString msg="", int quiet=0);
    void start   (int quiet=0, int sleep_time=0);
    void exit    (int quiet=0);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Particle class
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // This class holds every information about a "particle".
    // The class has a generic data structure using std::map<TString, 'data'>.
    // User can use it for whatever reason.
    class Particle
    {

        enum VarType
        {
            kINT = 0,
            kFLT = 1,
            kTLV = 2,
            kSTR = 3,
            kNON = -1
        };

        typedef std::vector<TString> vecStr;
        typedef std::map<TString, VarType> mapTyp;
        typedef std::map<TString, float> mapFlt;
        typedef std::map<TString, int> mapInt;
        typedef std::map<TString, TLorentzVector> mapTLV;
        typedef std::map<TString, TString> mapStr;

        private:

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Members
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // Vector to hold the types
        vecStr var_Nam;
        mapStr var_Ttl;
        mapTyp var_Typ;
        mapFlt var_Flt;
        mapInt var_Int;
        mapTLV var_TLV;
        mapStr var_Str;

        public:

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Functions
        //////////////////////////////////////////////////////////////////////////////////////////////////

        Particle();
        ~Particle();

        bool varExists(TString name);
        void createFltVar(TString name, TString title="");
        void createIntVar(TString name, TString title="");
        void createTLVVar(TString name, TString title="");
        void createStrVar(TString name, TString title="");

        static void printError             (TString name, TString msg, TString action, VarType type=kNON);
        static void printAlreadyExistsError(TString name, TString action, VarType type=kNON);
        static void printOutOfRangeError   (TString name, TString action, VarType type=kNON);

        void setFltVar(TString name, float          var);
        void setIntVar(TString name, int            var);
        void setStrVar(TString name, TString        var);
        void setTLVVar(TString name, TLorentzVector var);

        const float&          getFltVar(TString name) const;
        const int&            getIntVar(TString name) const;
        const TLorentzVector& getTLVVar(TString name) const;
        const TString&        getStrVar(TString name) const;

        VarType getType(TString name);

        void print();
        void printFltVar(TString name);
        void printIntVar(TString name);
        void printTLVVar(TString name);
        void printStrVar(TString name);

    };
    typedef std::vector<Particle> Particles;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Data extractor
    ///////////////////////////////////////////////////////////////////////////////////////////////
    Particles get(std::function<int()>, std::function<bool(int)>, std::function<Particle(int)>);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Math functions
    ///////////////////////////////////////////////////////////////////////////////////////////////
    namespace Math
    {
        // Short hand for float LorentzVector
        typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
        ///////////////////////////////////////////////////////////////////////////////////////////
        // TLorentzVector related operations
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Create TLVs
        TLorentzVector TLV(LorentzVector);
        TLorentzVector TLVXYZE(float, float, float, float);
        TLorentzVector TLVPtEtaPhiM(float, float, float, float);
        TLorentzVector TLVPtEtaPhiE(float, float, float, float);
        TLorentzVector TLVMET(float, float);
        // The custom MT function, not from the TLorentzVector member function
        float MT  (TLorentzVector, TLorentzVector);
        // TLorentzVector
        float DEta(TLorentzVector a, TLorentzVector b);
        float DPhi(TLorentzVector a, TLorentzVector b);
        float DR  (TLorentzVector a, TLorentzVector b);
        float DPt (TLorentzVector a, TLorentzVector b);
        float Mass(TLorentzVector a, TLorentzVector b);
        float Pt  (TLorentzVector a, TLorentzVector b);
        float MT  (TLorentzVector p4  , float met, float met_phi);
        ///////////////////////////////////////////////////////////////////////////////////////////
        // LorentzVector related operations
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Creating LorentzVector
        LorentzVector LV(TLorentzVector);
        LorentzVector LVXYZE(float, float, float, float);
        LorentzVector METLV(float, float);
        // LorentzVector
        float DEta(LorentzVector a, float         b);
        float DPhi(LorentzVector a, float         b);
        float DEta(LorentzVector a, LorentzVector b);
        float DPhi(LorentzVector a, LorentzVector b);
        float DR  (LorentzVector a, LorentzVector b);
        float DPt (LorentzVector a, LorentzVector b);
        float Mass(LorentzVector a, LorentzVector b);
        float Pt  (LorentzVector a, LorentzVector b);
        float MT  (LorentzVector a, LorentzVector b);
        float MT  (LorentzVector a, float met, float met_phi);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Looper class
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // NOTE: This class assumes accessing TTree in the SNT style which uses the following,
    // https://github.com/cmstas/Software/blob/master/makeCMS3ClassFiles/makeCMS3ClassFiles.C
    // It is assumed that the "template" class passed to this class will have
    // 1. "Init(TTree*)"
    // 2. "GetEntry(uint)"
    // 3. "progress(nevtProc'ed, total)"
    template <class TREECLASS>
    class Looper
    {
        // Members
        TChain* tchain;
        TObjArray *listOfFiles;
        TObjArrayIter* fileIter;
        TFile* tfile;
        TTree* ttree;
        unsigned int nEventsTotalInChain;
        unsigned int nEventsTotalInTree;
        int nEventsToProcess;
        unsigned int nEventsProcessed;
        unsigned int indexOfEventInTTree;
        bool fastmode;
        TREECLASS* treeclass;
        public:
        // Functions
        Looper(TChain* chain=0, TREECLASS* treeclass=0, int nEventsToProcess=-1);
        ~Looper();
        void setTChain(TChain* c);
        void setTreeClass(TREECLASS* t);
        bool allEventsInTreeProcessed();
        bool allEventsInChainProcessed();
        bool nextEvent();
        TTree* getTree() { return ttree; }
        unsigned int getNEventsProcessed() { return nEventsProcessed; }
        private:
        void setFileList();
        void setNEventsToProcess();
        bool nextTree();
        bool nextEventInTree();
    };

}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Event Looper (Looper) class template implementation
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// It's easier to put the implementation in the header file to avoid compilation issues.

//_________________________________________________________________________________________________
template <class TREECLASS>
TasUtil::Looper<TREECLASS>::Looper(TChain* c, TREECLASS* t, int nevtToProc) :
    tchain(0),
    listOfFiles(0),
    fileIter(0),
    tfile(0),
    ttree(0),
    nEventsTotalInChain(0),
    nEventsTotalInTree(0),
    nEventsToProcess(nevtToProc),
    nEventsProcessed(0),
    indexOfEventInTTree(0),
    fastmode(true),
    treeclass(0)
{
    print("Start EventLooping");
    start();
    if (c) setTChain(c);
    if (t) setTreeClass(t);
}

//_________________________________________________________________________________________________
template <class TREECLASS>
TasUtil::Looper<TREECLASS>::~Looper()
{
    print("Finished EventLooping");
    exit();
    if (fileIter) delete fileIter;
    if (tfile) delete tfile;
}

//_________________________________________________________________________________________________
template <class TREECLASS>
void TasUtil::Looper<TREECLASS>::setTChain(TChain* c)
{
    if (c)
    {
        tchain = c;
        setNEventsToProcess();
        setFileList();
    }
    else
    {
        error("You provided a null TChain pointer!", __FUNCTION__);
    }
}

//_________________________________________________________________________________________________
template <class TREECLASS>
void TasUtil::Looper<TREECLASS>::setTreeClass(TREECLASS* t)
{
    if (t)
    {
        treeclass = t;
    }
    else
    {
        error("You provided a null TreeClass pointer!", __FUNCTION__);
    }
}

//_________________________________________________________________________________________________
template <class TREECLASS>
bool TasUtil::Looper<TREECLASS>::nextTree()
{
    if (!fileIter)
        error("fileIter not set but you are trying to access the next file", __FUNCTION__);
    // Get the TChainElement from TObjArrayIter.
    // If no more to run over, Next returns 0.
    TChainElement* chainelement = (TChainElement*) fileIter->Next();
    if (chainelement)
    {
        // If there is already a TFile opened from previous iteration, close it.
        if (tfile) tfile->Close();
        // Open up a new file
        tfile = new TFile(chainelement->GetTitle());
        // Get the ttree
        ttree = (TTree*) tfile->Get(tchain->GetName());
        if (!ttree)
            error("TTree is null!??", __FUNCTION__);
        // Set some fast mode stuff
        if (fastmode) TTreeCache::SetLearnEntries(10);
        if (fastmode) ttree->SetCacheSize(128*1024*1024);
        // Print some info to stdout
        announce("Working on " +
                 TString(tfile->GetName()) +
                 "/TTree:" +
                 TString(ttree->GetName()));
        // Reset the nEventsTotalInTree in this tree
        nEventsTotalInTree = ttree->GetEntries();
        // Reset the event index as we got a new ttree
        indexOfEventInTTree = 0;
        // Set the ttree to the TREECLASS
        treeclass->Init(ttree);
        // Return that I got a good one
        return true;
    }
    else
    {
        // Announce that we are done with this chain
        print("Done with all trees in this chain", __FUNCTION__);
        return false;
    }
}

//_________________________________________________________________________________________________
template <class TREECLASS>
bool TasUtil::Looper<TREECLASS>::allEventsInTreeProcessed()
{
    if (indexOfEventInTTree >= nEventsTotalInTree) return true;
    else return false;
}

//_________________________________________________________________________________________________
template <class TREECLASS>
bool TasUtil::Looper<TREECLASS>::allEventsInChainProcessed()
{
    if (nEventsProcessed >= (unsigned int) nEventsToProcess) return true;
    else return false;
}

//_________________________________________________________________________________________________
template <class TREECLASS>
bool TasUtil::Looper<TREECLASS>::nextEventInTree()
{
    // Sanity check before loading the next event.
    if (!ttree) error("current ttree not set!", __FUNCTION__);
    if (!tfile) error("current tfile not set!", __FUNCTION__);
    if (!fileIter) error("fileIter not set!", __FUNCTION__);
    // Check whether I processed everything
    if (allEventsInTreeProcessed()) return false;
    if (allEventsInChainProcessed()) return false;
    // if fast mode do some extra
    if (fastmode) ttree->LoadTree(indexOfEventInTTree);
    // Set the event index in TREECLASS
    treeclass->GetEntry(indexOfEventInTTree);
    // Print progress
    treeclass->progress(nEventsProcessed, nEventsToProcess);
    // Increment the counter for this ttree
    ++indexOfEventInTTree;
    // Increment the counter for the entire tchain
    ++nEventsProcessed;
    // If all fine return true
    return true;
}

//_________________________________________________________________________________________________
template <class TREECLASS>
bool TasUtil::Looper<TREECLASS>::nextEvent()
{
    // If no tree it means this is the beginning of the loop.
    if (!ttree)
    {
        // Load the next tree if it returns true, then proceed to next event in tree.
        while (nextTree())
        {
            // If the next event in tree was successfully loaded return true, that it's good.
            if (nextEventInTree())
                return true;
            // If the first event in this tree was not good, continue to the next tree
            else
                continue;
        }
        // If looping over all trees, we fail to find first event that's good,
        // return false and call it quits.
        // At this point it will exit the loop without processing any events.
        return false;
    }
    // If tree exists, it means that we're in the middle of a loop
    else
    {
        // If next event is successfully loaded proceed.
        if (nextEventInTree())
        {
            return true;
        }
        // If next event is not loaded then check why.
        else
        {
            // If failed because it was the last event in the whole chain to process, exit the loop.
            // You're done!
            if (allEventsInChainProcessed())
            {
                return false;
            }
            // If failed because it's last in the tree then load the next tree and the event
            else if (allEventsInTreeProcessed())
            {
                // Load the next tree if it returns true, then proceed to next event in tree.
                while (nextTree())
                {
                    // If the next event in tree was successfully loaded return true, that it's good.
                    if (nextEventInTree())
                        return true;
                    // If the first event in this tree was not good, continue to the next tree
                    else
                        continue;
                }
                // If looping over all trees, we fail to find first event that's good,
                // return false and call it quits.
                // Again you're done!
                return false;
            }
            else
            {
                // Why are you even here?
                // spit error and return false to avoid warnings
                error("You should not be here! Please contact philip@physics.ucsd.edu", __FUNCTION__);
                return false;
            }
        }
    }
}

//_________________________________________________________________________________________________
template <class TREECLASS>
void TasUtil::Looper<TREECLASS>::setFileList()
{
    if (!fileIter)
    {
        listOfFiles = tchain->GetListOfFiles();
        fileIter = new TObjArrayIter(listOfFiles);
    }
}

//_________________________________________________________________________________________________
template <class TREECLASS>
void TasUtil::Looper<TREECLASS>::setNEventsToProcess()
{
    if (tchain)
    {
        nEventsTotalInChain = tchain->GetEntries();
        if (nEventsToProcess < 0)
            nEventsToProcess = nEventsTotalInChain;
    }
}

#endif
