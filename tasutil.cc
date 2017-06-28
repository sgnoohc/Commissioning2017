//  .
// ..: P. Chang, philip@physics.ucsd.edu

// N.B.: I kept a very compact style in listing the functions in order to easily move around with {,} in vim.
//       Therefore, if one wants to add new feature please make the function extremely short as possible.
//       If the function is too long, maybe it doesn't belong here!

#include "tasutil.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Particle class implementation
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

//__________________________________________________________________________________________________
TasUtil::Particle::Particle()
{
}

//__________________________________________________________________________________________________
TasUtil::Particle::~Particle()
{
}

//__________________________________________________________________________________________________
bool TasUtil::Particle::varExists(TString name)
{
    if (std::find(var_Nam.begin(), var_Nam.end(), name) != var_Nam.end())
        return true;
    return false;
}

//__________________________________________________________________________________________________
void TasUtil::Particle::createFltVar(TString name, TString title)
{
    if (varExists(name)) printAlreadyExistsError(name, __FUNCTION__, kFLT);
    var_Nam.push_back(name);
    var_Ttl[name] = title;
    var_Typ[name] = kFLT;
    var_Flt[name] = -999;
}

//__________________________________________________________________________________________________
void TasUtil::Particle::createIntVar(TString name, TString title)
{
    if (varExists(name)) printAlreadyExistsError(name, __FUNCTION__, kINT);
    var_Nam.push_back(name);
    var_Ttl[name] = title;
    var_Typ[name] = kINT;
    var_Int[name] = -999;
}

//__________________________________________________________________________________________________
void TasUtil::Particle::createTLVVar(TString name, TString title)
{
    if (varExists(name)) printAlreadyExistsError(name, __FUNCTION__, kTLV);
    var_Nam.push_back(name);
    var_Ttl[name] = title;
    var_Typ[name] = kTLV;
    var_TLV[name].SetXYZT(0, 0, 0, 0);
}

//__________________________________________________________________________________________________
void TasUtil::Particle::createStrVar(TString name, TString title)
{
    if (varExists(name)) printAlreadyExistsError(name, __FUNCTION__, kSTR);
    var_Nam.push_back(name);
    var_Ttl[name] = title;
    var_Typ[name] = kSTR;
    var_Str[name] = "";
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printError(TString name, TString msg, TString action, TasUtil::Particle::VarType type)
{
    std::cerr << "ERROR - " << msg.Data() << ", ";
    std::cerr << "name = " << name.Data() << ", ";
    std::cerr << "what was tried = " << action.Data() << ", ";
    switch (type)
    {
        case kINT: std::cerr << "type = int"            << " "; break;
        case kFLT: std::cerr << "type = float"          << " "; break;
        case kTLV: std::cerr << "type = TLorentzVector" << " "; break;
        case kSTR: std::cerr << "type = TString"        << " "; break;
        default: break;
    }
    std::cerr << std::endl;
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printAlreadyExistsError(TString name, TString action, TasUtil::Particle::VarType type)
{
    printError(name, "already exists!", action, type);
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printOutOfRangeError(TString name, TString action, TasUtil::Particle::VarType type)
{
    printError(name, "out-of-range error", action, type);
}

//__________________________________________________________________________________________________
void TasUtil::Particle::setFltVar(TString name, float var)
{
    try { var_Flt.at(name) = var; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kFLT); exit(1); }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::setIntVar(TString name, int var)
{
    try { var_Int.at(name) = var; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kINT); exit(1); }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::setStrVar(TString name, TString var)
{
    try { var_Str.at(name) = var; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kSTR); exit(1); }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::setTLVVar(TString name, TLorentzVector var)
{
    try { var_TLV.at(name) = var; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kTLV); exit(1); }
}

//__________________________________________________________________________________________________
const float& TasUtil::Particle::getFltVar(TString name) const
{
    try { return var_Flt.at(name); }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kFLT); exit(2); }
    return var_Flt.at(name);
}

//__________________________________________________________________________________________________
const int& TasUtil::Particle::getIntVar(TString name) const
{
    try { return var_Int.at(name); }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kINT); exit(2); }
    return var_Int.at(name);
}

//__________________________________________________________________________________________________
const TLorentzVector& TasUtil::Particle::getTLVVar(TString name) const
{
    try { return var_TLV.at(name); }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kTLV); exit(2); }
    return var_TLV.at(name);
}

//__________________________________________________________________________________________________
const TString& TasUtil::Particle::getStrVar(TString name) const
{
    try { return var_Str.at(name); }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kSTR); exit(2); }
    return var_Str.at(name);
}

//__________________________________________________________________________________________________
TasUtil::Particle::VarType TasUtil::Particle::getType(TString name)
{
    try { return var_Typ.at(name); }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__); exit(4); }
    return var_Typ.at(name);
}

//__________________________________________________________________________________________________
void TasUtil::Particle::print()
{
    for (auto& name : var_Nam)
    {
        switch (getType(name))
        {
            case kFLT: printFltVar(name); break;
            case kINT: printIntVar(name); break;
            case kTLV: printTLVVar(name); break;
            case kSTR: printStrVar(name); break;
            case kNON: printOutOfRangeError(name, __FUNCTION__); exit(5); break;
        }
    }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printFltVar(TString name)
{
    try { std::cout << name.Data() << " kFLT " << " " << var_Flt.at(name) << std::endl; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kFLT); exit(3); }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printIntVar(TString name)
{
    try { std::cout << name.Data() << " kINT " << " " << var_Int.at(name) << std::endl; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kINT); exit(3); }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printTLVVar(TString name)
{
    try {
        std::cout << name.Data() << " kTLV ";
        std::cout << var_TLV.at(name).Pt () << " ";
        std::cout << var_TLV.at(name).Eta() << " ";
        std::cout << var_TLV.at(name).Phi() << " ";
        std::cout << var_TLV.at(name).M  () << " ";
        std::cout << std::endl;
    }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kTLV); exit(3); }
}

//__________________________________________________________________________________________________
void TasUtil::Particle::printStrVar(TString name)
{
    try { std::cout << name.Data() << " kSTR " << " " << var_Str.at(name) << std::endl; }
    catch (const std::out_of_range& oor) { printOutOfRangeError(name, __FUNCTION__, kSTR); exit(3); }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Data extractor
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////
TasUtil::Particles TasUtil::get(
        std::function<int()> size,
        std::function<bool(int)> pass,
        std::function<TasUtil::Particle(int)> extractor)
{
    // The return object
    Particles particles;
    for (unsigned int idx_in_vec_ttree = 0; idx_in_vec_ttree < (unsigned int) size(); ++idx_in_vec_ttree)
    {
        if (pass(idx_in_vec_ttree))
        {
            particles.push_back(extractor(idx_in_vec_ttree));
        }
    }
    return particles;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Math functions
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

//_________________________________________________________________________________________________
TLorentzVector TasUtil::Math::TLV(TasUtil::Math::LorentzVector p)
{
    TLorentzVector tlv;
    tlv.SetPxPyPzE(p.x(), p.y(), p.z(), p.e());
    return tlv;
}

//_________________________________________________________________________________________________
TasUtil::Math::LorentzVector TasUtil::Math::LV(TLorentzVector p)
{
    TasUtil::Math::LorentzVector lv;
    lv.SetPxPyPzE(p.X(), p.Y(), p.X(), p.E());
    return lv;
}

//_________________________________________________________________________________________________
TLorentzVector TasUtil::Math::TLVXYZE(float x, float y, float z, float e)
{
    TLorentzVector tlv;
    tlv.SetPxPyPzE(x, y, z, e);
    return tlv;
}

//_________________________________________________________________________________________________
TasUtil::Math::LorentzVector TasUtil::Math::LVXYZE (float x, float y, float z, float e)
{
    TasUtil::Math::LorentzVector lv;
    lv.SetPxPyPzE(x, y, z, e);
    return lv;
}

//_________________________________________________________________________________________________
TLorentzVector TasUtil::Math::TLVPtEtaPhiM(float x, float y, float z, float e)
{
    TLorentzVector tlv;
    tlv.SetPtEtaPhiM(x, y, z, e);
    return tlv;
}

//_________________________________________________________________________________________________
TLorentzVector TasUtil::Math::TLVPtEtaPhiE(float x, float y, float z, float e)
{
    TLorentzVector tlv;
    tlv.SetPtEtaPhiE(x, y, z, e);
    return tlv;
}

//_________________________________________________________________________________________________
float TasUtil::Math::MT(TLorentzVector p4, TLorentzVector metp4)
{
    /// Compute MT
    float mt = sqrt(2 * p4.Pt() * metp4.Pt() * (1 - cos(metp4.DeltaPhi(p4))));
    return mt;
}

//_________________________________________________________________________________________________
TLorentzVector TasUtil::Math::TLVMET(float met, float met_phi)
{
    TLorentzVector tlv;
    tlv.SetPtEtaPhiE(met, 0, met_phi, met);
    return tlv;
}

//_________________________________________________________________________________________________
TasUtil::Math::LorentzVector  TasUtil::Math::METLV(float met, float met_phi)
{
    // For some odd reason LorentzVector of float that SNT uses only work with SetPxPyPzE
    TLorentzVector tlv = TLVMET(met, met_phi);
    LorentzVector lv;
    lv.SetPxPyPzE(lv.x(), lv.y(), lv.z(), lv.e());
    return lv;
}

float TasUtil::Math::DEta(TLorentzVector a, TLorentzVector b) { return fabs(a.Eta() - b.Eta()); }
float TasUtil::Math::DPhi(TLorentzVector a, TLorentzVector b) { return fabs(a.DeltaPhi(b));     }
float TasUtil::Math::DR  (TLorentzVector a, TLorentzVector b) { return fabs(a.DeltaR(b));       }
float TasUtil::Math::DPt (TLorentzVector a, TLorentzVector b) { return fabs(a.Pt()-b.Pt());     }
float TasUtil::Math::Mass(TLorentzVector a, TLorentzVector b) { return (a + b).M();             }
float TasUtil::Math::Pt  (TLorentzVector a, TLorentzVector b) { return (a + b).Pt();            }
float TasUtil::Math::MT  (TLorentzVector p4  , float met, float met_phi) { return MT(p4, TLVMET(met, met_phi)); }

float TasUtil::Math::DEta(TasUtil::Math::LorentzVector a, float                        b) { return DEta(TLV(a), TLVPtEtaPhiE(1,b,0,1)); }
float TasUtil::Math::DPhi(TasUtil::Math::LorentzVector a, float                        b) { return DPhi(TLV(a), TLVPtEtaPhiE(1,0,b,1)); }
float TasUtil::Math::DEta(TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return DEta(TLV(a), TLV(b)               ); }
float TasUtil::Math::DPhi(TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return DPhi(TLV(a), TLV(b)               ); }
float TasUtil::Math::DR  (TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return DR  (TLV(a), TLV(b)               ); }
float TasUtil::Math::DPt (TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return DPt (TLV(a), TLV(b)               ); }
float TasUtil::Math::Mass(TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return Mass(TLV(a), TLV(b)               ); }
float TasUtil::Math::Pt  (TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return Pt  (TLV(a), TLV(b)               ); }
float TasUtil::Math::MT  (TasUtil::Math::LorentzVector a, TasUtil::Math::LorentzVector b) { return MT  (TLV(a), TLV(b)               ); }
float TasUtil::Math::MT  (TasUtil::Math::LorentzVector a, float met, float met_phi) { return MT(TLV(a), TLVMET(met, met_phi)); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Printing utilities
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

//_________________________________________________________________________________________________
void TasUtil::print(TString msg, const char* fname, int flush_before, int flush_after)
{
    /// printf replacement
    if (flush_before) printf("\n");
    if (strlen(fname) == 0) printf("TasUtil:: %s\n", msg.Data());
    else printf("TasUtil:: [in func %s()] %s\n", fname, msg.Data());
    if (flush_after) printf("\n");
    fflush(stdout);
}

//_________________________________________________________________________________________________
void TasUtil::warning(TString msg, const char* fname)
{
    /// warning message. Does not exit the program.
    print("WARNING - "+msg+".\n", fname);
}

//_________________________________________________________________________________________________
void TasUtil::error(TString msg, const char* fname, int is_error)
{
    /// Error & exit
    if (!is_error)
        return;
    //exit();
    print("ERROR - "+msg+", exiting.\n", fname);
    abort();
}

//_________________________________________________________________________________________________
void TasUtil::start(int q, int sleep_time)
{
    /// Fun start (from TM Hong's BaBar days)
    if (q)
        return;
    print(" _");
    print("/\\\\");
    print("\\ \\\\  \\__/ \\__/");
    print(" \\ \\\\ (oo) (oo)  Here we come!");
    print("  \\_\\\\/~~\\_/~~\\_");
    print(" _.-~===========~-._");
    print("(___________________)");
    print("      \\_______/");
    print();
    print("System info:");
    gSystem->Exec("hostname");
    gSystem->Exec("uname -a");
    gSystem->Exec("date");
    gSystem->Exec("whoami");
    gSystem->Exec("pwd");
    fflush(stdout);
    if (sleep_time>0)
        sleep(sleep_time);
}

//_________________________________________________________________________________________________
void TasUtil::announce(TString msg, Int_t q)
{
    /// Fun message presented by the moose
    if (q)
        return;
    // Random
    srand(time(NULL));      // Init rand seed
    Int_t   r = rand()%10+1;// Generate rand number
    Int_t   moose = r>4 ? 1 : 0;
    // Moose type
    TString eyes  = "open";
    if      (r==9) eyes = "closed";
    else if (r==8) eyes = "dead";
    else if (r==7) eyes = "small";
    else if (r==7) eyes = "sunny";
    else if (r==6) eyes = "calc";
    else if (r==4) eyes = "crazy";
    else if (r==3) eyes = "vampire";
    else if (r==2) eyes = "rich";
    else if (r==1) eyes = "sick";
    print();
    if      (msg.Length()>0) print("________________________________________");
    if      (msg.Length()>0) print(msg);
    if      (msg.Length()>0) print("--O-------------------------------------");
    if      (moose)          print("  O    \\_\\_    _/_/");
    if      (moose)          print("   O       \\__/");
    else                     print("   O       ^__^");
    if      (eyes=="open")   print("    o      (oo)\\_______");
    else if (eyes=="closed") print("    o      (==)\\_______");
    else if (eyes=="dead")   print("    o      (xx)\\_______");
    else if (eyes=="small")  print("    o      (..)\\_______");
    else if (eyes=="sunny")  print("    o      (66)\\_______");
    else if (eyes=="calc")   print("    o      (00)\\_______");
    else if (eyes=="crazy")  print("    o      (**)\\_______");
    else if (eyes=="vampire")print("    o      (@@)\\_______");
    else if (eyes=="rich")   print("    o      ($$)\\_______");
    else if (eyes=="sick")   print("    o      (++)\\_______");
    if (true)                print("     o     (__)\\       )\\/\\");
    if      (eyes=="dead")   print("            U  ||----w |");
    else if (eyes=="crazy")  print("            U  ||----w |");
    else if (eyes=="sick")   print("            U  ||----w |");
    else if (eyes=="vampire")print("            VV ||----w |");
    else                     print("               ||----w |");
    if (true)                print("               ||     ||");
    print();
    //sleep(0);
}

//_________________________________________________________________________________________________
void TasUtil::exit(int q)
{
    /// Fun exit (from TM Hong's BaBar days)
    if (q)
        return;
    print();
    print("   \\__/    \\__/");
    print("   (oo)    (oo)");
    print("(\\//~~\\\\  //~~\\\\");
    print(" \\/\\__//  \\\\__//\\");
    print("   ||||    ||\\\\ Who cares!");
    print("__ |||| __ |||| ___");
    print("  (_)(_)  (_)(_)");
    print();
}
//eof
