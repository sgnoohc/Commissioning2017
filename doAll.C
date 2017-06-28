#include "ScanChain.C"

void doAll(TString filepath="", TString output_prefix="")
{
    if (filepath.IsNull()) return;
    TChain *ch = new TChain("Events");
    ch->Add(filepath);
    ScanChain(ch, output_prefix);
}
