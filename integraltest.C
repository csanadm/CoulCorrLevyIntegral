#include <iostream>
#include <fstream>
#include <limits>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLatex.h>

using namespace std;

struct DataEntry {
    int NMaxIter;
    int LogTolerance;
    double Integral;
    int NFuncCalls;
};

void readDataFromFile(const char* inputfile, const char* outputfile) {
    ifstream infile(inputfile);
    if (!infile.is_open()) {
        cerr << "Error: Unable to open the input file." << endl;
        return;
    }

    TFile* file = TFile::Open(outputfile, "RECREATE");
    TTree* tree = new TTree("tree", "Data Tree");

    DataEntry data;
    tree->Branch("NMaxIter", &data.NMaxIter, "NMaxIter/I");
    tree->Branch("LogTolerance", &data.LogTolerance, "LogTolerance/I");
    tree->Branch("Integral", &data.Integral, "Integral/D");
    tree->Branch("NFuncCalls", &data.NFuncCalls, "NFuncCalls/I");

    while (infile >> data.NMaxIter >> data.LogTolerance >> data.Integral >> data.NFuncCalls) {
        tree->Fill();
    }

    tree->Write();
    file->Close();
}

void createMatrixPlot(const char* datafile, const char* outputfile) {
    TFile* file = TFile::Open(datafile);
    if (!file || file->IsZombie()) {
        cerr << "Error: Unable to open the input file." << endl;
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("tree"));
    if (!tree) {
        cerr << "Error: Tree 'tree' not found in the file." << endl;
        return;
    }

    DataEntry data;
    tree->SetBranchAddress("NMaxIter", &data.NMaxIter);
    tree->SetBranchAddress("LogTolerance", &data.LogTolerance);
    tree->SetBranchAddress("Integral", &data.Integral);
    tree->SetBranchAddress("NFuncCalls", &data.NFuncCalls);

    double minNMaxIter = tree->GetMinimum("NMaxIter")-0.5;
    double maxNMaxIter = tree->GetMaximum("NMaxIter")+0.5;
    double minLogTolerance = tree->GetMinimum("LogTolerance")-0.5;
    double maxLogTolerance = tree->GetMaximum("LogTolerance")+0.5;

    TH2D* hIntegral = new TH2D("hIntegral", "", maxLogTolerance - minLogTolerance, minLogTolerance, maxLogTolerance, maxNMaxIter - minNMaxIter, minNMaxIter, maxNMaxIter);
    TH2D* hNFuncCalls = new TH2D("hNFuncCalls", "", maxLogTolerance - minLogTolerance, minLogTolerance, maxLogTolerance, maxNMaxIter - minNMaxIter, minNMaxIter, maxNMaxIter);

    Long64_t numEntries = tree->GetEntries();
    for (Long64_t i = 0; i < numEntries; ++i) {
        tree->GetEntry(i);
        hIntegral->SetBinContent(hIntegral->FindBin(data.LogTolerance, data.NMaxIter), data.Integral);
        hNFuncCalls->SetBinContent(hNFuncCalls->FindBin(data.LogTolerance, data.NMaxIter), data.NFuncCalls);
    }

    TCanvas* canvas = new TCanvas("canvas", "Matrix Plots", 1200, 600);
    canvas->Divide(2, 1);

    canvas->cd(1);
    hIntegral->SetTitle("; -log_{10}(tolerance); # of max iterations");
    hIntegral->SetStats(false);
    hIntegral->Draw("colz");

    TLatex* latex1 = new TLatex();
    latex1->SetTextAlign(12);
    latex1->SetTextSize(0.04);
    latex1->SetNDC();
    latex1->DrawLatex(0.15, 0.9, "C(Q)-1 for");
    latex1->DrawLatex(0.15, 0.85, Form("#lambda=%0.1f, R=%0.1f fm, #alpha=%0.1f, Q=%0.2f GeV/c", 1.0, 6.0, 1.4, 0.20));

    TLegend* colorBar = new TLegend(0.0, 0.1, 0.12, 0.9);
    colorBar->SetFillColor(kWhite);
    colorBar->SetLineColor(kWhite);
    colorBar->SetShadowColor(kWhite);
    colorBar->SetFillStyle(0);
    colorBar->SetBorderSize(0);
    colorBar->SetTextAlign(12);
    colorBar->SetTextSize(0.03);
    colorBar->AddEntry("", "Min", "");
    colorBar->AddEntry("", "Max", "");
    colorBar->Draw();

    canvas->cd(2);
    hNFuncCalls->SetTitle("; -log_{10}(tolerance); # of max iterations");
    hNFuncCalls->SetStats(false);
    hNFuncCalls->Draw("colz");

    TLatex* latex2 = new TLatex();
    latex2->SetTextAlign(12);
    latex2->SetTextSize(0.04);
    latex2->SetNDC();
    latex2->DrawLatex(0.15, 0.9, "# of function calls for");
    latex2->DrawLatex(0.15, 0.85, Form("#lambda=%0.1f, R=%0.1f fm, #alpha=%0.1f, Q=%0.2f GeV/c", 1.0, 6.0, 1.4, 0.20));

    canvas->SaveAs(outputfile);

    delete hIntegral;
    delete hNFuncCalls;
    delete latex1;
    delete latex2;
    delete canvas;
    file->Close();
}

int integraltest() {
    const char* inputfile = "integraltest.Q200.out";
    const char* rootfile = "integraltest.Q200.root";
    const char* imagefile = "integraltest.Q200.pdf";

    readDataFromFile(inputfile, rootfile);

    createMatrixPlot(rootfile, imagefile);

    return 0;
}
