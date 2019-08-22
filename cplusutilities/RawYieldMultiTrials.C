#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include "AliHFMassFitter.h"
#include "AliHFMultiTrials.h"
#endif

const std::string histNamePrefix("hmasspt_cand");

/////// Analyses ///////

/*
const std::string histNameMiddle("pt_jet");
const std::string inputFile("../LckINT7HighMultwithJets/vAN-20190810_ROOT6-1/pp_mc_prodD2H/resultsMBjetvspt/masshisto.root");
const unsigned int nPtBins = 6;
const unsigned int nSecondBins = 1;
const std::vector<std::pair<unsigned int,unsigned int>> ptBins = { {2,4}, {4,5}, {5,6}, {6,8}, {8,12}, {12,24} };
const std::vector<std::pair<std::string,std::string>> secondBins = { {"5.00", "15.00"} };
// Probabilities
const std::vector<std::string> probsString = {"0.40","0.40","0.40","0.40","0.30","0.30"};
// signal ranges
const std::vector<std::pair<double,double>> massRanges = { {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}};
// initial means and sigmas
const std::vector<double> means = {2.2864,2.2864,2.2864,2.2864,2.2864,2.2864};
const std::vector<double> sigmas = {0.006,0.009,0.009,0.01,0.011,0.014};
// Use a given number of RMS of the histo to fit
const bool useRMSFitRange = true;
const std::vector<double> nRMSForFitRange = {3., 3., 3., 3., 2., 2.};
*/
const std::string histNameMiddle("n_tracklets_corr");
const std::string inputFile("../LckINT7HighMultwithJets/vAN-20190810_ROOT6-1/pp_data/resultsMBvspt/masshisto.root");
const unsigned int nPtBins = 8;
const unsigned int nSecondBins = 1;
const std::vector<std::pair<unsigned int,unsigned int>> ptBins = { {1,2}, {2,3}, {3,4}, {4,5}, {5,6}, {6,8}, {8,12}, {12,24} };
const std::vector<std::pair<std::string,std::string>> secondBins = { {"0.00","9999.00"}, {"0.00", "30.00"}, {"30.00", "60.00"}, {"60.00", "100.00"} };
// Probabilities
const std::vector<std::string> probsString = {"0.40","0.40","0.40","0.40","0.40","0.40","0.30","0.30"};
// signal ranges
const std::vector<std::pair<double,double>> massRanges = { {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}, {2.14,2.436}};
// initial means and sigmas
const std::vector<double> means = {2.2864,2.2864,2.2864,2.2864,2.2864,2.2864,2.2864,2.2864};
const std::vector<double> sigmas = {0.007790, 0.007760, 0.008215, 0.009110, 0.009685, 0.010916, 0.013211, 0.018401};
// Use a given number of RMS of the histo to fit
bool fitDouble = false;
//bool fitDouble = true;
const bool useRMSFitRange = true;
//const bool useRMSFitRange = false;
const std::vector<double> nRMSForFitRange = {3., 3., 3., 3., 3., 3., 3., 3.};
//const std::vector<double> nRMSForFitRange = {2., 2., 2., 2., 2., 2., 2., 2.};
const double rmsStep = 0.2;
const int rmsScanSteps = 11;
const double rmsStart = 2;

//////// TRIAL SETTINGS ////////

const Int_t nRebinSteps=2;
Int_t rebinStep[nRebinSteps]={2,4};
const Int_t nMinMassSteps=4;
Double_t minMassStep[nMinMassSteps]={2.10, 2.14, 2.18};
const Int_t nMaxMassSteps=4;
Double_t maxMassStep[nMinMassSteps]={2.40, 2.436, 2.48};

const Int_t nStepsBC=3;
Double_t nSigmasBC[nStepsBC]={2.5,3.0,4.0};


void RawYieldMultiTrials() {
 
    TFile file(inputFile.c_str(), "READ");
    for(unsigned int i = 0; i < nPtBins; i++) {
        for(unsigned int j = 0; j < nSecondBins; j++) {
            const std::string histName = histNamePrefix + std::to_string(ptBins[i].first) + "_" + std::to_string(ptBins[i].second) + "_" + probsString[i] + histNameMiddle + "_" + secondBins[j].first + "_" + secondBins[j].second;
            TH1F* hmerge = dynamic_cast<TH1F*>(file.Get(histName.c_str()));
            if(!hmerge) {
                std::cerr << "Cannot find histogram " << histName << std::endl;
                exit(1);
            }

            AliHFMultiTrials* mt=new AliHFMultiTrials();
            mt->SetSuffixForHistoNames("");
            mt->SetMass(means[i]);
            mt->SetSigmaGaussMC(sigmas[i]);
            mt->SetUseExpoBackground(kTRUE);
            mt->SetUseLinBackground(kFALSE);
            mt->SetUsePol2Background(kTRUE);
            mt->SetUsePol3Background(kFALSE);
            mt->SetUsePol4Background(kFALSE);
            mt->SetUsePol5Background(kFALSE);
            mt->ConfigureRebinSteps(nRebinSteps,rebinStep);
            mt->ConfigureLowLimFitSteps(nMinMassSteps,minMassStep);
            mt->ConfigureUpLimFitSteps(nMaxMassSteps,maxMassStep);
            mt->SetDrawIndividualFits(kFALSE);

            TH1D hmergeAsDouble;
            hmerge->Copy(hmergeAsDouble);
            TCanvas canvasMassFit("canvasMassFit","MassFit");
            Bool_t success = mt->DoMultiTrials(&hmergeAsDouble, &canvasMassFit);

            if(success) {
                TCanvas canvasTrials("canvasTrials","All Trials");
                TString outfilnam=Form("RawSyst_pt_%u_%u_%s_%s_%s.root", ptBins[i].first, ptBins[i].second, histNameMiddle.c_str(), secondBins[j].first.c_str(), secondBins[j].second.c_str());

                mt->DrawHistos(&canvasTrials);
                mt->SaveToRoot(outfilnam.Data(),"recreate");
            }
        }
    }
}
