/*
 Macro to fit invariant mass distributions
 1) input files:
    masshistoLctopK0sPbPbCen[010,3050]data1_[prob_value]_[18r,18q].root
 2) output file (containig h_raw_signal_prob[prob], h_invmass_[ptmin]_[ptmax]_prob[prob], h_invmasstot_[ptmin]_[ptmax]_prob[prob]):
    raw_yields_[010,3050].root
 
 .L merge_and_fit_histo.C
 merge_and_fit_histo()
 */

#include "AliHFInvMassFitter.h"
/////// COMMON ////////
Double_t mass = 2.2864;
Int_t signal_fit_f = 0;    // kGaus=0, k2Gaus=1, k2GausSigmaRatioPar=2
Int_t background_fit_f = 2;// kExpo=0, kLin=1, kPol2=2, kNoBk=3, kPow=4, kPowEx=5
Int_t rebin = 6;

bool fixsigma = false;

bool signalOnly = true;

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
const std::string inputFile("../LckINT7HighMultwithJets/vAN-20190810_ROOT6-1/pp_mc/prodD2H/resultsMBvspt/masshisto.root");
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
const std::vector<double> sigmas = {0.005,0.0051,0.006,0.009,0.009,0.01,0.011,0.014};
// Use a given number of RMS of the histo to fit
bool fitDouble = false;
//bool fitDouble = true;
//const bool useRMSFitRange = true;
const bool useRMSFitRange = false;
//const std::vector<double> nRMSForFitRange = {3., 3., 3., 3., 3., 3., 3., 3.};
const std::vector<double> nRMSForFitRange = {2., 2., 2., 2., 2., 2., 2., 2.};
const double rmsStep = 0.2;
const int rmsScanSteps = 11;
const double rmsStart = 2;

double gausFitSingle(double* x, double* par)
{
    return par[0] / TMath::Sqrt(2 * TMath::Pi()) / par[2] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]);
}

double gausFitDouble(double* x, double* par)
{
    return par[0] / TMath::Sqrt(2 * TMath::Pi()) / par[2] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2])
           + par[3] / TMath::Sqrt(2 * TMath::Pi()) / par[5] * TMath::Exp(-(x[0] - par[4]) * (x[0] - par[4]) / 2. / par[5] / par[5]);
}


void mass_fit()
{
    std::vector<std::vector<double>> extractedSigmas(nPtBins);
    extractedSigmas.resize(nPtBins);
    TFile file(inputFile.c_str(), "READ");
    for(unsigned int i = 0; i < nPtBins; i++) {
        for(unsigned int j = 0; j < nSecondBins; j++) {
            const std::string histName = histNamePrefix + std::to_string(ptBins[i].first) + "_" + std::to_string(ptBins[i].second) + "_" + probsString[i] + histNameMiddle + "_" + secondBins[j].first + "_" + secondBins[j].second;
            TH1F* hmerge = dynamic_cast<TH1F*>(file.Get(histName.c_str()));
            if(!hmerge) {
                std::cerr << "Cannot find histogram " << histName << std::endl;
                exit(1);
            }
            hmerge->Rebin(rebin);
            Float_t bin_width = hmerge->GetXaxis()->GetBinWidth(3);
            TString histo_title=Form("p_{T}/ (GeV/c) #in [%.0u,%.0u], prob_{ML} > %s; #it{M} (GeV/#it{c}^{2}); Counts/%.0f MeV/#it{c}^{2}", ptBins[i].first,
                                                                                                                                  ptBins[i].second,
                                                                                                                                  probsString[i].c_str(),
                                                                                                                                  bin_width*1000.);
            /*TString histo_title=Form("p_{T}: %.0u-%.0u, %s: %s-%s, prob>%s; #it{M} (GeV/#it{c}^{2}); Counts/%.0f MeV/#it{c}^{2}", ptBins[i].first,
                                                                                                                                  ptBins[i].second,
                                                                                                                                  histNameMiddle.c_str(),
                                                                                                                                  secondBins[j].first.c_str(),
                                                                                                                                  secondBins[j].second.c_str(),
                                                                                                                                  probsString[i].c_str(),
                                                                                                                                  bin_width*1000.);*/
            TCanvas canvas("name","name", 500, 500);
            canvas.SetTicks();
            canvas.SetLogy();
            std::string canvasSave = histName;
            hmerge->SetTitle(histo_title.Data());
            auto rms = hmerge->GetRMS();
            auto histMax = hmerge->GetMaximum();
            if(signalOnly) {
                std::vector<double> sigmaChi2s;
                canvasSave += "_signalOnlyFit";
                TVirtualFitter::SetDefaultFitter("Minuit");
                TF1* signalFit = nullptr;
                if(fitDouble) {
                    signalFit = new TF1("signalFit", &gausFitDouble, massRanges[i].first, massRanges[i].second, 6);
                    signalFit->SetParName(0, "Int1");
                    signalFit->SetParName(1, "Mean1");
                    signalFit->SetParName(2, "Sigma1");
                    signalFit->SetParName(3, "Int2");
                    signalFit->SetParName(4, "Mean2");
                    signalFit->SetParName(5, "Sigma2");
                    signalFit->SetParameter(0, rms * histMax);
                    signalFit->SetParameter(1, means[i]);
                    signalFit->SetParameter(2, sigmas[i]);
                    signalFit->SetParameter(3, 0.01 * rms * histMax);
                    signalFit->SetParameter(4, means[i]);
                    signalFit->SetParameter(5, 4 * sigmas[i]);
                } else {
                    signalFit = new TF1("signalFit", &gausFitSingle, massRanges[i].first, massRanges[i].second, 3);
                    signalFit->SetParName(0, "Int");
                    signalFit->SetParName(1, "Mean");
                    signalFit->SetParName(2, "Sigma");
                    signalFit->SetParameter(0, rms * histMax);
                    signalFit->SetParameter(1, means[i]);
                    signalFit->SetParameter(2, sigmas[i]);
                }
                
                if(useRMSFitRange) {
                    hmerge->Fit("signalFit", "R,M,L,E,0", "", means[i] - nRMSForFitRange[i] * rms, means[i] + nRMSForFitRange[i] * rms);
                } else {
                    hmerge->Fit("signalFit", "R,M,L,E,0");
                }

                /*
                for(int n = 0; n < rmsScanSteps; n++) {
                    signalFit->SetParameters(0, rms * hmerge->GetMaximum());
                    signalFit->SetParameter(1, means[i]);
                    signalFit->SetParameter(2, sigmas[i]);
                    hmerge->Fit("signalFit", "M,L,E,0", "", means[i] - (rmsStart + n * rmsStep) * rms, means[i] + (rmsStart + n * rmsStep) * rms);
                    auto chi2 = signalFit->GetChisquare();
                    if(chi2 <= 0) {
                        chi2 = DBL_MAX;
                    }
                    sigmaChi2s.push_back(chi2);
                }
                std::vector<double>::iterator it = std::min_element(std::begin(sigmaChi2s), std::end(sigmaChi2s));
                int pos = std::distance(std::begin(sigmaChi2s), it);
                sigmaChi2s.clear();
                signalFit->SetParameters(0, rms * hmerge->GetMaximum());
                signalFit->SetParameter(1, means[i]);
                signalFit->SetParameter(2, sigmas[i]);
                hmerge->Fit("signalFit", "M,L,E,0", "", means[i] - (rmsStart + pos * rmsStep) * rms, means[i] + (rmsStart + pos * rmsStep) * rms);
                */
                extractedSigmas[i].push_back(signalFit->GetParameter(2));
                canvas.cd();
                hmerge->SetMarkerStyle(20);
                hmerge->SetMarkerSize(1);
                hmerge->SetStats(false);
                hmerge->Draw("PE");
                signalFit->SetLineColor(kBlue);
                signalFit->Draw("same");
                TPaveText* infoBox = new TPaveText(0.12, 0.6, 0.5, 0.85, "NDC");
                infoBox->SetTextSize(0.03);
                infoBox->SetBorderSize(0);
                infoBox->SetFillStyle(0);
                infoBox->SetTextAlign(11);
                for(int i = 0; i < signalFit->GetNpar(); i++) {
                    infoBox->AddText(Form("%s: %f", signalFit->GetParName(i), signalFit->GetParameter(i)));
                }
                infoBox->AddText(Form("#chi^{2} / NDF of fit: %f", signalFit->GetChisquare() / signalFit->GetNDF()));
                infoBox->Draw();
                canvas.Update();
                if(useRMSFitRange) {
                    infoBox->AddText(Form("fit in #RMS range: %.1f", nRMSForFitRange[i]));
                    TLine* lineLeft = new TLine(means[i] - nRMSForFitRange[i] * rms, 0, means[i] - nRMSForFitRange[i] * rms, hmerge->GetMaximum());
                    TLine* lineRight = new TLine(means[i] + nRMSForFitRange[i] * rms, 0, means[i] + nRMSForFitRange[i] * rms, hmerge->GetMaximum());
                    lineLeft->SetLineColor(kRed);
                    lineLeft->SetLineWidth(1);
                    lineLeft->SetLineStyle(7);
                    lineLeft->Draw();
                    lineRight->SetLineColor(kRed);
                    lineRight->SetLineWidth(1);
                    lineRight->SetLineStyle(7);
                    lineRight->Draw();
                }

            } else {
                AliHFInvMassFitter *fitter = new AliHFInvMassFitter(hmerge,massRanges[i].first,massRanges[i].second,background_fit_f,signal_fit_f);
                fitter->SetNSigma4SideBands(2.);
                fitter->SetUseLikelihoodFit();
                fitter->SetInitialGaussianMean(mass);
                fitter->SetInitialGaussianSigma(sigmas[i]);
                if(fixsigma)fitter->SetFixGaussianSigma(sigmas[i]);
                Bool_t out=fitter->MassFitter(0);
                if(!out) {
                    fitter->GetHistoClone()->Draw();
                }
                fitter->DrawHere(&canvas);
                extractedSigmas[i].push_back(fitter->GetSigma());
            }
            canvasSave += ".eps";
            canvas.SaveAs(canvasSave.c_str());
        }
    }
    file.Close();
    std::cout << "The sigmas after the fit are:\n";
    for(unsigned int i = 0; i < nPtBins; i++) {
        for(unsigned int j = 0; j < nSecondBins; j++) {
            std::cout << "p_T: " << ptBins[i].first << " - " << ptBins[i].second << "; " << histNameMiddle << ": " << secondBins[j].first << "-" << secondBins[j].second << ": " << extractedSigmas[i][j] << "\n";
        }
    }
}

