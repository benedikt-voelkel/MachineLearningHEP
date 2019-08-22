//0 1-2
//1 2-3
//2 3-4

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

void PlotSystRawYield(Int_t iPtBin=0, TString bkgTreat="",TString esesel="Largeq2"){
    
    for(unsigned int i = 0; i < nPtBins; i++) {
        for(unsigned int j = 0; j < nSecondBins; j++) {

            TString inputFilename = Form("RawSyst_pt_%u_%u_%s_%s_%s.root", ptBins[i].first, 
                                                                           ptBins[i].second, 
                                                                           histNameMiddle.c_str(), 
                                                                           secondBins[j].first.c_str(), 
                                                                           secondBins[j].second.c_str());
            // aka fil1
            TFile inputTrialsFile(inputFilename.Data(), "READ");


            /*
            TString filenameref="../D0_yeldRatios_3050.root";
            TFile *fref=new TFile(filenameref.Data());
            TH1D    *hRawYieldsUnb=(TH1D*)fref->Get("hRawYieldsUnbiased");
            TH1D    *hRawYieldsESESel=(TH1D*)fref->Get(Form("hRawYields%s",esesel.Data()));
            TH1D    *hRawYieldsRatio_notNorm=(TH1D*)hRawYieldsESESel->Clone("hRawYieldsRatioNotNotm");
            hRawYieldsRatio_notNorm->Divide(hRawYieldsUnb);
            Double_t ratioref=hRawYieldsRatio_notNorm->GetBinContent(iPtBin+1);
            */
            
            
            const Int_t nBackFuncCases=1;
            const Int_t nConfigCases=6;
            Int_t colorBC0=kGreen+2;
            Int_t colorBC1=kAzure-8;
            Int_t minBCrange=3;
            Int_t maxBCrange=5;
            Int_t nBCranges=maxBCrange-minBCrange+1;
            TString confCase[nConfigCases]={"FixedS","FixedSp20","FixedSm20","FreeS","FixedMeanFreeS","FixedMeanFixedS"};
            TString bkgFunc[nBackFuncCases]={"Pol2"};
            //TString bkgFunc[nBackFuncCases]={"Expo","Lin","Pol2","Pol3","Pol4","Pol5"};
            
            const Int_t totCases=nConfigCases*nBackFuncCases;
            
            // 0= not used; 1 = used for fit; 2= used also for bin count0, 3=use also bin count1, 4=use both binc
            Int_t mask[totCases]={0,0,1,0,0,0,   // fixed sigma (Expo, Lin, Pol2,Pol3,Pol4)
                0,0,0,0,0,0,   // fixed sigma upper
                0,0,0,0,0,0,   // fixed sigma lower
                0,0,0,0,0,0,   // free sigma, free mean
                0,0,0,0,0,0,   // free sigma, fixed mean
                0,0,0,0,0,0,   // fixed mean, fixed sigma
            };
            
            // Require chi2 cut and ignore everything above
            Double_t chi2Cut=2.;
            
            // Extracting all required histograms
            TH1F* histo1[totCases];
            Int_t jh=0;
            std::cout << "nconfigcases " << nConfigCases << "\t nbackgfunccases " << nBackFuncCases << std::endl;
            for(Int_t iConf=0; iConf<nConfigCases; iConf++) {
                for(Int_t iType=0; iType<nBackFuncCases; iType++) {
                    histo10[jh++]=(TH1F*)inputTrialsFile.Get(Form("hRawYieldTrial%s%s%s",bkgFunc[iType].Data(),confCase[iConf].Data(),bkgTreat.Data()));
                    if (!histo10[jh]) {
                        std::cerr << "Histo10 " 
                                  << Form("hRawYieldTrial%s%s%s",bkgFunc[iType].Data(),confCase[iConf].Data(),bkgTreat.Data()) 
                                  << " not found \n";
                    }
                }
            }
            
            /*
            TH1F* histo10[totCases];
            Int_t jhi=0;
            for(Int_t iConf=0; iConf<nConfigCases; iConf++){
                for(Int_t iType=0; iType<nBackFuncCases; iType++){
                    histo1[jhi++]=(TH1F*)fil1->Get(Form("hRawYieldTrial%s%s%s",bkgFunc[iType].Data(),confCase[iConf].Data(),bkgTreat.Data()));
                    if(!histo1[jhi])cout<<" Histo1 " << Form("hRawYieldTrial%s%s%s",bkgFunc[iType].Data(),confCase[iConf].Data(),bkgTreat.Data()) << " not found " <<endl;
                
                }
            }
            */
            
            Int_t totTrials=0;
            Int_t totTrialsBC0=0;
            Int_t totTrialsBC1=0;
            Int_t totHistos=0;
            Int_t first[totCases];
            Int_t last[totCases];
            Int_t firstBC0[totCases];
            Int_t lastBC0[totCases];
            Int_t firstBC1[totCases];
            Int_t lastBC1[totCases];
            Double_t minyd=9e9;
            Double_t maxyd=0;
            for(Int_t j=0; j<totCases; j++){
                if(mask[j]){
                    if(histo1[j]){
                        first[j]=totTrials;
                        totTrials+=histo1[j]->GetNbinsX();
                        last[j]=totTrials;
                        Double_t thisMin=histo1[j]->GetBinContent(histo1[j]->GetMinimumBin());
                        Double_t thisMax=histo1[j]->GetBinContent(histo1[j]->GetMaximumBin());
                        if(thisMin<minyd) minyd=thisMin;
                        if(thisMax>maxyd) maxyd=thisMax;
                        ++totHistos;
                        if(mask[j]==2 || mask[j]==4){
                            TString hbcname=histo1[j]->GetName();
                            hbcname.ReplaceAll("Trial","TrialBinC0");
                            //cout<< " name bc " << hbcname.Data() << endl;
                            TH2F* hbc2dt=(TH2F*)fil1->Get(hbcname.Data());
                            Int_t bnx,bny,bnz;
                            hbc2dt->GetMinimumBin(bnx,bny,bnz);
                            thisMin=hbc2dt->GetBinContent(bnx,bny);
                            hbc2dt->GetMaximumBin(bnx,bny,bnz);
                            thisMax=hbc2dt->GetBinContent(bnx,bny);
                            if(thisMin<minyd) minyd=thisMin;
                            if(thisMax>maxyd) maxyd=thisMax;
                            firstBC0[j]=totTrialsBC0;
                            totTrialsBC0+=hbc2dt->GetNbinsX();
                            lastBC0[j]=totTrialsBC0;
                        }
                        if(mask[j]==3 || mask[j]==4){
                            TString hbcname=histo1[j]->GetName();
                            hbcname.ReplaceAll("Trial","TrialBinC1");
                            TH2F* hbc2dt=(TH2F*)fil1->Get(hbcname.Data());
                            Int_t bnx,bny,bnz;
                            hbc2dt->GetMinimumBin(bnx,bny,bnz);
                            thisMin=hbc2dt->GetBinContent(bnx,bny);
                            hbc2dt->GetMaximumBin(bnx,bny,bnz);
                            thisMax=hbc2dt->GetBinContent(bnx,bny);
                            if(thisMin<minyd) minyd=thisMin;
                            if(thisMax>maxyd) maxyd=thisMax;
                            firstBC1[j]=totTrialsBC1;
                            totTrialsBC1+=hbc2dt->GetNbinsX();
                            lastBC1[j]=totTrialsBC1;
                        }
                    }else{
                        mask[j]=0;
                    }
                }
            }
            TLine **vlines = new TLine*[totCases];
            TLatex **tlabels = new TLatex*[totCases+1];
    
            printf("Histos merged = %d    totTrials=%d\n",totHistos,totTrials);
            for(Int_t ja=0; ja<totCases; ja++){
                if(mask[ja]){
                    printf("  %d) %s  -- %d \n",ja,histo1[ja]->GetName(),first[ja]);
                    vlines[ja]=new TLine(last[ja],0.,last[ja],50000.);
                    vlines[ja]->SetLineColor(kMagenta+2);
                    vlines[ja]->SetLineStyle(2);
                    TString ttt=histo1[ja]->GetName();
                    ttt.ReplaceAll("hRawYieldTrial","");
                    if(ttt.Contains("FixedMean")) ttt="Fix #mu";
                    if(ttt.Contains("FixedSp20")) ttt="#sigma+";
                    if(ttt.Contains("FixedSm20")) ttt="#sigma-";
                    if(ttt.Contains("FreeS")) ttt="Free #sigma";
                    ttt.ReplaceAll("FixedS","");
                    if(bkgTreat!="" && ttt.Contains(bkgTreat.Data())) ttt.ReplaceAll(bkgTreat.Data(),"");
                    tlabels[ja]=new TLatex(first[ja]+0.02*totTrials,10,ttt.Data());
                    tlabels[ja]->SetTextColor(kMagenta+2);
                    tlabels[ja]->SetTextColor(kMagenta+2);
                }
            }
            tlabels[totCases]=new TLatex(totTrials+30,10,"BinCnt");
            tlabels[totCases]->SetTextColor(kMagenta+2);
            tlabels[totCases]->SetTextColor(kMagenta+2);
            
            Printf("tottrials \t tottrialsBC0 \t \t tottrialsBC1 \t nBC ranges \t ");
            Printf("%d \t %d \t %d \t %d",totTrials,totTrialsBC0,totTrialsBC1,nBCranges);
   
    TH1F* hRawYieldAll=new TH1F("hRawYieldAll"," ; Trial # ; q_{2}-sel / unbiased",totTrials+totTrialsBC0*nBCranges,0.,totTrials+totTrialsBC0*nBCranges);
    TH1F* hRawYieldAllBC0=new TH1F("hRawYieldAllBC0"," ; Trial # ; q_{2}-sel / unbiased",totTrialsBC0*nBCranges,totTrials,totTrials+totTrialsBC0*nBCranges);
    TH1F* hRawYieldAllBC1=new TH1F("hRawYieldAllBC1"," ; Trial # ; q_{2}-sel / unbiased",totTrialsBC1*nBCranges,totTrials,totTrials+totTrialsBC1*nBCranges);
    //  TH1F* hRawYieldAll=new TH1F("hRawYieldAll"," ; Trial # ; Raw Yield",totTrials+totTrialsBC0+totTrialsBC1,0.3,1.5);
    //  TH1F* hRawYieldAllBC0=new TH1F("hRawYieldAllBC0"," ; Trial # ; Raw Yield",totTrialsBC0*nBCranges,1.5,10.5);
    //  TH1F* hRawYieldAllBC1=new TH1F("hRawYieldAllBC1"," ; Trial # ; Raw Yield",totTrialsBC1*nBCranges,1.5,10.5);
    TH1F* hMeanAll=new TH1F("hMeanAll"," ; Trial # ; Gaussian mean",totTrials,0.,totTrials);
    TH1F* hSigmaAll=new TH1F("hSigmaAll"," ; Trial # ; Gaussian #sigma",totTrials,0,totTrials);
    TH1F* hChi2All=new TH1F("hChi2All"," ; Trial # ; #chi^{2}",totTrials,0.,totTrials);
    hMeanAll->SetMarkerColor(kBlue+1);
    hSigmaAll->SetMarkerColor(kBlue+1);
    hChi2All->SetMarkerColor(kBlue+1);
    hMeanAll->SetLineColor(kBlue+1);
    hSigmaAll->SetLineColor(kBlue+1);
    hChi2All->SetLineColor(kBlue+1);
    TH1F* hMeanAll6=new TH1F("hMeanAll6"," ; Trial # ; Gaussian mean",totTrials,0.,totTrials);
    TH1F* hSigmaAll6=new TH1F("hSigmaAll6"," ; Trial # ; Gaussian #sigma",totTrials,0,totTrials);
    TH1F* hChi2All6=new TH1F("hChi2All6"," ; Trial # ; #chi^{2}",totTrials,0.,totTrials);
    hMeanAll6->SetMarkerColor(kRed+1);
    hSigmaAll6->SetMarkerColor(kRed+1);
    hChi2All6->SetMarkerColor(kRed+1);
    hMeanAll6->SetLineColor(kRed+1);
    hSigmaAll6->SetLineColor(kRed+1);
    hChi2All6->SetLineColor(kRed+1);
    TH1F* hRawYieldDistAll=new TH1F("hRawYieldDistAll","  ; q_{2}-sel / unbiased",300,0.,1.5);
    //hRawYieldDistAll->GetXaxis()->SetRangeUser(0.2,0.6);
    hRawYieldDistAll->SetFillStyle(3003);
    hRawYieldDistAll->SetFillColor(kBlue+1);
    TH1F* hRawYieldDistAllBC0=new TH1F("hRawYieldDistAllBC0","  ; q_{2}-sel / unbiased",300,0.,1.5);
    //hRawYieldDistAllBC0->GetXaxis()->SetRangeUser(0.2,0.6);
    hRawYieldDistAllBC0->SetFillStyle(3004);
    TH1F* hRawYieldDistAllBC1=new TH1F("hRawYieldDistAllBC1","  ; q_{2}-sel / unbiased",300,0,1.5);
    //hRawYieldDistAllBC1->GetXaxis()->SetRangeUser(0.2,0.6);
    hRawYieldDistAllBC1->SetFillStyle(3005);
    TH1F* hStatErrDistAll=new TH1F("hStatErrDistAll","  ; Stat Unc on Yield",300,0,1.5);
    TH1F* hRelStatErrDistAll=new TH1F("hRelStatErrDistAll","  ; Rel Stat Unc on Yield",100,0.,1.);
    
    Double_t minYield=999999.;
    Double_t maxYield=0.;
    Double_t sumy[4]={0.,0.,0.,0.};
    Double_t sumwei[4]={0.,0.,0.,0.};
    Double_t sumerr[4]={0.,0.,0.,0.};
    Double_t counts=0.;
    Double_t wei[4];
    Double_t maxFilled=-1;
    Printf("trial \t first \t last \t firstBC0 \t lastBC0 \t firstBC1 \t lastBC1");
    for(Int_t j=0; j<totCases; j++){
        if(mask[j]){
            Printf("%d \t %d \t %d \t %d \t %d \t %d \t %d",j,first[j],last[j],firstBC0[j],lastBC0[j],firstBC1[j],lastBC1[j]);
            
            TString hmeanname=histo1[j]->GetName();
            hmeanname.ReplaceAll("RawYield","Mean");
            TH1F* hmeant=(TH1F*)fil1->Get(hmeanname.Data());
            //TH1F* hmeant6=(TH1F*)fil6->Get(hmeanname.Data());
            
            TString hsigmaname=histo1[j]->GetName();
            hsigmaname.ReplaceAll("RawYield","Sigma");
            TH1F* hsigmat=(TH1F*)fil1->Get(hsigmaname.Data());
            //TH1F* hsigmat6=(TH1F*)fil6->Get(hsigmaname.Data());
            
            TString hchi2name=histo1[j]->GetName();
            hchi2name.ReplaceAll("RawYield","Chi2");
            TH1F* hchi2t=(TH1F*)fil1->Get(hchi2name.Data());
            //TH1F* hchi2t6=(TH1F*)fil6->Get(hchi2name.Data());
            
            TString hbcname=histo1[j]->GetName();
            hbcname.ReplaceAll("Trial","TrialBinC0");
            TH2F* hbc2dt010=(TH2F*)fil1->Get(hbcname.Data());
            //TH2F* hbc2dt060=(TH2F*)fil6->Get(hbcname.Data());
            //TH2F* hbc2dt0=(TH2F*)hbc2dt010->Clone(hbcname.Data());
            //hbc2dt0->Divide(hbc2dt010,hbc2dt060,1,1,"B");
            hbcname.ReplaceAll("BinC0","BinC1");
            TH2F* hbc2dt10=(TH2F*)fil1->Get(hbcname.Data());
            //TH2F* hbc2dt60=(TH2F*)fil6->Get(hbcname.Data());
            //TH2F* hbc2dt1=(TH2F*)hbc2dt10->Clone(hbcname.Data());
            //hbc2dt1->Divide(hbc2dt10,hbc2dt60,1,1,"B");
            for(Int_t ib=1; ib<=histo1[j]->GetNbinsX(); ib++){
                Double_t ry=histo1[j]->GetBinContent(ib);
                //cout<< " ry " << ry <<endl;
                Double_t ery=histo1[j]->GetBinError(ib);
                
                Double_t pos=hmeant->GetBinContent(ib);
                Double_t epos=hmeant->GetBinError(ib);
                
                Double_t sig=hsigmat->GetBinContent(ib);
                Double_t esig=hsigmat->GetBinError(ib);
                
                Double_t chi2=hchi2t->GetBinContent(ib);
                
                Double_t pos6=hmeant6->GetBinContent(ib);
                Double_t epos6=hmeant6->GetBinError(ib);
                
                Double_t sig6=hsigmat6->GetBinContent(ib);
                Double_t esig6=hsigmat6->GetBinError(ib);
                
                Double_t chi26=hchi2t->GetBinContent(ib);
                
                if(ry>0.001 && ery>(0.01*ry) && ery<(0.5*ry) && chi2<chi2Cut && chi26<chi2Cut){
                    hRawYieldDistAll->Fill(ry);
                    hStatErrDistAll->Fill(ery);
                    hRelStatErrDistAll->Fill(ery/ry);
                    hRawYieldAll->SetBinContent(first[j]+ib,ry);
                    hRawYieldAll->SetBinError(first[j]+ib,ery);
                    if(ry<minYield) minYield=ry;
                    if(ry>maxYield) maxYield=ry;
                    wei[0]=1.;
                    wei[1]=1./(ery*ery);
                    wei[2]=1./(ery*ery/(ry*ry));
                    wei[3]=1./(ery*ery/ry);
                    for(Int_t kw=0; kw<4; kw++){
                        sumy[kw]+=wei[kw]*ry;
                        sumerr[kw]+=wei[kw]*wei[kw]*ery*ery;
                        sumwei[kw]+=wei[kw];
                    }
                    counts+=1.;
                    hSigmaAll->SetBinContent(first[j]+ib,sig);
                    hSigmaAll->SetBinError(first[j]+ib,esig);
                    hMeanAll->SetBinContent(first[j]+ib,pos);
                    hMeanAll->SetBinError(first[j]+ib,epos);
                    hChi2All->SetBinContent(first[j]+ib,chi2);
                    hChi2All->SetBinError(first[j]+ib,0.000001);
                    hSigmaAll6->SetBinContent(first[j]+ib,sig6);
                    hSigmaAll6->SetBinError(first[j]+ib,esig6);
                    hMeanAll6->SetBinContent(first[j]+ib,pos6);
                    hMeanAll6->SetBinError(first[j]+ib,epos6);
                    hChi2All6->SetBinContent(first[j]+ib,chi26);
                    hChi2All6->SetBinError(first[j]+ib,0.000001);
                    if(mask[j]==2 || mask[j]==4){
                        for(Int_t iy=minBCrange; iy<=maxBCrange;iy++){
                            Double_t bc=hbc2dt0->GetBinContent(ib,iy);
                            Double_t ebc=hbc2dt0->GetBinError(ib,iy);
                            if(bc>0.001 && ebc<0.5*bc && bc<5.*ry){
                                Int_t theBin=iy+(firstBC0[j]+ib-1)*nBCranges;
                                cout<< " bin content " << bc << " the bin " << theBin << " BCrange " << iy << endl;
                                hRawYieldAllBC0->SetBinContent(theBin-2,bc);
                                hRawYieldAllBC0->SetBinError(theBin-2,ebc);
                                hRawYieldDistAllBC0->Fill(bc);
                                if(hRawYieldAllBC0->GetBinCenter(theBin-2)>maxFilled) maxFilled=hRawYieldAllBC0->GetBinCenter(theBin-2);
                            }
                        }
                    }
                    if(mask[j]==3 || mask[j]==4){
                        for(Int_t iy=minBCrange; iy<=maxBCrange;iy++){
                            Double_t bc=hbc2dt1->GetBinContent(ib,iy);
                            Double_t ebc=hbc2dt1->GetBinError(ib,iy);
                            if(bc>0.001 && ebc<0.5*bc && bc<5.*ry){
                                cout<< " bc " << bc << endl;
                                Int_t theBin=iy+(firstBC1[j]+ib-1)*nBCranges;
                                hRawYieldAllBC1->SetBinContent(theBin-2,bc);
                                hRawYieldAllBC1->SetBinError(theBin-2,ebc);
                                hRawYieldDistAllBC1->Fill(bc);
                                if(hRawYieldAllBC1->GetBinCenter(theBin-2)>maxFilled) maxFilled=hRawYieldAllBC1->GetBinCenter(theBin-2);
                            }
                        }
                    }
                }
            }
        }
    }
    
    Double_t weiav[4]={0.,0.,0.,0.};
    Double_t eweiav[4]={0.,0.,0.,0.};
    for(Int_t kw=0; kw<4; kw++){
        if(sumwei[kw]>0.){
            weiav[kw]=sumy[kw]/sumwei[kw];
            eweiav[kw]=TMath::Sqrt(sumerr[kw])/sumwei[kw];
        }
    }
    
    hRawYieldAll->SetStats(0);
    hMeanAll->SetStats(0);
    hSigmaAll->SetStats(0);
    hChi2All->SetStats(0);
    hChi2All->SetMarkerStyle(7);
    hMeanAll6->SetStats(0);
    hSigmaAll6->SetStats(0);
    hChi2All6->SetStats(0);
    hChi2All6->SetMarkerStyle(7);
    hMeanAll->SetMinimum(1.85);
    hMeanAll->SetMaximum(1.88);
    hSigmaAll->SetMinimum(0.);
    if(hSigmaAll->GetMaximum()<0.018) hSigmaAll->SetMaximum(0.018);
    hMeanAll6->SetMinimum(1.85);
    hMeanAll6->SetMaximum(1.88);
    hSigmaAll6->SetMinimum(0.);
    if(hSigmaAll6->GetMaximum()<0.018) hSigmaAll6->SetMaximum(0.018);
    hRawYieldAllBC0->SetStats(0);
    hRawYieldAllBC1->SetStats(0);
    // hRawYieldAllBC->SetMarkerStyle(7);
    // hRawYieldAllBC->SetMarkerSize(0.1);
    hRawYieldAllBC0->SetMarkerColor(colorBC0);
    hRawYieldAllBC0->SetLineColor(colorBC0);
    hRawYieldDistAllBC0->SetLineColor(colorBC0);
    hRawYieldDistAllBC0->SetFillColor(colorBC0);
    hRawYieldDistAllBC0->SetLineWidth(2);
    hRawYieldDistAllBC0->SetLineStyle(7);
    hRawYieldDistAllBC0->Scale(hRawYieldDistAll->GetEntries()/hRawYieldDistAllBC0->GetEntries());
    hRawYieldAllBC1->SetMarkerColor(colorBC1);
    hRawYieldAllBC1->SetLineColor(colorBC1);
    hRawYieldDistAllBC1->SetLineColor(colorBC1);
    hRawYieldDistAllBC1->SetFillColor(colorBC1);
    hRawYieldDistAllBC1->SetLineWidth(2);
    hRawYieldDistAllBC1->SetLineStyle(2);
    hRawYieldDistAllBC1->Scale(hRawYieldDistAll->GetEntries()/hRawYieldDistAllBC1->GetEntries());
    hRawYieldDistAll->SetLineWidth(2);
    
    
    TLine *l=new TLine(ratioref,0.,ratioref,hRawYieldDistAll->GetMaximum());
    l->SetLineColor(kRed);
    l->SetLineWidth(2);
    
    TLine *ll=new TLine(0.,ratioref,totTrials+totTrialsBC0*nBCranges,ratioref);
    ll->SetLineColor(kRed);
    ll->SetLineWidth(2);
    
    TCanvas* call=new TCanvas("call","All",1400,800);
    call->Divide(3,2);
    call->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.06);
    hSigmaAll->GetYaxis()->SetTitleOffset(1.7);
    hSigmaAll->Draw();
    hSigmaAll6->Draw("same");
    call->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.06);
    hMeanAll->GetYaxis()->SetTitleOffset(1.7);
    hMeanAll->Draw();
    hMeanAll6->Draw("same");
    call->cd(3);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.06);
    hChi2All->GetYaxis()->SetTitleOffset(1.7);
    hChi2All->Draw();
    hChi2All6->Draw("same");
    call->cd(4);
    hRawYieldAll->SetTitle(Form("p_{T} bin %d",iPtBin));
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.06);
    Double_t newmax=1.25*(hRawYieldAll->GetMaximum()+hRawYieldAll->GetBinError(1));
    hRawYieldAll->GetYaxis()->SetTitleOffset(1.7);
    hRawYieldAll->SetMaximum(newmax);
    if(maxFilled>0) hRawYieldAll->GetXaxis()->SetRangeUser(0.,maxFilled);
    hRawYieldAll->Draw();
    hRawYieldAllBC0->Draw("same");
    hRawYieldAllBC1->Draw("same");
    ll->Draw("same");
    TLatex* tweimean[4];
    for(Int_t kw=0; kw<4; kw++){
        tweimean[kw]=new TLatex(0.16,0.84-0.06*kw,Form("<Yield>_{wei%d} = %.1f #pm %.1f\n",kw,weiav[kw],eweiav[kw]*sqrt(counts)));
        tweimean[kw]->SetNDC();
        tweimean[kw]->SetTextColor(4);
        //    tweimean[kw]->Draw();
    }
    
    for(Int_t j=0; j<totCases; j++){
        if(mask[j]){
            vlines[j]->SetY2(newmax);
            vlines[j]->Draw("same");
            tlabels[j]->SetY(0.05*newmax);
            tlabels[j]->Draw();
        }
    }
    tlabels[totCases]->SetY(0.05*newmax);
    tlabels[totCases]->Draw();
    
    call->cd(5);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    hRawYieldDistAll->SetTitle(Form("p_{T} bin %d",iPtBin));
    hRawYieldDistAll->Draw();
    hRawYieldDistAll->GetXaxis()->SetRangeUser(minYield*0.8,maxYield*1.2);
    hRawYieldDistAllBC0->Draw("sameshist");
    hRawYieldDistAllBC1->Draw("sameshist");
    l->Draw("same");
    gPad->Update();
    TPaveStats* st=(TPaveStats*)hRawYieldDistAll->GetListOfFunctions()->FindObject("stats");
    st->SetY1NDC(0.71);
    st->SetY2NDC(0.9);
    TPaveStats* stb0=(TPaveStats*)hRawYieldDistAllBC0->GetListOfFunctions()->FindObject("stats");
    stb0->SetY1NDC(0.51);
    stb0->SetY2NDC(0.7);
    stb0->SetTextColor(hRawYieldDistAllBC0->GetLineColor());
    TPaveStats* stb1=(TPaveStats*)hRawYieldDistAllBC1->GetListOfFunctions()->FindObject("stats");
    stb1->SetY1NDC(0.31);
    stb1->SetY2NDC(0.5);
    stb1->SetTextColor(hRawYieldDistAllBC1->GetLineColor());
    Double_t perc[3]={0.15,0.5,0.85}; // quantiles for +-1 sigma
    Double_t lim70[3];
    hRawYieldDistAll->GetQuantiles(3,lim70,perc);
    call->cd(6);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    Double_t aver=hRawYieldDistAll->GetMean();
    TLatex* tmean=new TLatex(0.15,0.93,Form("mean=%.3f",aver));
    tmean->SetNDC();
    tmean->Draw();
    TLatex* tmedian=new TLatex(0.15,0.86,Form("median=%.3f",lim70[1]));
    tmedian->SetNDC();
    tmedian->Draw();
    Double_t averBC0=hRawYieldDistAllBC0->GetMean();
    TLatex* tmeanBC0=new TLatex(0.15,0.79,Form("mean(BinCount0)=%.3f",averBC0));
    tmeanBC0->SetNDC();
    tmeanBC0->SetTextColor(hRawYieldDistAllBC0->GetLineColor());
    tmeanBC0->Draw();
    Double_t averBC1=hRawYieldDistAllBC1->GetMean();
    TLatex* tmeanBC1=new TLatex(0.15,0.72,Form("mean(BinCount1)=%.3f",averBC1));
    tmeanBC1->SetNDC();
    tmeanBC1->SetTextColor(hRawYieldDistAllBC1->GetLineColor());
    tmeanBC1->Draw();
    Double_t val=hRawYieldDistAll->GetRMS();
    TLatex* thrms=new TLatex(0.15,0.62,Form("rms=%.3f  (%.2f%%)",val,val/aver*100.));
    thrms->SetNDC();
    thrms->Draw();
    val=hRawYieldDistAllBC0->GetRMS();
    TLatex* thrmsBC0=new TLatex(0.15,0.55,Form("rms(BinCount0)=%.3f  (%.2f%%)",val,val/averBC0*100.));
    thrmsBC0->SetNDC();
    thrmsBC0->SetTextColor(hRawYieldDistAllBC0->GetLineColor());
    thrmsBC0->Draw();
    val=hRawYieldDistAllBC1->GetRMS();
    TLatex* thrmsBC1=new TLatex(0.15,0.48,Form("rms(BinCount1)=%.3f  (%.2f%%)",val,val/averBC1*100.));
    thrmsBC1->SetNDC();
    thrmsBC1->SetTextColor(hRawYieldDistAllBC1->GetLineColor());
    thrmsBC1->Draw();
    TLatex* tmin=new TLatex(0.15,0.38,Form("min=%.3f      max=%.2f",minYield,maxYield));
    tmin->SetNDC();
    tmin->Draw();
    val=(maxYield-minYield)/sqrt(12);
    TLatex* trms=new TLatex(0.15,0.31,Form("(max-min)/sqrt(12)=%.3f  (%.2f%%)",val,val/aver*100.));
    trms->SetNDC();
    trms->Draw();
    TLatex* meanRef=new TLatex(0.15,0.24,Form("mean(ref)=%.3f",ratioref));
    meanRef->SetNDC();
    meanRef->SetTextColor(kRed);
    meanRef->Draw();
    TLatex* meanRefDiff=new TLatex(0.15,0.17,Form("mean(ref)-mean(fit)=%.3f  (%.2f%%)",ratioref-aver,100.*(ratioref-aver)/ratioref));
    meanRefDiff->SetNDC();
    meanRefDiff->SetTextColor(kBlack);
    meanRefDiff->Draw();
    TLatex* meanRefDiffBC0=new TLatex(0.15,0.10,Form("mean(ref)-mean(BC0)=%.3f  (%.2f%%)",ratioref-averBC0,100.*(ratioref-averBC0)/ratioref));
    meanRefDiffBC0->SetNDC();
    meanRefDiffBC0->SetTextColor(hRawYieldDistAllBC0->GetLineColor());
    meanRefDiffBC0->Draw();
    TLatex* meanRefDiffBC1=new TLatex(0.15,0.03,Form("mean(ref)-mean(BC1)=%.3f  (%.2f%%)",ratioref-averBC1,100.*(ratioref-averBC1)/ratioref));
    meanRefDiffBC1->SetNDC();
    meanRefDiffBC1->SetTextColor(hRawYieldDistAllBC1->GetLineColor());
    meanRefDiffBC1->Draw();
    /*val=(maxYield-aver)/sqrt(3);
    TLatex* tup=new TLatex(0.15,0.39,Form("(max-mean)/sqrt(3)=%.3f  (%.2f%%)",val,val/aver*100.));
    tup->SetNDC();
    tup->Draw();
    val=(aver-minYield)/sqrt(3);
    TLatex* tdw=new TLatex(0.15,0.32,Form("(mean-min)/sqrt(3)=%.3f  (%.2f%%)",val,val/aver*100.));
    tdw->SetNDC();
    tdw->Draw();
    TLatex* tl15=new TLatex(0.15,0.22,Form("15 percentile=%.3f",lim70[0]));
    tl15->SetNDC();
    tl15->Draw();
    TLatex* tl85=new TLatex(0.15,0.13,Form("85 percentile=%.3f",lim70[2]));
    tl85->SetNDC();
    tl85->Draw();
    val=(lim70[2]-lim70[0])/2.;
    TLatex* t1s=new TLatex(0.15,0.06,Form("70%% range =%.3f  (%.2f%%)",val,val/aver*100.));
    t1s->SetNDC();
    t1s->Draw();*/
    call->SaveAs(Form("%s/6pad_ptbin%d_%s.eps",esesel.Data(),iPtBin,esesel.Data()));
    for(Int_t kw=0; kw<4; kw++){
        printf("Weight %d: %.1f +- %.1f(stat) +- %.1f (syst)\n",kw,weiav[kw],eweiav[kw]*sqrt(counts),(maxYield-minYield)/sqrt(12));
    }
    
    
    
    
    TCanvas* c3p=new TCanvas("c3p","3Pad",1400,400);
    c3p->Divide(3,1);
    c3p->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.06);
    hRawYieldAll->Draw("PFC");
    hRawYieldAllBC0->Draw("PFCsame");
    hRawYieldAllBC1->Draw("PFCsame");
    ll->Draw("same");
    for(Int_t j=0; j<totCases; j++){
        if(mask[j]){
            vlines[j]->Draw("same");
            tlabels[j]->Draw();
        }
    }
    tlabels[totCases]->Draw();
    
    c3p->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    hRawYieldDistAll->Draw();
    hRawYieldDistAllBC0->Draw("sameshist");
    hRawYieldDistAllBC1->Draw("sameshist");
    l->Draw("same");
    c3p->cd(3);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    tmean->Draw();
    tmeanBC0->Draw();
    tmeanBC1->Draw();
    tmedian->Draw();
    thrms->Draw();
    thrmsBC0->Draw();
    thrmsBC1->Draw();
    tmin->Draw();
    trms->Draw();
     meanRef->Draw();
    meanRefDiff->Draw();
     meanRefDiffBC0->Draw();
     meanRefDiffBC1->Draw();
    //tup->Draw();
    //tdw->Draw();
    //tl15->Draw();
    //tl85->Draw();
    //t1s->Draw();
    c3p->SaveAs(Form("%s/3pad_ptbin%d_%s.eps",esesel.Data(),iPtBin,esesel.Data()));

    //   c3p->SaveAs(Form("MultiTrial_3pad.eps",iPtBin));
    //   c3p->SaveAs(Form("MultiTrial_3pad.gif",iPtBin));
    
    TCanvas* c2p=new TCanvas("c2p","2Pad",933,400);
    c2p->Divide(2,1);
    c2p->cd(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    hRawYieldDistAll->Draw();
    hRawYieldDistAllBC0->Draw("sameshist");
    hRawYieldDistAllBC1->Draw("sameshist");
    l->Draw("same");
    c2p->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    tmean->Draw();
    tmeanBC0->Draw();
    tmeanBC1->Draw();
    tmedian->Draw();
    thrms->Draw();
    thrmsBC0->Draw();
    thrmsBC1->Draw();
    tmin->Draw();
    trms->Draw();
    meanRef->Draw();
    meanRefDiff->Draw();
    meanRefDiffBC0->Draw();
    meanRefDiffBC1->Draw();
    //tup->Draw();
    //tdw->Draw();
    //tl15->Draw();
    //tl85->Draw();
    //t1s->Draw();
    //   c2p->SaveAs(Form("MultiTrial_2pad.eps",iPtBin));
    c2p->SaveAs(Form("%s/2pad_ptbin%d_%s.eps",esesel.Data(),iPtBin,esesel.Data()));

    
    //  TString outn = filnam1.Data();
    //  outn.Prepend("Comb");
    //  cout << outn.Data() << endl;
    //  TFile *outf = new TFile(outn.Data(), "RECREATE");
    ///  call->Write();
    //  outf->Write();
    // outf->Close();
    //   TString callname = outn.Data();
    //   callname.ReplaceAll(".root",".pdf");
    //   call->Print(callname);
    //   callname.ReplaceAll(".root",".png");
    //   call->Print(callname);
}

