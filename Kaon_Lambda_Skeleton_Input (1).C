#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TGraph.h>
#include <TLine.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TStyle.h>

//// Add acceptance 5-45 degrees ////
//To get the angle of the kaon you use the TLorentzvector method called .Theta. After line 102 you type
//double angle=Kaon_Plus.Theta()*TMath::RadToPi(); Then you have an if statement for the angle
void Kaon_Lambda_Skeleton_Input(){
  TFile *fileOutput1= new TFile("./Kaon_Lambda_Output.root","recreate");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating histograms here

  // Angular dependence histograms
  TH2D* h_el_Scattered=new TH2D("h_el_Scattered","Angular dependence of scattered electron",200,0,11,200,0,180);
  TH2D* h_Kaon_Plus=new TH2D("h_Kaon_Plus","Angular dependence of K^{+}",200,0,11,200,0,180);
  TH2D* h_Lamba=new TH2D("h_Lamba","Angular dependence of lambda",200,0,11,200,0,180);
  TH2D* h_Lambda_Proton=new TH2D("h_Lambda_Proton","Angular dependence of p (#Lambda)",200,0,11,200,0,180);
  TH2D* h_Lambda_Pi_Minus=new TH2D("h_Lambda_Pi_Minus","Angular dependence of #pi^{-} (#Lambda)",200,0,11,200,0,180);

  // Invariant mass histograms
  TH1F* h_inv_Lambda=new TH1F("h_inv_Lambda","Invariant mass of #Lambda",10,0,2);
  TH1F* h_inv_Lambda_e_FD=new TH1F("h_inv_Lambda_e_FD","Invariant mass of #Lambda",10,0,2);


  // Acceptance histograms forward detector
  TH1F* h_acceptance_Lambda_e_FD=new TH1F("h_acceptance_Lambda_e_FD","Acceptance with e' in FD ",10,0,2);
  TH1F* h_acceptance_Lambda_eK_FD=new TH1F("h_acceptance_Lambda_eK_FD","Acceptance with e' and K in FD ",10,0,2);
  TH1F* h_acceptance_Lambda_eKp_FD=new TH1F("h_acceptance_Lambda_eKp_FD","Acceptance with e' , K and p in FD ",10,0,2);
  TH1F* h_acceptance_Lambda_eKpi_FD=new TH1F("h_acceptance_Lambda_eKpi_FD","Acceptance with e', K, and pi in FD",10,0,2);
  TH1F* h_acceptance_Lambda_eKppi_FD=new TH1F("h_acceptance_Lambda_eKppi_FD","Acceptance with e', K , p and pi",10,0,2);

  // Kaon event histograms
  //h_inv_Lambda_ek_FD h_Lambda_Proton_ek  h_Lambda_Pi_Minus_ek  h_Lambda_Proton_ekppi

  TH1F* h_inv_Lambda_ek_FD=new TH1F("h_inv_Lambda_ek_FD","Kaon Acceptance ",10,0,2); //for forward detector?
  TH1F* h_Lambda_Proton_ek=new TH1F("h_Lambda_Proton_ek","Acceptance for e k protons",10,0,2);
  TH1F* h_Lambda_Pi_Minus_ek=new TH1F("h_Lambda_Pi_Minus_ek","Acceptance of e k pions",10,0,2);
  TH1F* h_Lambda_Proton_ekppi=new TH1F("h_Lambda_Proton_ekppi","Acceptance of e k p pi",10,0,2);

  TH1F* h_Kaon_Angle=new TH1F("h_Kaon_Angle","Kaon Angle",100,-1,1);


  // central detector histograms

  //TH1F* h_inv_Lambda_ek_FD_Pi=new TH1F("h_inv_Lambda_ek_FD_Pi","Kaon Acceptance for central detector",10,0,2);
  TH1F* h_Numerator=new TH1F("h_Numerator","Numerator for; Case Number ",10,1,11); // case 6
  TH1F* h_Denominator=new TH1F("h_Denominator","Denomator for; Case Number;  ",10,1,11); // case 6

  TH1F* h_eFD_KFD_PiCD=new TH1F("h_eFD_KFD_PiC","Acceptance for ",10,0,2); // case 6
  TH1F* h_eFD_KFD_PiCD_PFD=new TH1F("h_eFD_KFD_PiCD_PFD","Acceptance for ",10,0,2); //case 7
  TH1F* h_eFD_KFD_PCD=new TH1F("h_eFD_KFD_PCD","Acceptance for ",10,0,2); //case 8
  TH1F* h_eFD_KFD_PCD_PiFD=new TH1F("h_eFD_KFD_PCD_PiFD","Acceptance for ",10,0,2); //case 9
  TH1F* h_eFD_KFD_PCD_PiCD=new TH1F("h_eFD_KFD_PCD_PiCD","Acceptance for ",10,0,2); //case 10

  TH1F* h_Lambda_Proton_Angle=new TH1F("h_Lambda_Proton_Angle","Lambda Proton Angle",100,-1,1);

  // Acceptance histograms central detector

  TH1F* h_acceptance_Lambda_Pion_CD=new TH1F("h_acceptance_Lambda_Pion_CD","Pion Acceptance CD",10,0,2);
  TH1F* h_acceptance_Lambda_PionCD_Proton_FD=new TH1F("h_acceptance_Lambda_PionCD_Proton_FD","Proton FD Pion Acceptance CD ",10,0,2);
  TH1F* h_acceptance_Lambda_Proton_CD=new TH1F("h_acceptance_Lambda_Proton_CD","Proton Acceptance CD",10,0,2);
  TH1F* h_acceptance_Lambda_ProtonCD_Pion_FD=new TH1F("h_acceptance_Lambda_ProtonCD_Pion_FD","Pion FD Proton Acceptance",10,0,2);
  TH1F* h_acceptance_Lambda_Proton_Pion_CD=new TH1F("h_acceptance_Lambda_Proton_Pion_CD","Proton Pion CD Acceptance",10,0,2);

  // 2d acceptance histograms

  TH1F* h_AcceptanceAll=new TH1F("h_AcceptanceAll","Acceptance; Case Number;  ",10,1,11);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating particles and terms

  // How many events to simulate and percentage completed
  Int_t nevents=100000;
  Int_t Percentage=nevents/100;

  // Creating TLorentzVectors of particles
  TLorentzVector Beam, Target; // Beam and target
  TLorentzVector *el_Scattered, *Lambda, *Kaon_Plus; // First vertex particles
  TLorentzVector *Lambda_Proton, *Lambda_Pi_Minus; // Second vertex particles
  // q is momentum transfer from beam to target, used to calculate q^2 weight
  TLorentzVector q;

  // Making Weights
  Double_t Phasespace_Weight_1;
  Double_t qWeight;

  // Setting TLorentzVectors for beam and target
  Beam.SetXYZM(0,0,10.6,0.00051);
  Target.SetXYZM(0,0,0,0.93827); // change to have deutoron mass (in new code)

  // Defining initial vertex and masses of particles
  TLorentzVector V1 = Beam + Target; // total energy
  Double_t Masses_1[3] = {0.00051,1.11568,0.49368}; // e^-, lambda, K^+  // new definition in new code = 6 variables = e , lambda, cascade , K + , K+ , K0
  Double_t Masses_2[2] = {0.93827,0.13957}; // p, pi^-
  TLorentzVector Inv_Lambda; // 4-vector of lambda from proton + pion

  // Creating decay vertices
  TGenPhaseSpace Vertex_1, Vertex_2;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over events

  // Looping over simulated events
  for (Long64_t i=0;i<nevents;){

    // Counter, shows percentage of simulated events completed
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }

    // SetDecay(total energy, no. particles, mass array)
    // Setting the first decay
    Vertex_1.SetDecay(V1,3,Masses_1); // 3 becomes 6 in new code

    // Defining the phasespace for this decay
    Phasespace_Weight_1=Vertex_1.Generate();

    // Assigning the decay particles to the correct mass from the array above
    el_Scattered=Vertex_1.GetDecay(0);
    Lambda=Vertex_1.GetDecay(1);
    Kaon_Plus=Vertex_1.GetDecay(2);

    // Defining q^2 weight
    q=Beam-(TLorentzVector)*el_Scattered; // Subtracting e^- out from e^- in
    qWeight=1/(q.Rho()*q.Rho()); // .Rho() gives momentum of the 4-vector


    TLorentzVector TargetProton(0,0,0, 0.93827);
    TLorentzVector CMToTal=q+TargetProton;
    TVector3 BoostVec=CMToTal.BoostVector();
    TLorentzVector BoostedKaon=(TLorentzVector)*Kaon_Plus;
    BoostedKaon.Boost(-BoostVec);
    double CMThetaAngleOfKaon=BoostedKaon.Theta();
    double PwaveWeight=1+3*TMath::Power(TMath::Cos(CMThetaAngleOfKaon),2);  // fill Kaon Angle
    //double Theta_Kaon = TMath::RadToDeg()CMThetaAngleOfKaon;
    // Getting total energy for second decay from lambda
    TLorentzVector V2 = (TLorentzVector)*Lambda;

    // Setting the second decay
    Vertex_2.SetDecay(V2,2,Masses_2);
    // Defining the phasespace for this decay
    Double_t Phasespace_Weight_2=Vertex_2.Generate();
    Lambda_Proton=Vertex_2.GetDecay(0);
    Lambda_Pi_Minus=Vertex_2.GetDecay(1);

    // Defining LambdaDecay weight

    TLorentzVector CMToTalLambda=(TLorentzVector)*Lambda_Proton +(TLorentzVector)*Lambda_Pi_Minus;

    TVector3 BoostVecLambda=CMToTalLambda.BoostVector();
    TLorentzVector BoostedProton=(TLorentzVector)*Lambda_Proton;
    BoostedProton.Boost(-BoostVecLambda);
    double LambdaAngleofProton=BoostedProton.Angle(Lambda->Vect());
    double LambdaDecayWeight=1+0.75*TMath::Cos(LambdaAngleofProton);


    // total weight
    double TotalWeight=Phasespace_Weight_1*Phasespace_Weight_2*qWeight; //*PwaveWeight*LambdaDecayWeight;

    h_Kaon_Angle->Fill(TMath::Cos(CMThetaAngleOfKaon),TotalWeight);
    h_Lambda_Proton_Angle->Fill(TMath::Cos(LambdaAngleofProton),TotalWeight);

    // Adding 4-vector of proton and pion, use to find invariant mass of lambda with .M()
    Inv_Lambda = (TLorentzVector)*Lambda_Proton + (TLorentzVector)*Lambda_Pi_Minus;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Filling histograms (x, weights applied)
    h_inv_Lambda->Fill(Inv_Lambda.M(),TotalWeight);
    for (int i=1;i<11;i++) h_Denominator->Fill(i,TotalWeight);

    // Setting acceptance cut for scattered electron
    if(el_Scattered->Theta()*TMath::RadToDeg()>5 && el_Scattered->Theta()*TMath::RadToDeg()<35){ // always 5 - 35
      h_Numerator->Fill(1,TotalWeight);
      // Filling angular distribution histograms (x,y,weights applied)
      h_el_Scattered->Fill(el_Scattered->Rho(),el_Scattered->Theta()*TMath::RadToDeg(),TotalWeight);
      h_Kaon_Plus->Fill(Lambda->Rho(),Lambda->Theta()*TMath::RadToDeg(),TotalWeight);
      h_Lamba->Fill(Kaon_Plus->Rho(),Kaon_Plus->Theta()*TMath::RadToDeg(),TotalWeight);
      h_Lambda_Proton->Fill(Lambda_Proton->Rho(),Lambda_Proton->Theta()*TMath::RadToDeg(), TotalWeight);
      h_Lambda_Pi_Minus->Fill(Lambda_Pi_Minus->Rho(),Lambda_Pi_Minus->Theta()*TMath::RadToDeg(), TotalWeight);

      // Filling invariant mass histogram
      h_inv_Lambda_e_FD->Fill(Inv_Lambda.M(), TotalWeight);

//case 2 Just kaon cut {}
//case 6 if kaon if pions
      if(Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<45
        && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
            h_eFD_KFD_PiCD->Fill(Inv_Lambda.M(), TotalWeight);
            h_Numerator->Fill(6,TotalWeight);
        }

//case 7 if kaon if pion if proton cd
      if(Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<45
      && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45
      && Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Lambda_Proton->Theta()*TMath::RadToDeg()<45){
        h_eFD_KFD_PiCD_PFD->Fill(Inv_Lambda.M(), TotalWeight);
        h_Numerator->Fill(7,TotalWeight);
      }

//case 8 if kaon if proton fd
      if(Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<45
      && Lambda_Proton->Theta()*TMath::RadToDeg()>45 && Lambda_Proton->Theta()*TMath::RadToDeg()<120){
        h_eFD_KFD_PCD->Fill(Inv_Lambda.M(), TotalWeight);
        h_Numerator->Fill(8,TotalWeight);
      }


//case 9
      if(Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<45
       && Lambda_Proton->Theta()*TMath::RadToDeg()>45 && Lambda_Proton->Theta()*TMath::RadToDeg()<120
       && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
//FILL HISTOGRAM here
        h_eFD_KFD_PCD_PiFD->Fill(Inv_Lambda.M(), TotalWeight);
        h_Numerator->Fill(9,TotalWeight);
      }

//case 10
      if(Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<45
      && Lambda_Proton->Theta()*TMath::RadToDeg()>45 && Lambda_Proton->Theta()*TMath::RadToDeg()<120
      && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>45 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<120){
        h_eFD_KFD_PCD_PiCD->Fill(Inv_Lambda.M(), TotalWeight);
        h_Numerator->Fill(10,TotalWeight);
      }


      //if your kaon falls between 5 and 45 degrees -> Fill another histograms  _*h_inv_Lambda_ek_FD
      if(Kaon_Plus->Theta()*TMath::RadToDeg()>5 && Kaon_Plus->Theta()*TMath::RadToDeg()<45){ //case 2
        h_inv_Lambda_ek_FD->Fill(Inv_Lambda.M(), TotalWeight);
        h_Numerator->Fill(2,TotalWeight);
        if(Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Lambda_Proton->Theta()*TMath::RadToDeg()<45){ // forward detector only //case 3
          //Acceptance for e k protons
          h_Lambda_Proton_ek->Fill(Inv_Lambda.M(), TotalWeight);
          h_Numerator->Fill(3,TotalWeight);
        }
        if(Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){ // forward detector only // case 4
          //acceptance of e k pions
          h_Lambda_Pi_Minus_ek->Fill(Inv_Lambda.M(), TotalWeight);
          h_Numerator->Fill(4,TotalWeight);

          if(Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Lambda_Proton->Theta()*TMath::RadToDeg()<45){ // forward detector only //case 5
            //acceptance of e k p pi
            h_Lambda_Proton_ekppi->Fill(Inv_Lambda.M(), TotalWeight);
            h_Numerator->Fill(5,TotalWeight);


          }
        }
      }
    }

    i++;

  }
  // Find acceptance by dividing histogram after angle cut by original histogram
  // -> Divide (A,B) gives A/B
  // acceptance with and without pweight **
  h_acceptance_Lambda_e_FD->Divide(h_inv_Lambda_e_FD,h_inv_Lambda);
  h_acceptance_Lambda_eK_FD->Divide(h_inv_Lambda_ek_FD,h_inv_Lambda);
  h_acceptance_Lambda_eKp_FD->Divide(h_Lambda_Proton_ek,h_inv_Lambda);
  h_acceptance_Lambda_eKpi_FD->Divide(h_Lambda_Pi_Minus_ek,h_inv_Lambda);
  h_acceptance_Lambda_eKppi_FD->Divide(h_Lambda_Proton_ekppi,h_inv_Lambda);
  h_acceptance_Lambda_Proton_CD->Divide(h_eFD_KFD_PCD,h_inv_Lambda);
  h_acceptance_Lambda_Proton_Pion_CD->Divide(h_eFD_KFD_PCD_PiCD,h_inv_Lambda);
  h_acceptance_Lambda_Pion_CD->Divide(h_eFD_KFD_PiCD,h_inv_Lambda);
  h_acceptance_Lambda_PionCD_Proton_FD->Divide(h_eFD_KFD_PiCD_PFD,h_inv_Lambda);
  h_acceptance_Lambda_ProtonCD_Pion_FD->Divide(h_eFD_KFD_PCD_PiFD,h_inv_Lambda);
  TCanvas *c1=new TCanvas("c1","c1",800,900);
  c1->Divide(2,2);
  c1->cd(1);
  h_acceptance_Lambda_e_FD->Draw();
  c1->cd(2);
  h_acceptance_Lambda_eK_FD->Draw();

  h_AcceptanceAll->Divide(h_Numerator, h_Denominator);
  TCanvas *c2=new TCanvas("c2","c2",800,900);
  h_AcceptanceAll->Draw();

  //print:(Theta_Kaon);
  fileOutput1->Write();
}
