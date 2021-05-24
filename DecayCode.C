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
void DecayCode(){
  TFile *fileOutput1= new TFile("./Kaon_Decay_Output.root","recreate");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating histograms here
  // acceptance ratio histograms
  TH1F* h_Numerator=new TH1F("h_Numerator","Numerator for; Case Number ",10,1,11); //
  TH1F* h_Denominator=new TH1F("h_Denominator","Denomator for; Case Number;  ",10,1,11); //

  // Angular dependence histograms

  TH2D* h_el_Scattered=new TH2D("h_el_Scattered","Angular dependence of scattered electron",200,0,11,200,0,180);
  TH2D* h_Kaon_Plus1=new TH2D("h_Kaon_Plus1","Angular dependence of first K^{+}",200,0,11,200,0,180);
  TH2D* h_Kaon_Plus2=new TH2D("h_Kaon_Plus2","Angular dependence of second K^{+}",200,0,11,200,0,180);
  TH2D* h_Kaon0=new TH2D("h_Kaon0","Angular dependence of K^{0}",200,0,11,200,0,180);
  TH2D* h_Cascade=new TH2D("h_Cascade","Angular dependence of Cascade",200,0,11,200,0,180);
  TH2D* h_Lambda=new TH2D("h_Lambda","Angular dependence of Lambda",200,0,11,200,0,180);
  TH2D* h_Kaon0_Pion=new TH2D("h_Kaon0_Pion","Angular dependence of K^{0} Pion",200,0,11,200,0,180);
  TH2D* h_Kaon0_Pi_Minus=new TH2D("h_Kaon0_Pi_Minus","Angular dependence of K^{0} Pion^{-}",200,0,11,200,0,180);
  TH2D* h_Lambda_Proton=new TH2D("h_Lambda_Proton","Angular dependence of initial Lambda proton",200,0,11,200,0,180);
  TH2D* h_Lambda_Pi_Minus=new TH2D("h_Lambda_Pi_Minus","Angular dependence of initial Lambda pion^{-}",200,0,11,200,0,180);
  TH2D* h_Cascade_Lambda=new TH2D("h_Cascade_Lambda","Angular dependence of cascade lambda",200,0,11,200,0,180);
  TH2D* h_Cascade_Pi_Minus=new TH2D("h_Cascade_Pi_Minus","Angular dependence of cascade Pion^{-}",200,0,11,200,0,180);
  TH2D* h_Cascade_Lambda_Proton=new TH2D("h_Cascade_Lambda_Proton","Angular dependence of cascade lambda proton",200,0,11,200,0,180);
  TH2D* h_Cascade_Lambda_Pi_Minus=new TH2D("h_Cascade_Lambda_Pi_Minus","Angular dependence of cascade lambda Pion^{-}",200,0,11,200,0,180);

  // invariant Masses

  TH1F* h_inv_Lambda=new TH1F("h_inv_Lambda","Invariant mass of initial #Lambda",10,0,2);


  // Acceptances for Cascade, Kaon0, Cascade_Lambda, Cascade and Cascade_Lambda
  // all are eFD_PiMinusFD_PiPlusFD

  TH1F* h_eFD_PiMinusFD_PiPlusFD=new TH1F("h_eFD_PiMinusFD_PiPlusFD","Acceptance of e K+ K+ with Pion^{-}  Pion^{+} FD",10,0,2); // case 1
  TH1F* h_LProton_LPiMinus=new TH1F("h_LProton_LPiMinus","Acceptance of e Pion^{-} Pion^{+}, Lambda p Pion^{-}",10,0,2); // case 2
  TH1F* h_CLProton_CLPiMinus_CPiMinus=new TH1F("h_CLProton_CLPiMinus_CPiMinus","Acceptance of e Pion^{-} Pion^{+}, Cascade Lambda p and Pion^{-}, Cascade Pion^{-}",10,0,2);  // case 3
  TH1F* h_LProton_LPiMinus_CLProton_CLPiMinus=new TH1F("h_LProton_LPiMinus_CLProton_CLPiMinus","Acceptance of e Pion^{-} Pion^{+}, Lambda p and Pion^{-}, Cascade Lambda Pion^{-} and p",10,0,2); // case 4
  TH1F* h_LProton_LPiMinus_CLProton_CLPiMinus_CPiMinus=new TH1F("h_LProton_LPiMinus_CLProton_CLPiMinus_CPiMinus","Acceptance of e Pion^{-} Pion^{+}, Lambda p and Pion^{-}, Cascade Lambda Pion^{-} and p, and Cascade Pion^{-}",10,0,2); // case 5
  TH1F* h_eFD_KFD_NK0=new TH1F("h_eFD_KFD_NK0","Acceptance of e K+ K+ without Pion^{-} Pion^{+} in FD",10,0,2); // case 6
  TH1F* h_LProton_LPiMinus_NK0=new TH1F("h_LProton_LPiMinus_NK0","Acceptance of e K+ K+, Lambda Proton and Lambda Pion^{-}",10,0,2); // case 7
  TH1F* h_CPiMinus_CLProton_CLPiMinus_NK0=new TH1F("h_CPiMinus_CLProton_CLPiMinus_NK0","Acceptance of e K+ K+, Lambda p and Pion^{-}, Cascade Lambda p and Pion^{-}",10,0,2); // case 8
  TH1F* h_LProton_LPiMinus_CLProton_CLPiMinus_CPiMinus_NK0=new TH1F("h_LProton_LPiMinus_CLProton_CLPiMinus_CPiMinus_NK0","Acceptance of e K+ K+, Lambda p and Pion^{-}, Cascade Lambda p and Pion^{-} and Cascade Pion^{-}",10,0,2); // case 9
  TH1F* h_LProton_LPiMinus_CPiMinus_NK0=new TH1F("h_LProton_LPiMinus_CPiMinus_NK0","Acceptance of e K+ K+, Lambda p and Pion^{-} and Cascade Pion^{-}",10,0,2); // case 10

  TH1F* h_AcceptanceAll=new TH1F("h_AcceptanceAll","Acceptances of each case listed",10,1,11);

  // Angles

  TH1F* h_Lambda_Proton_Angle=new TH1F("h_Lambda_Proton_Angle","Lambda Proton Angle",100,-1,1);
  TH1F* h_Kaon_Angle=new TH1F("h_Kaon_Angle","Kaon Angle",100,-1,1);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating particles and terms

  // How many events to simulate and percentage completed
  Int_t nevents=100000;
  Int_t Percentage=nevents/100;

  // Creating TLorentzVectors of particles
  TLorentzVector Beam, Target; // Beam and target
  TLorentzVector *el_Scattered, *Cascade, *Kaon_Plus1, *Kaon_Plus2, *Kaon_0, *Lambda; // First vertex particles
  TLorentzVector *Kaon0_Pion , *Kaon0_Pi_Minus; // Second vertex particles
  TLorentzVector *Lambda_Proton, *Lambda_Pi_Minus; // third
  TLorentzVector *Cascade_Lambda, *Cascade_Pi_Minus; // fourth
  TLorentzVector *Cascade_Lambda_Proton, *Cascade_Lambda_Pi_Minus; // fifth

  // q is momentum transfer from beam to target, used to calculate q^2 weight
  TLorentzVector q;

  // Making Weights
  Double_t Phasespace_Weight_1;  // define other weights
  Double_t Phasespace_Weight_2;
  Double_t Phasespace_Weight_3;
  Double_t Phasespace_Weight_4;
  Double_t Phasespace_Weight_5;
  Double_t qWeight;

  // Setting TLorentzVectors for beam and target
  Beam.SetXYZM(0,0,10.6,0.00051);
  Target.SetXYZM(0,0,0,1.87516); // change to have deutoron mass (in new code)

  // Defining initial vertex and masses of particles
  TLorentzVector V1 = Beam + Target; // total energy

  Double_t Masses_Init[6] = {0.00051,0.493677,0.493677,0.49765,1.315,1.11568}; // new definition in new code = 6 variables = e , , K + , K+ , K0, cascade , lambda
  Double_t Masses_K0[2] = {0.13957,0.13957}; // pi^+, pi^- -> kaon0 decay
  Double_t Masses_Lambda[2] = {0.93827,0.13957}; // p, pi^- -> cascade_lambda decay p, pi^- -> cascade_lambda decay
  Double_t Masses_Cascade[2] = {1.11568,0.13957}; // lambda, pi^- -> cascade decay


  TLorentzVector Inv_Lambda; // 4-vector of lambda from proton + pion

  // Creating decay vertices
  TGenPhaseSpace Vertex_Initial, Vertex_K0, Vertex_Cascade, Vertex_Cascade_Lambda, Vertex_Lambda;  // edit these into new code below

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over events

  // Looping over simulated events
  for (Long64_t i=0;i<nevents;i++){

    // Counter, shows percentage of simulated events completed
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }

    // SetDecay(total energy, no. particles, mass array)
    // Setting the first decay
    Vertex_Initial.SetDecay(V1,6,Masses_Init); // 3 becomes 6 in new code

    // Defining the phasespace for this decay
    Phasespace_Weight_1=Vertex_Initial.Generate();


    // Assigning the decay particles to the correct mass from the array above
    el_Scattered=Vertex_Initial.GetDecay(0);
    Kaon_Plus1=Vertex_Initial.GetDecay(1);
    Kaon_Plus2=Vertex_Initial.GetDecay(2);
    Kaon_0=Vertex_Initial.GetDecay(3);
    Cascade=Vertex_Initial.GetDecay(4);
    Lambda=Vertex_Initial.GetDecay(5);

    // Defining q^2 weight
    q=Beam-(TLorentzVector)*el_Scattered; // Subtracting e^- out from e^- in
    qWeight=1/(q.Rho()*q.Rho()); // .Rho() gives momentum of the 4-vector

    // Getting total energy for second decay from kaon0
    TLorentzVector V2 = (TLorentzVector)*Kaon_0;

    // Setting the second decay -> Kaon_0 decay
    Vertex_K0.SetDecay(V2,2,Masses_K0);
    // Defining the phasespace for this decay
    Double_t Phasespace_Weight_2=Vertex_K0.Generate();
    Kaon0_Pion=Vertex_K0.GetDecay(0);
    Kaon0_Pi_Minus=Vertex_K0.GetDecay(1);


    // Getting total energy for fourth decay from cascade_lambda
    TLorentzVector V3 = (TLorentzVector)*Lambda;

    // 4th Decay -> lambda_1 decay
    Vertex_Lambda.SetDecay(V3,2,Masses_Lambda);
    Double_t Phasespace_Weight_3=Vertex_Lambda.Generate();
    Lambda_Proton=Vertex_Lambda.GetDecay(0);
    Lambda_Pi_Minus=Vertex_Lambda.GetDecay(1);

    // Getting total energy for fifth decay from Cascade
    TLorentzVector V4 = (TLorentzVector)*Cascade;

    // 5th Decay -> cascade
    Vertex_Cascade.SetDecay(V4,2,Masses_Cascade);
    Double_t Phasespace_Weight_4=Vertex_Cascade.Generate();
    Cascade_Lambda=Vertex_Cascade.GetDecay(0);
    Cascade_Pi_Minus=Vertex_Cascade.GetDecay(1);

    // Getting total energy for second decay from lambda
    TLorentzVector V5 = (TLorentzVector)*Cascade_Lambda;

    // 6th Decay -> lambda_2 decay -> cascade_lambda
    Vertex_Cascade_Lambda.SetDecay(V5,2,Masses_Lambda);
    Double_t Phasespace_Weight_5=Vertex_Cascade_Lambda.Generate();
    Cascade_Lambda_Proton=Vertex_Cascade_Lambda.GetDecay(0);
    Cascade_Lambda_Pi_Minus=Vertex_Cascade_Lambda.GetDecay(1);

    // Adding 4-vector of proton and pion, use to find invariant mass of lambda with .M()
    Inv_Lambda = (TLorentzVector)*Lambda_Proton + (TLorentzVector)*Lambda_Pi_Minus;
    //Inv_Kaon0 =
    //Inv_Cascade

    // PwaveWeight
    TLorentzVector TargetDeuteron(0,0,0, 1.87516); //
    TLorentzVector CMToTal=q+TargetDeuteron;
    TVector3 BoostVec=CMToTal.BoostVector();
    TLorentzVector BoostedKaon=(TLorentzVector)*Kaon_Plus1 + (TLorentzVector)*Kaon_Plus2 + (TLorentzVector)*Kaon_0; // which Kaon?
    BoostedKaon.Boost(-BoostVec);
    double CMThetaAngleOfKaons=BoostedKaon.Theta();
    double PwaveWeight=1+3*TMath::Power(TMath::Cos(CMThetaAngleOfKaons),2);  // fill Kaon Angle

    //Initial Lambda Decay Weight

    TLorentzVector CMToTalLambda=(TLorentzVector)*Lambda_Proton +(TLorentzVector)*Lambda_Pi_Minus;
    TVector3 BoostVecLambda=CMToTalLambda.BoostVector();
    TLorentzVector BoostedProton=(TLorentzVector)*Lambda_Proton;
    BoostedProton.Boost(-BoostVecLambda);
    double LambdaAngleofProton=BoostedProton.Angle(Lambda->Vect());
    double LambdaDecayWeight=1+0.75*TMath::Cos(LambdaAngleofProton);

    //Cascade_Minus Weight
    TLorentzVector CMToTalCascade=(TLorentzVector)*Cascade_Lambda +(TLorentzVector)*Cascade_Pi_Minus;
    TVector3 BoostVecCascade=CMToTalCascade.BoostVector();
    TLorentzVector BoostedLambda=(TLorentzVector)*Cascade_Lambda;
    BoostedLambda.Boost(-BoostVecCascade);
    double CascadeAngleofLambda=BoostedLambda.Angle(Cascade->Vect());
    double CascadeMinusDecayWeight=1-0.392*TMath::Cos(CascadeAngleofLambda);

    //Cascade Lambda Decay
    TLorentzVector CMToTalCascadeLambda=(TLorentzVector)*Cascade_Lambda_Proton +(TLorentzVector)*Cascade_Lambda_Pi_Minus;
    TVector3 BoostVecCascadeLambda=CMToTalCascadeLambda.BoostVector();
    TLorentzVector BoostedCascadeProton=(TLorentzVector)*Cascade_Lambda_Proton;
    BoostedCascadeProton.Boost(-BoostVecCascadeLambda);
    double CascadeLambdaAngleofProton=BoostedCascadeProton.Angle(Cascade_Lambda->Vect());
    double CascadeLambdaDecayWeight=1+0.75*TMath::Cos(CascadeLambdaAngleofProton);





    // total weight
    double TotalWeight=Phasespace_Weight_1*Phasespace_Weight_2*qWeight*Phasespace_Weight_3*Phasespace_Weight_4*Phasespace_Weight_5;//*CascadeMinusDecayWeight*LambdaDecayWeight*CascadeLambdaDecayWeight*PwaveWeight;

    //with and without pwaveweight
    // Proton momentum angle and Kaon momentum angle

    // Angle histograms

    h_Kaon_Angle->Fill(TMath::Cos(CMThetaAngleOfKaons),TotalWeight);
    h_Lambda_Proton_Angle->Fill(TMath::Cos(LambdaAngleofProton),TotalWeight);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Filling histograms (x, weights applied)
    // K+ K+ K0 in FD always -> case0 , rest can be alternated
    // write down all alternatives for FD and CD
    h_inv_Lambda->Fill(Inv_Lambda.M(),TotalWeight);
    for (int i=1;i<11;i++) h_Denominator->Fill(i,TotalWeight); // can be used in new code

    // check googledoc lab notes

    // Setting acceptance cut for scattered electron
    //if(el_Scattered->Theta()*TMath::RadToDeg()>5 && el_Scattered->Theta()*TMath::RadToDeg()<35){ // always 5 - 35 - forward detector
      // Filling angular distribution histograms (x,y,weights applied)
    h_el_Scattered->Fill(el_Scattered->Rho(),el_Scattered->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Kaon_Plus1->Fill(Kaon_Plus1->Rho(),Kaon_Plus1->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Kaon_Plus2->Fill(Kaon_Plus2->Rho(),Kaon_Plus2->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Kaon0->Fill(Kaon_0->Rho(),Kaon_0->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Lambda->Fill(Lambda->Rho(),Lambda->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Cascade->Fill(Cascade->Rho(),Cascade->Theta()*TMath::RadToDeg(), TotalWeight);
    h_Kaon0_Pion->Fill(Kaon0_Pion->Rho(),Kaon0_Pion->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Kaon0_Pi_Minus->Fill(Kaon0_Pi_Minus->Rho(),Kaon0_Pi_Minus->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Lambda_Proton->Fill(Lambda_Proton->Rho(),Lambda_Proton->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Lambda_Pi_Minus->Fill(Lambda_Pi_Minus->Rho(),Lambda_Pi_Minus->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Cascade_Lambda->Fill(Cascade_Lambda->Rho(),Cascade_Lambda->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Cascade_Pi_Minus->Fill(Cascade_Pi_Minus->Rho(),Cascade_Pi_Minus->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Cascade_Lambda_Proton->Fill(Cascade_Lambda_Proton->Rho(),Cascade_Lambda_Proton->Theta()*TMath::RadToDeg(),TotalWeight);
    h_Cascade_Lambda_Pi_Minus->Fill(Cascade_Lambda_Pi_Minus->Rho(),Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg(),TotalWeight);
    //}

    if(Kaon_Plus1->Theta()*TMath::RadToDeg()>5 && Kaon_Plus1->Theta()*TMath::RadToDeg()<45 && Kaon_Plus2->Theta()*TMath::RadToDeg()>5 && Kaon_Plus2->Theta()*TMath::RadToDeg()<45
    && el_Scattered->Theta()*TMath::RadToDeg()>5 && el_Scattered->Theta()*TMath::RadToDeg()<35){

      h_eFD_KFD_NK0->Fill(Inv_Lambda.M(), TotalWeight);
      h_Numerator->Fill(6,TotalWeight);


      if(Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Lambda_Proton->Theta()*TMath::RadToDeg()<45
      && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
        h_Numerator->Fill(7,TotalWeight);
        h_LProton_LPiMinus_NK0->Fill(Inv_Lambda.M(), TotalWeight);

        if(Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()<45
        && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
          h_Numerator->Fill(9,TotalWeight);
          h_LProton_LPiMinus_CLProton_CLPiMinus_CPiMinus_NK0->Fill(Inv_Lambda.M(), TotalWeight);

          if(Cascade_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()<45){

            h_Numerator->Fill(10,TotalWeight);
            h_LProton_LPiMinus_CPiMinus_NK0->Fill(Inv_Lambda.M(), TotalWeight);
          }
        }

      }
      if(Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()<45
      && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45
      && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()<45){
        h_Numerator->Fill(8,TotalWeight);
        h_CPiMinus_CLProton_CLPiMinus_NK0->Fill(Inv_Lambda.M(), TotalWeight);
      }


      if(Kaon0_Pion->Theta()*TMath::RadToDeg()>5 && Kaon0_Pion->Theta()*TMath::RadToDeg()<45
      && Kaon0_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Kaon0_Pi_Minus->Theta()*TMath::RadToDeg()<45){

        h_eFD_PiMinusFD_PiPlusFD->Fill(Inv_Lambda.M(),TotalWeight);
        h_Numerator->Fill(1,TotalWeight);


        if(Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Lambda_Proton->Theta()*TMath::RadToDeg()<45
        && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
          h_Numerator->Fill(2,TotalWeight);
          h_LProton_LPiMinus->Fill(Inv_Lambda.M(),TotalWeight);


          if(Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()<45
          && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
            h_LProton_LPiMinus_CLProton_CLPiMinus->Fill(Inv_Lambda.M(), TotalWeight);
            h_Numerator->Fill(4,TotalWeight);

            if(Cascade_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()<45){
              h_LProton_LPiMinus_CLProton_CLPiMinus_CPiMinus->Fill(Inv_Lambda.M(), TotalWeight);
              h_Numerator->Fill(5,TotalWeight);
            }
          }



        }
        if(Cascade_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()<45
        && Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Proton->Theta()*TMath::RadToDeg()<45
        && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45){
          h_CLProton_CLPiMinus_CPiMinus->Fill(Inv_Lambda.M(),TotalWeight);
          h_Numerator->Fill(3,TotalWeight);
        }
      }

      // case 3 cascade pion, cascade lambda proton



    }
    h_AcceptanceAll->Divide(h_Numerator, h_Denominator);

    i++;
  }



  //if(Kaon_Plus1->Theta()*TMath::RadToDeg()>5 && Kaon_Plus1->Theta()*TMath::RadToDeg()<45 && Kaon_Plus2->Theta()*TMath::RadToDeg()>5 && Kaon_Plus2->Theta()*TMath::RadToDeg()<45
  //&& Kaon_0->Theta()*TMath::RadToDeg()>5 && Kaon_0->Theta()*TMath::RadToDeg()<45 && Lambda_Proton->Theta()*TMath::RadToDeg()>5 && Lambda_Proton->Theta()*TMath::RadToDeg()<45
  //&& Lambda_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Lambda_Pi_Minus->Theta()*TMath::RadToDeg()<45 && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()>5 && Cascade_Pi_Minus->Theta()*TMath::RadToDeg()<45){
  //h_LProton_LPiMinus_CPiMinus->Fill(Inv_Lambda.M(),TotalWeight);
  //h_Numerator->Fill(3,TotalWeight);






  // if K0 and K+ K+ in FD, electron, pion- , pion+  case 1
  // if K0 and K+ K+ in FD, p , pi - from lambda initial case 2
  // if K0 and K+ K+ in FD, p , pi-, pi - from cascade decay and lambda decay case 3
  // if p, pi - from lambda initial from case 2, p pi - from cascade_lambda case 4
  // case 4 with pi -
  // everything except K0

  // Find acceptance by dividing histogram after angle cut by original histogram
  // -> Divide (A,B) gives A/B
  // acceptance with and without pweight **


  //print:(Theta_Kaon);
  fileOutput1->Write();
}
