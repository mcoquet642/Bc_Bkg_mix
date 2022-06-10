#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
using namespace std;

bool IsinAcc(double pt, double eta, double ptcut){
	if (pt > ptcut && eta > 2.5 && eta < 3.6)
	{
		return true;
	}else{
		return false;
		}
}

bool PrevDecay(double decayZ){
	double MFT_f0Z = 460*1000;
	return (decayZ < MFT_f0Z);
}

void muon_pairs_bkg(int nb_events, char * output, char * inputPip, char * inputPim, char * inputDp, char * inputDm, char * inputKp, char * inputKm, double m_lim_inf, double m_lim_sup, double pt_cut, char * inputJPsi, double proba_Jpsi){


   TFile *fjpsi = new TFile(inputJPsi);
   TTree *tjpsi = (TTree*)fjpsi->Get("DecayTree");

   Double_t Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_pt, Jpsi_y, Jpsi_sigip, Jpsi_ip, Jpsi_vtxX, Jpsi_vtxY, Jpsi_vtxZ, Jpsi_m, Jpsi_e;
   Double_t mum_jpsi_vtxX, mum_jpsi_vtxY, mum_jpsi_vtxZ, mum_jpsi_PX, mum_jpsi_PY, mum_jpsi_PZ, mum_jpsi_IP;
   Double_t mup_jpsi_vtxX, mup_jpsi_vtxY, mup_jpsi_vtxZ, mup_jpsi_PX, mup_jpsi_PY, mup_jpsi_PZ, mup_jpsi_IP;

   tjpsi->SetBranchAddress("Jpsi_PT",&Jpsi_pt);
   tjpsi->SetBranchAddress("Jpsi_PX",&Jpsi_px);
   tjpsi->SetBranchAddress("Jpsi_PY",&Jpsi_py);
   tjpsi->SetBranchAddress("Jpsi_PZ",&Jpsi_pz);
   tjpsi->SetBranchAddress("Jpsi_M",&Jpsi_m);
   tjpsi->SetBranchAddress("Jpsi_y",&Jpsi_y);
   tjpsi->SetBranchAddress("Jpsi_origX",&Jpsi_vtxX);
   tjpsi->SetBranchAddress("Jpsi_origY",&Jpsi_vtxY);
   tjpsi->SetBranchAddress("Jpsi_origZ",&Jpsi_vtxZ);
   tjpsi->SetBranchAddress("Jpsi_E",&Jpsi_e);
   tjpsi->SetBranchAddress("Jpsi_IP",&Jpsi_ip);
   tjpsi->SetBranchAddress("Jpsi_SIGMAIP_TRUE",&Jpsi_sigip);
   tjpsi->SetBranchAddress("mum_jpsi_origX",&mum_jpsi_vtxX);
   tjpsi->SetBranchAddress("mum_jpsi_origY",&mum_jpsi_vtxY);
   tjpsi->SetBranchAddress("mum_jpsi_origZ",&mum_jpsi_vtxZ);
   tjpsi->SetBranchAddress("mum_jpsi_PX",&mum_jpsi_PX);
   tjpsi->SetBranchAddress("mum_jpsi_PY",&mum_jpsi_PY);
   tjpsi->SetBranchAddress("mum_jpsi_PZ",&mum_jpsi_PZ);
   tjpsi->SetBranchAddress("mup_jpsi_origX",&mup_jpsi_vtxX);
   tjpsi->SetBranchAddress("mup_jpsi_origY",&mup_jpsi_vtxY);
   tjpsi->SetBranchAddress("mup_jpsi_origZ",&mup_jpsi_vtxZ);
   tjpsi->SetBranchAddress("mup_jpsi_PX",&mup_jpsi_PX);
   tjpsi->SetBranchAddress("mup_jpsi_PY",&mup_jpsi_PY);
   tjpsi->SetBranchAddress("mup_jpsi_PZ",&mup_jpsi_PZ);
   tjpsi->SetBranchAddress("mum_jpsi_IP",&mum_jpsi_IP);
   tjpsi->SetBranchAddress("mup_jpsi_IP",&mup_jpsi_IP);

	Int_t entries_Jpsi = (Int_t)tjpsi->GetEntries();

   TFile *f = new TFile(inputPip);
   TTree *t1 = (TTree*)f->Get("DecayTree");
   TFile *f_2 = new TFile(inputPim);
   TTree *t1_2 = (TTree*)f_2->Get("DecayTree");
   Double_t px_m, py_m, pz_m, e_m, ip_m, sig_m, origX_m, origY_m, origZ_m, Pim, pt_m, eta_m;
   Double_t px_p, py_p, pz_p, e_p, ip_p, sig_p, origX_p, origY_p, origZ_p, Pip, pt_p, eta_p;

   TBranch *br1_pxm = t1_2->GetBranch("mum_PX");
   TBranch *br1_pxp = t1->GetBranch("mup_PX");
   TBranch *br1_pym = t1_2->GetBranch("mum_PY");
   TBranch *br1_pyp = t1->GetBranch("mup_PY");
   TBranch *br1_pzm = t1_2->GetBranch("mum_PZ");
   TBranch *br1_pzp = t1->GetBranch("mup_PZ");
   TBranch *br1_em = t1_2->GetBranch("mum_E");
   TBranch *br1_ep = t1->GetBranch("mup_E");
   TBranch *br1_ipm = t1_2->GetBranch("mum_IP");
   TBranch *br1_ipp = t1->GetBranch("mup_IP");
   TBranch *br1_sigm = t1_2->GetBranch("mum_SIGMAIP_TRUE");
   TBranch *br1_sigp = t1->GetBranch("mup_SIGMAIP_TRUE");
   TBranch *br1_origXm = t1_2->GetBranch("mum_origX");
   TBranch *br1_origXp = t1->GetBranch("mup_origX");
   TBranch *br1_origYm = t1_2->GetBranch("mum_origY");
   TBranch *br1_origYp = t1->GetBranch("mup_origY");
   TBranch *br1_origZm = t1_2->GetBranch("mum_origZ");
   TBranch *br1_origZp = t1->GetBranch("mup_origZ");
   TBranch *br1_Pim = t1_2->GetBranch("pim_origZ");
   TBranch *br1_Pip = t1->GetBranch("pip_origZ");
   TBranch *br1_ptm = t1_2->GetBranch("mum_PT");
   TBranch *br1_ptp = t1->GetBranch("mup_PT");
   TBranch *br1_etam = t1_2->GetBranch("mum_eta");
   TBranch *br1_etap = t1->GetBranch("mup_eta");



   br1_pxm->SetAddress(&px_m);
   br1_pxp->SetAddress(&px_p);
   br1_pym->SetAddress(&py_m);
   br1_pyp->SetAddress(&py_p);
   br1_pzm->SetAddress(&pz_m);
   br1_pzp->SetAddress(&pz_p);
   br1_em->SetAddress(&e_m);
   br1_ep->SetAddress(&e_p);
   br1_ipm->SetAddress(&ip_m);
   br1_ipp->SetAddress(&ip_p);
   br1_sigm->SetAddress(&sig_m);
   br1_sigp->SetAddress(&sig_p);
   br1_origXm->SetAddress(&origX_m);
   br1_origXp->SetAddress(&origX_p);
   br1_origYm->SetAddress(&origY_m);
   br1_origYp->SetAddress(&origY_p);
   br1_origZm->SetAddress(&origZ_m);
   br1_origZp->SetAddress(&origZ_p);
   br1_Pim->SetAddress(&Pim);
   br1_Pip->SetAddress(&Pip);
   br1_ptm->SetAddress(&pt_m);
   br1_ptp->SetAddress(&pt_p);
   br1_etam->SetAddress(&eta_m);
   br1_etap->SetAddress(&eta_p);

	Int_t entries_m_Pip = (Int_t)br1_em->GetEntries();
	Int_t entries_p_Pip = (Int_t)br1_ep->GetEntries();

   TFile *f2 = new TFile(inputDp);
   TTree *t_Dp = (TTree*)f2->Get("DecayTree");
   TFile *f2_2 = new TFile(inputDm);
   TTree *t_Dm = (TTree*)f2_2->Get("DecayTree");
   Double_t px_m_2, py_m_2, pz_m_2, e_m_2, ip_m_2, sig_m_2, origX_m_2, origY_m_2, origZ_m_2, Dm, pt_m_2, eta_m_2;
   Double_t px_p_2, py_p_2, pz_p_2, e_p_2, ip_p_2, sig_p_2, origX_p_2, origY_p_2, origZ_p_2, Dp, pt_p_2, eta_p_2;


   TBranch *br2_pxm = t_Dm->GetBranch("mum_PX");
   TBranch *br2_pxp = t_Dp->GetBranch("mup_PX");
   TBranch *br2_pym = t_Dm->GetBranch("mum_PY");
   TBranch *br2_pyp = t_Dp->GetBranch("mup_PY");
   TBranch *br2_pzm = t_Dm->GetBranch("mum_PZ");
   TBranch *br2_pzp = t_Dp->GetBranch("mup_PZ");
   TBranch *br2_em = t_Dm->GetBranch("mum_E");
   TBranch *br2_ep = t_Dp->GetBranch("mup_E");
   TBranch *br2_ipm = t_Dm->GetBranch("mum_IP");
   TBranch *br2_ipp = t_Dp->GetBranch("mup_IP");
   TBranch *br2_sigm = t_Dm->GetBranch("mum_SIGMAIP_TRUE");
   TBranch *br2_sigp = t_Dp->GetBranch("mup_SIGMAIP_TRUE");
   TBranch *br2_origXm = t_Dm->GetBranch("mum_origX");
   TBranch *br2_origXp = t_Dp->GetBranch("mup_origX");
   TBranch *br2_origYm = t_Dm->GetBranch("mum_origY");
   TBranch *br2_origYp = t_Dp->GetBranch("mup_origY");
   TBranch *br2_origZm = t_Dm->GetBranch("mum_origZ");
   TBranch *br2_origZp = t_Dp->GetBranch("mup_origZ");
   TBranch *br2_Dm = t_Dm->GetBranch("Dm_origZ");
   TBranch *br2_Dp = t_Dp->GetBranch("Dp_origZ");
   TBranch *br2_ptm = t_Dm->GetBranch("mum_PT");
   TBranch *br2_ptp = t_Dp->GetBranch("mup_PT");
   TBranch *br2_etam = t_Dm->GetBranch("mum_eta");
   TBranch *br2_etap = t_Dp->GetBranch("mup_eta");

   br2_pxm->SetAddress(&px_m_2);
   br2_pxp->SetAddress(&px_p_2);
   br2_pym->SetAddress(&py_m_2);
   br2_pyp->SetAddress(&py_p_2);
   br2_pzm->SetAddress(&pz_m_2);
   br2_pzp->SetAddress(&pz_p_2);
   br2_em->SetAddress(&e_m_2);
   br2_ep->SetAddress(&e_p_2);
   br2_ipm->SetAddress(&ip_m_2);
   br2_ipp->SetAddress(&ip_p_2);
   br2_sigm->SetAddress(&sig_m_2);
   br2_sigp->SetAddress(&sig_p_2);
   br2_origXm->SetAddress(&origX_m_2);
   br2_origXp->SetAddress(&origX_p_2);
   br2_origYm->SetAddress(&origY_m_2);
   br2_origYp->SetAddress(&origY_p_2);
   br2_origZm->SetAddress(&origZ_m_2);
   br2_origZp->SetAddress(&origZ_p_2);
   br2_Dm->SetAddress(&Dm);
   br2_Dp->SetAddress(&Dp);
   br2_ptm->SetAddress(&pt_m_2);
   br2_ptp->SetAddress(&pt_p_2);
   br2_etam->SetAddress(&eta_m_2);
   br2_etap->SetAddress(&eta_p_2);


	Int_t entries_m_Dp = (Int_t)br2_em->GetEntries();
	Int_t entries_p_Dp = (Int_t)br2_ep->GetEntries();


   TFile *f3 = new TFile(inputKp);
   TTree *t_K = (TTree*)f3->Get("DecayTree");
   TFile *f3_2 = new TFile(inputKm);
   TTree *t_K_2 = (TTree*)f3_2->Get("DecayTree");
   Double_t px_m_3, py_m_3, pz_m_3, e_m_3, ip_m_3, sig_m_3, origX_m_3, origY_m_3, origZ_m_3, Km, pt_m_3, eta_m_3;
   Double_t px_p_3, py_p_3, pz_p_3, e_p_3, ip_p_3, sig_p_3, origX_p_3, origY_p_3, origZ_p_3, Kp, pt_p_3, eta_p_3;


   TBranch *br3_pxm = t_K_2->GetBranch("mum_PX");
   TBranch *br3_pxp = t_K->GetBranch("mup_PX");
   TBranch *br3_pym = t_K_2->GetBranch("mum_PY");
   TBranch *br3_pyp = t_K->GetBranch("mup_PY");
   TBranch *br3_pzm = t_K_2->GetBranch("mum_PZ");
   TBranch *br3_pzp = t_K->GetBranch("mup_PZ");
   TBranch *br3_em = t_K_2->GetBranch("mum_E");
   TBranch *br3_ep = t_K->GetBranch("mup_E");
   TBranch *br3_ipm = t_K_2->GetBranch("mum_IP");
   TBranch *br3_ipp = t_K->GetBranch("mup_IP");
   TBranch *br3_sigm = t_K_2->GetBranch("mum_SIGMAIP_TRUE");
   TBranch *br3_sigp = t_K->GetBranch("mup_SIGMAIP_TRUE");
   TBranch *br3_origXm = t_K_2->GetBranch("mum_origX");
   TBranch *br3_origXp = t_K->GetBranch("mup_origX");
   TBranch *br3_origYm = t_K_2->GetBranch("mum_origY");
   TBranch *br3_origYp = t_K->GetBranch("mup_origY");
   TBranch *br3_origZm = t_K_2->GetBranch("mum_origZ");
   TBranch *br3_origZp = t_K->GetBranch("mup_origZ");
   TBranch *br3_Km = t_K_2->GetBranch("Km_origZ");
   TBranch *br3_Kp = t_K->GetBranch("Kp_origZ");
   TBranch *br3_ptm = t_K_2->GetBranch("mum_PT");
   TBranch *br3_ptp = t_K->GetBranch("mup_PT");
   TBranch *br3_etam = t_K_2->GetBranch("mum_eta");
   TBranch *br3_etap = t_K->GetBranch("mup_eta");

   br3_pxm->SetAddress(&px_m_3);
   br3_pxp->SetAddress(&px_p_3);
   br3_pym->SetAddress(&py_m_3);
   br3_pyp->SetAddress(&py_p_3);
   br3_pzm->SetAddress(&pz_m_3);
   br3_pzp->SetAddress(&pz_p_3);
   br3_em->SetAddress(&e_m_3);
   br3_ep->SetAddress(&e_p_3);
   br3_ipm->SetAddress(&ip_m_3);
   br3_ipp->SetAddress(&ip_p_3);
   br3_sigm->SetAddress(&sig_m_3);
   br3_sigp->SetAddress(&sig_p_3);
   br3_origXm->SetAddress(&origX_m_3);
   br3_origXp->SetAddress(&origX_p_3);
   br3_origYm->SetAddress(&origY_m_3);
   br3_origYp->SetAddress(&origY_p_3);
   br3_origZm->SetAddress(&origZ_m_3);
   br3_origZp->SetAddress(&origZ_p_3);
   br3_Km->SetAddress(&Km);
   br3_Kp->SetAddress(&Kp);
   br3_ptm->SetAddress(&pt_m_3);
   br3_ptp->SetAddress(&pt_p_3);
   br3_etam->SetAddress(&eta_m_3);
   br3_etap->SetAddress(&eta_p_3);


	Int_t entries_m_K = (Int_t)br3_em->GetEntries();
	Int_t entries_p_K = (Int_t)br3_ep->GetEntries();


   Int_t entries_min_Pip = min(entries_m_Pip, entries_p_Pip); 
   cout << "Total number of Pi muons provided : " << entries_min_Pip << endl;


   Int_t entries_min_Dp = min(entries_m_Dp, entries_p_Dp); 
   cout << "Total number of Dp muons provided : " << entries_min_Dp << endl;

   Int_t entries_min_K = min(entries_m_K, entries_p_K); 
   cout << "Total number of K muons provided : " << entries_min_K << endl;

   cout << "Total number of Jpsi provided : " << entries_Jpsi << endl;

//Decay Probabilities_______________


//   Double_t proba_Jpsi=0.01404; // Nb bkg j/psi /event for (2.5 < eta <3.6) and pT(single mu)>1 GeV/c


   Double_t f_Dp = 0.232;

   Double_t BR_Pi = 0.9999;//PDG, inclusive ...
   Double_t BR_K = 0.6356;
   Double_t BR_Dp = 0.16;

//   Double_t sigma_pp = 3.744; //FONLL ... (mb) -10<y<10, pT>0
//   Double_t sigma_pp = 1.438; //FONLL ... (mb) 1<eta<6, pT>0
   Double_t sigma_pp = 0.385; //FONLL ... (mb) 2.5<eta<3.6, pT>0
//   Double_t sigma_pp = 0.3722; //FONLL ... (mb) 2.5<eta<3.6, pT>0.5
//   Double_t sigma_pp = 0.329; //FONLL ... (mb) 2.5<eta<3.6, pT>1
   Double_t t_AA = 23.26;//(mb^-1, 0-10% centrality)
   
   Double_t Npi=852;// 2.5<eta<3.6
   Double_t Nk=127;//2.5<eta<3.6

   Double_t Proba_pi=0.003239;//P(origZ<46cm)
   Double_t Proba_k=0.010590;//P(origZ<46cm) 

//   Double_t Npi=3072;// 1<eta<5
//   Double_t Nk=456;// 1<eta<5

   Int_t nb_charm=floor(sigma_pp*t_AA);
   cout << "Mean charm number /event : " << nb_charm << endl;

//____________________________

//Partie combinatoireJpsi_y

   TRandom3 r1, r2, r3, r4, r5, r6;

   cout << "Number of events : " << nb_events << endl;


//_________





   Double_t ev, dep2;
   Double_t jpsi_y, jpsi_pt, jpsi_ip, jpsi_sigip, jpsi_vtxX, jpsi_vtxY, jpsi_vtxZ, jpsi_m, jpsi_pseudo_tau, jpsi_pz;
   Double_t trimu_m, trimu_pseudo_tau, trimu_ip, trimu_sigip, trimu_vtxX, trimu_vtxY, trimu_vtxZ, trimu_pz, trimu_pt;
   Double_t costhetap, disp, diff_pseudo_tau;
   Int_t charge, signal;
   Double_t mu_IP, mu_sigIP;
   Double_t mum_jpsi_ip;
   Double_t mup_jpsi_ip;
   gROOT->cd();
   TFile * muon_out = new TFile(output, "recreate");
   TTree t2("BcDecay","dimuon Tree");

   t2.Branch("jpsi_ip",&jpsi_ip,"jpsi_ip/D");
   t2.Branch("trimu_ip",&trimu_ip,"trimu_ip/D");
   t2.Branch("jpsi_pseudo_tau",&jpsi_pseudo_tau,"jpsi_pseudo_tau/D");
   t2.Branch("cos_theta_p",&costhetap,"cos_theta_p/D");
   t2.Branch("disp",&disp,"disp/D");
   t2.Branch("diff_pseudo_tau",&diff_pseudo_tau,"diff_pseudo_tau/D");
   t2.Branch("trimu_pseudo_tau",&trimu_pseudo_tau,"trimu_pseudo_tau/D");

   t2.Branch("mum_jpsi_ip",&mum_jpsi_ip,"mum_jpsi_ip/D");
   t2.Branch("mup_jpsi_ip",&mup_jpsi_ip,"mup_jpsi_ip/D");

   t2.Branch("mu_IP",&mu_IP,"mu_IP/D");

   t2.Branch("jpsi_m",&jpsi_m,"jpsi_m/D");
   t2.Branch("trimu_m",&trimu_m,"trimu_m/D");

   t2.Branch("mu_sigIP",&mu_sigIP,"mu_sigIP/D");
   t2.Branch("jpsi_sigip",&jpsi_sigip,"jpsi_sigip/D");
   t2.Branch("jpsi_y",&jpsi_y,"jpsi_y/D");
   t2.Branch("jpsi_pt",&jpsi_pt,"jpsi_pt/D");
   t2.Branch("jpsi_PZ",&jpsi_pz,"jpsi_PZ/D");
   t2.Branch("jpsi_vtxX",&jpsi_vtxX,"jpsi_vtxX/D");
   t2.Branch("jpsi_vtxY",&jpsi_vtxY,"jpsi_vtxY/D");
   t2.Branch("jpsi_vtxZ",&jpsi_vtxZ,"jpsi_vtxZ/D");
   t2.Branch("trimu_vtxX",&trimu_vtxX,"trimu_vtxX/D");
   t2.Branch("trimu_vtxY",&trimu_vtxY,"trimu_vtxY/D");
   t2.Branch("trimu_vtxZ",&trimu_vtxZ,"trimu_vtxZ/D");
   t2.Branch("trimu_PZ",&trimu_pz,"trimu_PZ/D");
   t2.Branch("trimu_pt",&trimu_pt,"trimu_pt/D");
   t2.Branch("trimu_sigip",&trimu_sigip,"trimu_sigip/D");
   t2.Branch("events",&ev,"ev/D");
   t2.Branch("charge",&charge,"charge/I");
   t2.Branch("signal",&signal,"signal/I");

   TLorentzVector p1;
   TLorentzVector p2;
   TLorentzVector ptot;

   Int_t nb_muon_Pip=0;
   Int_t nb_muon_Dp=0;
   Int_t nb_muon_K=0;
   Int_t nb_muon_Pim=0;
   Int_t nb_muon_Dm=0;
   Int_t nb_muon_Kb=0;

   Int_t nb_Jpsi=0;

   Int_t somme = 0;

   Int_t list_hadrons[3*nb_events];
   Int_t blist_hadrons[3*nb_events];
   Int_t buff=0;

   Int_t list_jpsi[nb_events];

   for (Int_t k=0; k<nb_events; k++){

	list_jpsi[k]=r5.Binomial(1, proba_Jpsi);
 
	list_hadrons[3*k]=r1.Binomial(Npi, Proba_pi*BR_Pi);
	list_hadrons[3*k+1]=r2.Binomial(nb_charm, f_Dp*BR_Dp);
	list_hadrons[3*k+2]=r3.Binomial(Nk, Proba_k*BR_K);

	blist_hadrons[3*k]=r1.Binomial(Npi, Proba_pi*BR_Pi);
	blist_hadrons[3*k+1]=r2.Binomial(nb_charm, f_Dp*BR_Dp);
	blist_hadrons[3*k+2]=r3.Binomial(Nk, Proba_k*BR_K);

 	nb_muon_Pip+=list_hadrons[3*k];
 	nb_muon_Dp+=list_hadrons[3*k+1];
 	nb_muon_K+=list_hadrons[3*k+2];

 	nb_muon_Pim+=blist_hadrons[3*k];
 	nb_muon_Dm+=blist_hadrons[3*k+1];
 	nb_muon_Kb+=blist_hadrons[3*k+2];

	nb_Jpsi+=list_jpsi[k];


	buff = (list_hadrons[3*k]+list_hadrons[3*k+1]+list_hadrons[3*k+2]+blist_hadrons[3*k]+blist_hadrons[3*k+1]+blist_hadrons[3*k+2])*list_jpsi[k];
	somme += buff;
  	if (k%10000 == 0){
		cout << " Event : " << k << " , Number of Pip : " << list_hadrons[3*k] << " , Number of Dp : " << list_hadrons[3*k+1] <<  " , Number of K : " << list_hadrons[3*k+2] << " , Number of Pim : " << blist_hadrons[3*k] << " , Number of Dm : " << blist_hadrons[3*k+1] <<  " , Number of Kb : " << blist_hadrons[3*k+2] <<  " , Number of Jpsi : " << list_jpsi[k] << " , Number of trimuons : " << buff << endl;
  	}
   }

   cout << "Total number of Pip muons : " << nb_muon_Pip << ", Total number of Dp muons : " << nb_muon_Dp << ", Total number of K muons : " << nb_muon_K << "Total number of Pim muons : " << nb_muon_Pim << ", Total number of Dm muons : " << nb_muon_Dm << ", Total number of Kb muons : " << nb_muon_Kb <<  ", Total number of Jpsi : " << nb_Jpsi << ", Total number of trimuons : " << somme << endl;

   if (nb_muon_Pip > entries_min_Pip || nb_muon_Dp > entries_min_Dp || nb_muon_K > entries_min_K || nb_Jpsi > entries_Jpsi){
	cout << "Note enough muons provided." << endl;
   }

   Int_t count = 0;
   Int_t index_ev = 0;
   Int_t index_Pip = 0;
   Int_t index_Dp = 0;
   Int_t index_K = 0;

   Int_t index_Pim = 0;
   Int_t index_Dm = 0;
   Int_t index_Kb = 0;

   Int_t nb_final=0;
   Int_t nb_init=0;

   signal = 0;

   for (Int_t i=0; i<nb_Jpsi; i++){
	while (list_jpsi[index_ev]==0){
		index_Pim+=list_hadrons[3*index_ev];
		index_Dm+=list_hadrons[3*index_ev+1];
		index_Kb+=list_hadrons[3*index_ev+2];

		index_Pip+=blist_hadrons[3*index_ev];
		index_Dp+=blist_hadrons[3*index_ev+1];
		index_K+=blist_hadrons[3*index_ev+2];
		index_ev++;
	}
	tjpsi->GetEntry(i);

	p1.SetPxPyPzE(Jpsi_px,Jpsi_py,Jpsi_pz,Jpsi_e);
	double a = p1.P()*p1.P();
	jpsi_y=Jpsi_y;
	jpsi_pt=Jpsi_pt;
	jpsi_vtxX=Jpsi_vtxX;
	jpsi_vtxY=Jpsi_vtxY;
	jpsi_vtxZ=Jpsi_vtxZ;
	jpsi_m=Jpsi_m;
	jpsi_pz=Jpsi_pz;
	jpsi_sigip=1.2533*(30 + 60/jpsi_pt);
        double sigma1 = gRandom->Gaus(0.,jpsi_sigip);
	jpsi_ip=Jpsi_ip+sigma1;
	jpsi_pseudo_tau=Jpsi_vtxZ*Jpsi_m/Jpsi_pz;

	mum_jpsi_ip=mum_jpsi_IP;
	mup_jpsi_ip=mup_jpsi_IP;

   	for (Int_t j=0; j<blist_hadrons[3*index_ev]; j++) {
		t1->GetEntry(j+index_Pip);
		ev=index_ev;
		dep2=abs(Pip - origZ_p);

		p2.SetPxPyPzE(px_p,py_p,pz_p,e_p);
		ptot=p1+p2;
		trimu_m=ptot.M();
		trimu_pz=ptot.Pz();
		
		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(jpsi_vtxX - origX_p)+p1.Py()*(jpsi_vtxY - origY_p)+p1.Pz()*(jpsi_vtxZ - origZ_p);
		double e = p2.Px()*(jpsi_vtxX - origX_p)+p2.Py()*(jpsi_vtxY - origY_p)+p2.Pz()*(jpsi_vtxZ - origZ_p);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		trimu_vtxX = (jpsi_vtxX + origX_p + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxY = (jpsi_vtxY + origY_p + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxZ = (jpsi_vtxZ + origZ_p + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip=trimu_vtxX-ptot.Px()*trimu_vtxZ/ptot.Pz();
		double yip=trimu_vtxY-ptot.Py()*trimu_vtxZ/ptot.Pz();
		trimu_ip=sqrt(xip*xip+yip*yip);
		trimu_sigip=1.2533*(30 + 60/ptot.Pt());

		trimu_pseudo_tau=trimu_vtxZ*trimu_m/ptot.Pz();
		trimu_pt = ptot.Pt();

		double r2 = trimu_vtxX*trimu_vtxX +trimu_vtxY*trimu_vtxY+trimu_vtxZ*trimu_vtxZ;
		double r = sqrt(r2);
		costhetap = (trimu_vtxX*ptot.Px()+trimu_vtxY*ptot.Py()+trimu_vtxZ*ptot.Pz())/(r*ptot.P());
		
		double xdca1 = mum_jpsi_vtxX+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PX/mum_jpsi_PZ;
		double ydca1 = mum_jpsi_vtxY+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PY/mum_jpsi_PZ;

		double xdca2 = mup_jpsi_vtxX+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PX/mup_jpsi_PZ;
		double ydca2 = mup_jpsi_vtxY+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PY/mup_jpsi_PZ;;

		double xdca3 = origX_p + (trimu_vtxZ - origZ_p)*px_p/pz_p;
		double ydca3 = origY_p + (trimu_vtxZ - origZ_p)*py_p/pz_p;

		disp = sqrt(xdca1*xdca1 + ydca1*ydca1 + xdca2*xdca2 + ydca2*ydca2 + xdca3*xdca3 + ydca3*ydca3);
		diff_pseudo_tau = abs(jpsi_pseudo_tau-trimu_pseudo_tau);

		charge=1;
		mu_IP=ip_p;
		mu_sigIP=sig_p;

		if (trimu_m > m_lim_inf && trimu_m < m_lim_sup && IsinAcc(pt_p, eta_p, pt_cut) && PrevDecay(origZ_p)){
			t2.Fill();
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+1]; j++) {
		t_Dp->GetEntry(j+index_Dp);
		ev=index_ev;
		dep2=abs(Dp - origZ_p_2);

		p2.SetPxPyPzE(px_p_2,py_p_2,pz_p_2,e_p_2);
		ptot=p1+p2;
		trimu_m=ptot.M();
		trimu_pz=ptot.Pz();
		
		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(jpsi_vtxX - origX_p_2)+p1.Py()*(jpsi_vtxY - origY_p_2)+p1.Pz()*(jpsi_vtxZ - origZ_p_2);
		double e = p2.Px()*(jpsi_vtxX - origX_p_2)+p2.Py()*(jpsi_vtxY - origY_p_2)+p2.Pz()*(jpsi_vtxZ - origZ_p_2);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		trimu_vtxX = (jpsi_vtxX + origX_p_2 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxY = (jpsi_vtxY + origY_p_2 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxZ = (jpsi_vtxZ + origZ_p_2 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip=trimu_vtxX-ptot.Px()*trimu_vtxZ/ptot.Pz();
		double yip=trimu_vtxY-ptot.Py()*trimu_vtxZ/ptot.Pz();
		trimu_ip=sqrt(xip*xip+yip*yip);
		trimu_sigip=1.2533*(30 + 60/ptot.Pt());

		trimu_pseudo_tau=trimu_vtxZ*trimu_m/ptot.Pz();
		trimu_pt = ptot.Pt();

		double r2 = trimu_vtxX*trimu_vtxX +trimu_vtxY*trimu_vtxY+trimu_vtxZ*trimu_vtxZ;
		double r = sqrt(r2);
		costhetap = (trimu_vtxX*ptot.Px()+trimu_vtxY*ptot.Py()+trimu_vtxZ*ptot.Pz())/(r*ptot.P());

		double xdca1 = mum_jpsi_vtxX+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PX/mum_jpsi_PZ;
		double ydca1 = mum_jpsi_vtxY+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PY/mum_jpsi_PZ;

		double xdca2 = mup_jpsi_vtxX+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PX/mup_jpsi_PZ;
		double ydca2 = mup_jpsi_vtxY+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PY/mup_jpsi_PZ;;

		double xdca3 = origX_p_2 + (trimu_vtxZ - origZ_p_2)*px_p_2/pz_p_2;
		double ydca3 = origY_p_2 + (trimu_vtxZ - origZ_p_2)*py_p_2/pz_p_2;

		disp = sqrt(xdca1*xdca1 + ydca1*ydca1 + xdca2*xdca2 + ydca2*ydca2 + xdca3*xdca3 + ydca3*ydca3);
		diff_pseudo_tau = abs(jpsi_pseudo_tau-trimu_pseudo_tau);

		charge=1;
		mu_IP=ip_p_2;
		mu_sigIP=sig_p_2;

		if (trimu_m > m_lim_inf && trimu_m < m_lim_sup && IsinAcc(pt_p_2, eta_p_2, pt_cut) && PrevDecay(origZ_p_2)){
			t2.Fill();
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+2]; j++) {
		t_K->GetEntry(j+index_K);
		ev=index_ev;
		dep2=abs(Kp - origZ_p_3);

		p2.SetPxPyPzE(px_p_3,py_p_3,pz_p_3,e_p_3);
		ptot=p1+p2;
		trimu_m=ptot.M();
		trimu_pz=ptot.Pz();
		
		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(jpsi_vtxX - origX_p_3)+p1.Py()*(jpsi_vtxY - origY_p_3)+p1.Pz()*(jpsi_vtxZ - origZ_p_3);
		double e = p2.Px()*(jpsi_vtxX - origX_p_3)+p2.Py()*(jpsi_vtxY - origY_p_3)+p2.Pz()*(jpsi_vtxZ - origZ_p_3);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		trimu_vtxX = (jpsi_vtxX + origX_p_3 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxY = (jpsi_vtxY + origY_p_3 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxZ = (jpsi_vtxZ + origZ_p_3 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip=trimu_vtxX-ptot.Px()*trimu_vtxZ/ptot.Pz();
		double yip=trimu_vtxY-ptot.Py()*trimu_vtxZ/ptot.Pz();
		trimu_ip=sqrt(xip*xip+yip*yip);
		trimu_sigip=1.2533*(30 + 60/ptot.Pt());

		trimu_pseudo_tau=trimu_vtxZ*trimu_m/ptot.Pz();
		trimu_pt = ptot.Pt();

		double r2 = trimu_vtxX*trimu_vtxX +trimu_vtxY*trimu_vtxY+trimu_vtxZ*trimu_vtxZ;
		double r = sqrt(r2);
		costhetap = (trimu_vtxX*ptot.Px()+trimu_vtxY*ptot.Py()+trimu_vtxZ*ptot.Pz())/(r*ptot.P());

		double xdca1 = mum_jpsi_vtxX+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PX/mum_jpsi_PZ;
		double ydca1 = mum_jpsi_vtxY+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PY/mum_jpsi_PZ;

		double xdca2 = mup_jpsi_vtxX+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PX/mup_jpsi_PZ;
		double ydca2 = mup_jpsi_vtxY+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PY/mup_jpsi_PZ;;

		double xdca3 = origX_p_3 + (trimu_vtxZ - origZ_p_3)*px_p_3/pz_p_3;
		double ydca3 = origY_p_3 + (trimu_vtxZ - origZ_p_3)*py_p_3/pz_p_3;

		disp = sqrt(xdca1*xdca1 + ydca1*ydca1 + xdca2*xdca2 + ydca2*ydca2 + xdca3*xdca3 + ydca3*ydca3);
		diff_pseudo_tau = abs(jpsi_pseudo_tau-trimu_pseudo_tau);

		charge=1;
		mu_IP=ip_p_3;
		mu_sigIP=sig_p_3;

		if (trimu_m > m_lim_inf && trimu_m < m_lim_sup && IsinAcc(pt_p_3, eta_p_3, pt_cut) && PrevDecay(origZ_p_3)){
			t2.Fill();
		}
	}

	//------------------------------------------------------------------------------
	
   	for (Int_t j=0; j<list_hadrons[3*index_ev]; j++) {
		t1_2->GetEntry(j+index_Pim);
		ev=index_ev;
		dep2=abs(Pim - origZ_m);

		p2.SetPxPyPzE(px_m,py_m,pz_p,e_m);
		ptot=p1+p2;
		trimu_m=ptot.M();
		trimu_pz=ptot.Pz();
		
		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(jpsi_vtxX - origX_m)+p1.Py()*(jpsi_vtxY - origY_m)+p1.Pz()*(jpsi_vtxZ - origZ_m);
		double e = p2.Px()*(jpsi_vtxX - origX_m)+p2.Py()*(jpsi_vtxY - origY_m)+p2.Pz()*(jpsi_vtxZ - origZ_m);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		trimu_vtxX = (jpsi_vtxX + origX_m + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxY = (jpsi_vtxY + origY_m + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxZ = (jpsi_vtxZ + origZ_m + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip=trimu_vtxX-ptot.Px()*trimu_vtxZ/ptot.Pz();
		double yip=trimu_vtxY-ptot.Py()*trimu_vtxZ/ptot.Pz();
		trimu_ip=sqrt(xip*xip+yip*yip);
		trimu_sigip=1.2533*(30 + 60/ptot.Pt());

		trimu_pseudo_tau=trimu_vtxZ*trimu_m/ptot.Pz();
		trimu_pt = ptot.Pt();

		double r2 = trimu_vtxX*trimu_vtxX +trimu_vtxY*trimu_vtxY+trimu_vtxZ*trimu_vtxZ;
		double r = sqrt(r2);
		costhetap = (trimu_vtxX*ptot.Px()+trimu_vtxY*ptot.Py()+trimu_vtxZ*ptot.Pz())/(r*ptot.P());

		double xdca1 = mum_jpsi_vtxX+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PX/mum_jpsi_PZ;
		double ydca1 = mum_jpsi_vtxY+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PY/mum_jpsi_PZ;

		double xdca2 = mup_jpsi_vtxX+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PX/mup_jpsi_PZ;
		double ydca2 = mup_jpsi_vtxY+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PY/mup_jpsi_PZ;;

		double xdca3 = origX_m + (trimu_vtxZ - origZ_m)*px_m/pz_m;
		double ydca3 = origY_m + (trimu_vtxZ - origZ_m)*py_m/pz_m;

		disp = sqrt(xdca1*xdca1 + ydca1*ydca1 + xdca2*xdca2 + ydca2*ydca2 + xdca3*xdca3 + ydca3*ydca3);
		diff_pseudo_tau = abs(jpsi_pseudo_tau-trimu_pseudo_tau);

		charge=-1;
		mu_IP=ip_m;
		mu_sigIP=sig_m;

		if (trimu_m > m_lim_inf && trimu_m < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && PrevDecay(origZ_m)){
			t2.Fill();
		}
	}
	for (Int_t j=0; j<list_hadrons[3*index_ev+1]; j++) {
		t_Dm->GetEntry(j+index_Dm);
		ev=index_ev;
		dep2=abs(Dm - origZ_m_2);

		p2.SetPxPyPzE(px_m_2,py_m_2,pz_m_2,e_m_2);
		ptot=p1+p2;
		trimu_m=ptot.M();
		trimu_pz=ptot.Pz();
		
		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(jpsi_vtxX - origX_m_2)+p1.Py()*(jpsi_vtxY - origY_m_2)+p1.Pz()*(jpsi_vtxZ - origZ_m_2);
		double e = p2.Px()*(jpsi_vtxX - origX_m_2)+p2.Py()*(jpsi_vtxY - origY_m_2)+p2.Pz()*(jpsi_vtxZ - origZ_m_2);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		trimu_vtxX = (jpsi_vtxX + origX_m_2 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxY = (jpsi_vtxY + origY_m_2 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxZ = (jpsi_vtxZ + origZ_m_2 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip=trimu_vtxX-ptot.Px()*trimu_vtxZ/ptot.Pz();
		double yip=trimu_vtxY-ptot.Py()*trimu_vtxZ/ptot.Pz();
		trimu_ip=sqrt(xip*xip+yip*yip);
		trimu_sigip=1.2533*(30 + 60/ptot.Pt());

		trimu_pseudo_tau=trimu_vtxZ*trimu_m/ptot.Pz();
		trimu_pt = ptot.Pt();

		double r2 = trimu_vtxX*trimu_vtxX +trimu_vtxY*trimu_vtxY+trimu_vtxZ*trimu_vtxZ;
		double r = sqrt(r2);
		costhetap = (trimu_vtxX*ptot.Px()+trimu_vtxY*ptot.Py()+trimu_vtxZ*ptot.Pz())/(r*ptot.P());

		double xdca1 = mum_jpsi_vtxX+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PX/mum_jpsi_PZ;
		double ydca1 = mum_jpsi_vtxY+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PY/mum_jpsi_PZ;

		double xdca2 = mup_jpsi_vtxX+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PX/mup_jpsi_PZ;
		double ydca2 = mup_jpsi_vtxY+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PY/mup_jpsi_PZ;;

		double xdca3 = origX_m_2 + (trimu_vtxZ - origZ_m_2)*px_m_2/pz_m_2;
		double ydca3 = origY_m_2 + (trimu_vtxZ - origZ_m_2)*py_m_2/pz_m_2;

		disp = sqrt(xdca1*xdca1 + ydca1*ydca1 + xdca2*xdca2 + ydca2*ydca2 + xdca3*xdca3 + ydca3*ydca3);
		diff_pseudo_tau = abs(jpsi_pseudo_tau-trimu_pseudo_tau);

		charge=-1;
		mu_IP=ip_m_2;
		mu_sigIP=sig_m_2;

		if (trimu_m > m_lim_inf && trimu_m < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && PrevDecay(origZ_m_2)){
			t2.Fill();
		}
	}
	for (Int_t j=0; j<list_hadrons[3*index_ev+2]; j++) {
		t_K_2->GetEntry(j+index_Kb);
		ev=index_ev;
		dep2=abs(Km - origZ_m_3);

		p2.SetPxPyPzE(px_m_3,py_m_3,pz_m_3,e_m_3);
		ptot=p1+p2;
		trimu_m=ptot.M();
		trimu_pz=ptot.Pz();
		
		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(jpsi_vtxX - origX_m_3)+p1.Py()*(jpsi_vtxY - origY_m_3)+p1.Pz()*(jpsi_vtxZ - origZ_m_3);
		double e = p2.Px()*(jpsi_vtxX - origX_m_3)+p2.Py()*(jpsi_vtxY - origY_m_3)+p2.Pz()*(jpsi_vtxZ - origZ_m_3);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		trimu_vtxX = (jpsi_vtxX + origX_m_3 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxY = (jpsi_vtxY + origY_m_3 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		trimu_vtxZ = (jpsi_vtxZ + origZ_m_3 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip=trimu_vtxX-ptot.Px()*trimu_vtxZ/ptot.Pz();
		double yip=trimu_vtxY-ptot.Py()*trimu_vtxZ/ptot.Pz();
		trimu_ip=sqrt(xip*xip+yip*yip);
		trimu_sigip=1.2533*(30 + 60/ptot.Pt());

		trimu_pseudo_tau=trimu_vtxZ*trimu_m/ptot.Pz();
		trimu_pt = ptot.Pt();

		double r2 = trimu_vtxX*trimu_vtxX +trimu_vtxY*trimu_vtxY+trimu_vtxZ*trimu_vtxZ;
		double r = sqrt(r2);
		costhetap = (trimu_vtxX*ptot.Px()+trimu_vtxY*ptot.Py()+trimu_vtxZ*ptot.Pz())/(r*ptot.P());

		double xdca1 = mum_jpsi_vtxX+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PX/mum_jpsi_PZ;
		double ydca1 = mum_jpsi_vtxY+(trimu_vtxZ-mum_jpsi_vtxZ)*mum_jpsi_PY/mum_jpsi_PZ;

		double xdca2 = mup_jpsi_vtxX+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PX/mup_jpsi_PZ;
		double ydca2 = mup_jpsi_vtxY+(trimu_vtxZ-mup_jpsi_vtxZ)*mup_jpsi_PY/mup_jpsi_PZ;;

		double xdca3 = origX_m_3 + (trimu_vtxZ - origZ_m_3)*px_m_3/pz_m_3;
		double ydca3 = origY_m_3 + (trimu_vtxZ - origZ_m_3)*py_m_3/pz_m_3;

		disp = sqrt(xdca1*xdca1 + ydca1*ydca1 + xdca2*xdca2 + ydca2*ydca2 + xdca3*xdca3 + ydca3*ydca3);
		diff_pseudo_tau = abs(jpsi_pseudo_tau-trimu_pseudo_tau);

		charge=-1;
		mu_IP=ip_m_3;
		mu_sigIP=sig_m_3;

		if (trimu_m > m_lim_inf && trimu_m < m_lim_sup && IsinAcc(pt_m_3, eta_m_3, pt_cut) && PrevDecay(origZ_m_3)){
			t2.Fill();
		}
	}
	if ((count+1) >= list_hadrons[3*index_ev]){
		count=0;

		index_Pim+=list_hadrons[3*index_ev];
		index_Dm+=list_hadrons[3*index_ev+1];
		index_Kb+=list_hadrons[3*index_ev+2];

		index_Pip+=blist_hadrons[3*index_ev];
		index_Dp+=blist_hadrons[3*index_ev+1];
		index_K+=blist_hadrons[3*index_ev+2];
		index_ev++;
	}else{
		count++;
	}
   }

   t2.Print();
   muon_out->cd();
   t2.Write();
   muon_out->Close();

}

int main(int argc, char ** argv)
{
	if (argc > 13)
	{
		muon_pairs_bkg(stoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], stod(argv[9]), stod(argv[10]), stod(argv[11]), argv[12], stod(argv[13]));
	}
	return 0;
}
