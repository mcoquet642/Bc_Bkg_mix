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

int sgn(double d){
    double eps=1e-5;
    return d<-eps?-1:d>eps;
}

bool IsinAcc(double pt, double eta, double ptcut){
	if (pt > ptcut && eta > 2.5 && eta < 3.6)
	{
		return true;
	}else{
		return false;}
}

bool PrevDecay(double decayZ){
	double MFT_f0Z = 460*1000;
	return (decayZ < MFT_f0Z);
}

void muon_pairs_bkg(int nb_events, char * output, char * inputPip, char * inputPim, char * inputDp, char * inputDm, char * inputKp, char * inputKm, double m_lim_inf, double m_lim_sup, double pt_cut){
   TFile *f = new TFile(inputPip);
   TTree *t1 = (TTree*)f->Get("DecayTree");
   TFile *f_2 = new TFile(inputPim);
   TTree *t1_2 = (TTree*)f_2->Get("DecayTree");
   Double_t px_m, py_m, pz_m, e_m, ip_m, sig_m, origZ_m, Pim, pt_m, eta_m, origX_m, origY_m;
   Double_t px_p, py_p, pz_p, e_p, ip_p, sig_p, origZ_p, Pip, pt_p, eta_p, origX_p, origY_p;

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
   TBranch *br1_sigipm = t1_2->GetBranch("mum_SIGMAIP_TRUE");
   TBranch *br1_sigipp = t1->GetBranch("mup_SIGMAIP_TRUE");
   TBranch *br1_origZm = t1_2->GetBranch("mum_origZ");
   TBranch *br1_origZp = t1->GetBranch("mup_origZ");
   TBranch *br1_origYm = t1_2->GetBranch("mum_origY");
   TBranch *br1_origYp = t1->GetBranch("mup_origY");
   TBranch *br1_origXm = t1_2->GetBranch("mum_origX");
   TBranch *br1_origXp = t1->GetBranch("mup_origX");
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
   br1_sigipm->SetAddress(&sig_m);
   br1_sigipp->SetAddress(&sig_p);
   br1_origZm->SetAddress(&origZ_m);
   br1_origZp->SetAddress(&origZ_p);
   br1_origYm->SetAddress(&origY_m);
   br1_origYp->SetAddress(&origY_p);
   br1_origXm->SetAddress(&origX_m);
   br1_origXp->SetAddress(&origX_p);
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
   Double_t px_m_2, py_m_2, pz_m_2, e_m_2, ip_m_2, sig_m_2, origZ_m_2, Dm, pt_m_2, eta_m_2, origX_m_2, origY_m_2;
   Double_t px_p_2, py_p_2, pz_p_2, e_p_2, ip_p_2, sig_p_2, origZ_p_2, Dp, pt_p_2, eta_p_2, origX_p_2, origY_p_2;


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
   TBranch *br2_sigipm = t_Dm->GetBranch("mum_SIGMAIP_TRUE");
   TBranch *br2_sigipp = t_Dp->GetBranch("mup_SIGMAIP_TRUE");
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
   br2_sigipm->SetAddress(&sig_m_2);
   br2_sigipp->SetAddress(&sig_p_2);
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
   Double_t px_m_3, py_m_3, pz_m_3, e_m_3, ip_m_3, sig_m_3, origZ_m_3, Km, pt_m_3, eta_m_3, origX_m_3, origY_m_3;
   Double_t px_p_3, py_p_3, pz_p_3, e_p_3, ip_p_3, sig_p_3, origZ_p_3, Kp, pt_p_3, eta_p_3, origX_p_3, origY_p_3;


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
   TBranch *br3_sigipm = t_K_2->GetBranch("mum_SIGMAIP_TRUE");
   TBranch *br3_sigipp = t_K->GetBranch("mup_SIGMAIP_TRUE");
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
   br3_sigipm->SetAddress(&sig_m_3);
   br3_sigipp->SetAddress(&sig_p_3);
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


//Decay Probabilities_______________

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

//Partie combinatoire

   TRandom3 r1, r2, r3, r4, r5, r6;

   cout << "Number of events : " << nb_events << endl;


//_________





   Double_t m_inv, e, pt, y, ev, dep_min, dep_max, dep_av, px, py, pz, ip, sigip, origx, origy;
   Double_t mum_jpsi_pt, mum_jpsi_eta, mum_jpsi_sig, mum_jpsi_ip, mum_jpsi_vtxX, mum_jpsi_vtxY, mum_jpsi_vtxZ, mum_jpsi_px, mum_jpsi_py, mum_jpsi_pz;
   Double_t mup_jpsi_pt, mup_jpsi_eta, mup_jpsi_sig, mup_jpsi_ip, mup_jpsi_vtxX, mup_jpsi_vtxY, mup_jpsi_vtxZ, mup_jpsi_px, mup_jpsi_py, mup_jpsi_pz;
   Double_t dep1, dep2;
   gROOT->cd();
   TFile * muon_out = new TFile(output, "recreate");
   TTree t2("DecayTree","dimuon Tree");

   t2.Branch("Jpsi_M",&m_inv,"Jpsi_M/D");
   t2.Branch("Jpsi_PT",&pt,"Jpsi_PT/D");
   t2.Branch("Jpsi_PX",&px,"Jpsi_PX/D");
   t2.Branch("Jpsi_PY",&py,"Jpsi_PY/D");
   t2.Branch("Jpsi_PZ",&pz,"Jpsi_PZ/D");
   t2.Branch("Jpsi_E",&e,"Jpsi_E/D");
   t2.Branch("Jpsi_y",&y,"Jpsi_y/D");
   t2.Branch("Jpsi_origZ_max",&dep_max,"Jpsi_origZ_max/D");
   t2.Branch("Jpsi_origZ_min",&dep_min,"Jpsi_origZ_min/D");
   t2.Branch("Jpsi_origX",&origx,"Jpsi_origX/D");
   t2.Branch("Jpsi_origY",&origy,"Jpsi_origY/D");
   t2.Branch("Jpsi_origZ",&dep_av,"Jpsi_origZ/D");
   t2.Branch("Jpsi_IP",&ip,"Jpsi_IP/D");
   t2.Branch("Jpsi_SIGMAIP_TRUE",&sigip,"Jpsi_SIGMAIP_TRUE/D");
   t2.Branch("events",&ev,"ev/D");

   t2.Branch("mum_jpsi_PT",&mum_jpsi_pt);
   t2.Branch("mum_jpsi_eta",&mum_jpsi_eta);
   t2.Branch("mum_jpsi_SIGMAIP_TRUE",&mum_jpsi_sig);
   t2.Branch("mum_jpsi_IP",&mum_jpsi_ip);
   t2.Branch("mum_jpsi_origX",&mum_jpsi_vtxX);
   t2.Branch("mum_jpsi_origY",&mum_jpsi_vtxY);
   t2.Branch("mum_jpsi_origZ",&mum_jpsi_vtxZ);
   t2.Branch("mum_jpsi_PX",&mum_jpsi_px);
   t2.Branch("mum_jpsi_PY",&mum_jpsi_py);
   t2.Branch("mum_jpsi_PZ",&mum_jpsi_pz);

   t2.Branch("mup_jpsi_PT",&mup_jpsi_pt);
   t2.Branch("mup_jpsi_eta",&mup_jpsi_eta);
   t2.Branch("mup_jpsi_SIGMAIP_TRUE",&mup_jpsi_sig);
   t2.Branch("mup_jpsi_IP",&mup_jpsi_ip);
   t2.Branch("mup_jpsi_origX",&mup_jpsi_vtxX);
   t2.Branch("mup_jpsi_origY",&mup_jpsi_vtxY);
   t2.Branch("mup_jpsi_origZ",&mup_jpsi_vtxZ);
   t2.Branch("mup_jpsi_PX",&mup_jpsi_px);
   t2.Branch("mup_jpsi_PY",&mup_jpsi_py);
   t2.Branch("mup_jpsi_PZ",&mup_jpsi_pz);

   TLorentzVector p1;
   TLorentzVector p2;
   TLorentzVector ptot;

   Int_t nb_Pip=0;
   Int_t nb_Dp=0;
   Int_t nb_K=0;

   Int_t nb_muon_Pip=0;
   Int_t nb_muon_Dp=0;
   Int_t nb_muon_K=0;
   Int_t nb_muon_Pim=0;
   Int_t nb_muon_Dm=0;
   Int_t nb_muon_Kb=0;
   Int_t somme = 0;

   Int_t list_hadrons[3*nb_events];
   Int_t blist_hadrons[3*nb_events];
   Int_t buff=0;


   for (Int_t k=0; k<nb_events; k++){
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

	buff = (list_hadrons[3*k]+list_hadrons[3*k+1]+list_hadrons[3*k+2])*(blist_hadrons[3*k]+blist_hadrons[3*k+1]+blist_hadrons[3*k+2]);
	somme += buff;
  	if (k%10000 == 0){
		cout << " Event : " << k << " , Number of Pip : " << list_hadrons[3*k] << " , Number of Dp : " << list_hadrons[3*k+1] <<  " , Number of K : " << list_hadrons[3*k+2] << " , Number of Pim : " << blist_hadrons[3*k] << " , Number of Dm : " << blist_hadrons[3*k+1] <<  " , Number of Kb : " << blist_hadrons[3*k+2] << " , Number of pairs : " << buff << endl;
  	}
   }

   cout << "Total number of Pip muons : " << nb_muon_Pip << ", Total number of Dp muons : " << nb_muon_Dp << ", Total number of K muons : " << nb_muon_K << "Total number of Pim muons : " << nb_muon_Pim << ", Total number of Dm muons : " << nb_muon_Dm << ", Total number of Kb muons : " << nb_muon_Kb << ", Total number of pairs : " << somme << endl;

   if (nb_muon_Pip > entries_min_Pip || nb_muon_Dp > entries_min_Dp || nb_muon_K > entries_min_K){
	cout << "Note enough muons provided." << endl;
   }

   Int_t count = 0;
   Int_t index_ev = 0;
   Int_t index_Pip = 0;
   Int_t index_Dp = 0;
   Int_t index_K = 0;

   Int_t nb_final=0;
   Int_t nb_init=0;


   double l1, l2;
   double theta1, theta2;
   double alpha1, alpha2;
   double vtxr1, vtxr2;
   double dist1, dist2;


   for (Int_t i=0; i<nb_muon_Pip; i++){
	while (list_hadrons[3*index_ev]==0){
		index_Dp+=blist_hadrons[3*index_ev+1];
		index_K+=blist_hadrons[3*index_ev+2];
		index_ev++;
	}
	t1_2->GetEntry(i);
	dep1=abs(Pim - origZ_m);
	mum_jpsi_pt = pt_m;
	mum_jpsi_eta = eta_m;
	mum_jpsi_ip = ip_m;
	mum_jpsi_sig = sig_m;
	mum_jpsi_vtxX = origX_m;
	mum_jpsi_vtxY = origY_m;
	mum_jpsi_vtxZ = origZ_m;

	mum_jpsi_px = px_m;
	mum_jpsi_py = py_m;
	mum_jpsi_pz = pz_m;

	p1.SetPxPyPzE(px_m,py_m,pz_m,e_m);
	double a = p1.P()*p1.P();

   	for (Int_t j=0; j<blist_hadrons[3*index_ev]; j++) {
  		t1->GetEntry(j+i-count);
		p2.SetPxPyPzE(px_p,py_p,pz_p,e_p);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Pip - origZ_p);
		mup_jpsi_pt = pt_p;
		mup_jpsi_eta = eta_p;
		mup_jpsi_ip = ip_p;
		mup_jpsi_sig = sig_p;
		mup_jpsi_vtxX = origX_p;
		mup_jpsi_vtxY = origY_p;
		mup_jpsi_vtxZ = origZ_p;

		mup_jpsi_px = px_p;
		mup_jpsi_py = py_p;
		mup_jpsi_pz = pz_p;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m - origX_p)+p1.Py()*(origY_m - origY_p)+p1.Pz()*(origZ_m - origZ_p);
		double e = p2.Px()*(origX_m - origX_p)+p2.Py()*(origY_m - origY_p)+p2.Pz()*(origZ_m - origZ_p);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m + origX_p + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m + origY_p + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m + origZ_p + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && IsinAcc(pt_p, eta_p, pt_cut) && PrevDecay(origZ_p) && PrevDecay(origZ_m)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && IsinAcc(pt_p, eta_p, pt_cut)){
			nb_init++;
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+1]; j++) {
	  	t_Dp->GetEntry(j+index_Dp);
		p2.SetPxPyPzE(px_p_2,py_p_2,pz_p_2,e_p_2);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Dp - origZ_p_2);
		mup_jpsi_pt = pt_p_2;
		mup_jpsi_eta = eta_p_2;
		mup_jpsi_ip = ip_p_2;
		mup_jpsi_sig = sig_p_2;
		mup_jpsi_vtxX = origX_p_2;
		mup_jpsi_vtxY = origY_p_2;
		mup_jpsi_vtxZ = origZ_p_2;

		mup_jpsi_px = px_p_2;
		mup_jpsi_py = py_p_2;
		mup_jpsi_pz = pz_p_2;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m - origX_p_2)+p1.Py()*(origY_m - origY_p_2)+p1.Pz()*(origZ_m - origZ_p_2);
		double e = p2.Px()*(origX_m - origX_p_2)+p2.Py()*(origY_m - origY_p_2)+p2.Pz()*(origZ_m - origZ_p_2);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m + origX_p_2 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m + origY_p_2 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m + origZ_p_2 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && IsinAcc(pt_p_2, eta_p_2, pt_cut) && PrevDecay(origZ_p_2) && PrevDecay(origZ_m)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && IsinAcc(pt_p_2, eta_p_2, pt_cut)){
			nb_init++;
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+2]; j++) {
	  	t_K->GetEntry(j+index_K);
		p2.SetPxPyPzE(px_p_3,py_p_3,pz_p_3,e_p_3);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Kp - origZ_p_3);
		mup_jpsi_pt = pt_p_3;
		mup_jpsi_eta = eta_p_3;
		mup_jpsi_ip = ip_p_3;
		mup_jpsi_sig = sig_p_3;
		mup_jpsi_vtxX = origX_p_3;
		mup_jpsi_vtxY = origY_p_3;
		mup_jpsi_vtxZ = origZ_p_3;

		mup_jpsi_px = px_p_3;
		mup_jpsi_py = py_p_3;
		mup_jpsi_pz = pz_p_3;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m - origX_p_3)+p1.Py()*(origY_m - origY_p_3)+p1.Pz()*(origZ_m - origZ_p_3);
		double e = p2.Px()*(origX_m - origX_p_3)+p2.Py()*(origY_m - origY_p_3)+p2.Pz()*(origZ_m - origZ_p_3);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m + origX_p_3 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m + origY_p_3 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m + origZ_p_3 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && IsinAcc(pt_p_3, eta_p_3, pt_cut) && PrevDecay(origZ_p_3) && PrevDecay(origZ_m)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m, eta_m, pt_cut) && IsinAcc(pt_p_3, eta_p_3, pt_cut)){
			nb_init++;
		}
	}
	if ((count+1) >= list_hadrons[3*index_ev]){
		count=0;
		index_Dp+=blist_hadrons[3*index_ev+1];
		index_K+=blist_hadrons[3*index_ev+2];
		index_ev++;
	}else{
		count++;
	}
   }
   count=0;
   index_ev=0;
   index_K=0;

   for (Int_t i=0; i<nb_muon_Dp; i++){
	while (list_hadrons[3*index_ev+1]==0){
		index_Pip+=blist_hadrons[3*index_ev];
		index_K+=blist_hadrons[3*index_ev+2];
		index_ev++;
	}
	t_Dm->GetEntry(i);
	dep1=abs(Dm - origZ_m_2);
	mum_jpsi_pt = pt_m_2;
	mum_jpsi_eta = eta_m_2;
	mum_jpsi_ip = ip_m_2;
	mum_jpsi_sig = sig_m_2;
	mum_jpsi_vtxX = origX_m_2;
	mum_jpsi_vtxY = origY_m_2;
	mum_jpsi_vtxZ = origZ_m_2;

	mum_jpsi_px = px_m_2;
	mum_jpsi_py = py_m_2;
	mum_jpsi_pz = pz_m_2;

	p1.SetPxPyPzE(px_m_2,py_m_2,pz_m_2,e_m_2);
	double a = p1.P()*p1.P();

   	for (Int_t j=0; j<blist_hadrons[3*index_ev]; j++) {
  		t1->GetEntry(j+index_Pip);
		p2.SetPxPyPzE(px_p,py_p,pz_p,e_p);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Pip - origZ_p);
		mup_jpsi_pt = pt_p;
		mup_jpsi_eta = eta_p;
		mup_jpsi_ip = ip_p;
		mup_jpsi_sig = sig_p;
		mup_jpsi_vtxX = origX_p;
		mup_jpsi_vtxY = origY_p;
		mup_jpsi_vtxZ = origZ_p;

		mup_jpsi_px = px_p;
		mup_jpsi_py = py_p;
		mup_jpsi_pz = pz_p;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m_2 - origX_p)+p1.Py()*(origY_m_2 - origY_p)+p1.Pz()*(origZ_m_2 - origZ_p);
		double e = p2.Px()*(origX_m_2 - origX_p)+p2.Py()*(origY_m_2 - origY_p)+p2.Pz()*(origZ_m_2 - origZ_p);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m_2 + origX_p + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m_2 + origY_p + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m_2 + origZ_p + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && IsinAcc(pt_p, eta_p, pt_cut) && PrevDecay(origZ_p) && PrevDecay(origZ_m_2)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && IsinAcc(pt_p, eta_p, pt_cut)){
			nb_init++;
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+1]; j++) {
	  	t_Dp->GetEntry(j+i-count);
		p2.SetPxPyPzE(px_p_2,py_p_2,pz_p_2,e_p_2);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Dp - origZ_p_2);
		mup_jpsi_pt = pt_p_2;
		mup_jpsi_eta = eta_p_2;
		mup_jpsi_ip = ip_p_2;
		mup_jpsi_sig = sig_p_2;
		mup_jpsi_vtxX = origX_p_2;
		mup_jpsi_vtxY = origY_p_2;
		mup_jpsi_vtxZ = origZ_p_2;

		mup_jpsi_px = px_p_2;
		mup_jpsi_py = py_p_2;
		mup_jpsi_pz = pz_p_2;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m_2 - origX_p_2)+p1.Py()*(origY_m_2 - origY_p_2)+p1.Pz()*(origZ_m_2 - origZ_p_2);
		double e = p2.Px()*(origX_m_2 - origX_p_2)+p2.Py()*(origY_m_2 - origY_p_2)+p2.Pz()*(origZ_m_2 - origZ_p_2);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m_2 + origX_p_2 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m_2 + origY_p_2 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m_2 + origZ_p_2 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && IsinAcc(pt_p_2, eta_p_2, pt_cut) && PrevDecay(origZ_p_2) && PrevDecay(origZ_m_2)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && IsinAcc(pt_p_2, eta_p_2, pt_cut)){
			nb_init++;
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+2]; j++) {
	  	t_K->GetEntry(j+index_K);
		p2.SetPxPyPzE(px_p_3,py_p_3,pz_p_3,e_p_3);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Kp - origZ_p_3);
		mup_jpsi_pt = pt_p_3;
		mup_jpsi_eta = eta_p_3;
		mup_jpsi_ip = ip_p_3;
		mup_jpsi_sig = sig_p_3;
		mup_jpsi_vtxX = origX_p_3;
		mup_jpsi_vtxY = origY_p_3;
		mup_jpsi_vtxZ = origZ_p_3;

		mup_jpsi_px = px_p_3;
		mup_jpsi_py = py_p_3;
		mup_jpsi_pz = pz_p_3;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m_2 - origX_p_3)+p1.Py()*(origY_m_2 - origY_p_3)+p1.Pz()*(origZ_m_2 - origZ_p_3);
		double e = p2.Px()*(origX_m_2 - origX_p_3)+p2.Py()*(origY_m_2 - origY_p_3)+p2.Pz()*(origZ_m_2 - origZ_p_3);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m_2 + origX_p_3 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m_2 + origY_p_3 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m_2 + origZ_p_3 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && IsinAcc(pt_p_3, eta_p_3, pt_cut) && PrevDecay(origZ_p_3) && PrevDecay(origZ_m_2)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_2, eta_m_2, pt_cut) && IsinAcc(pt_p_3, eta_p_3, pt_cut)){
			nb_init++;
		}
	}
	if ((count+1) >= list_hadrons[3*index_ev+1]){
		count=0;
		index_Pip+=blist_hadrons[3*index_ev];
		index_K+=blist_hadrons[3*index_ev+2];
		index_ev++;
	}else{
		count++;
	}
  } 
   count=0;
   index_ev=0;
   index_Pip=0;
   index_Dp=0;

   for (Int_t i=0; i<nb_muon_K; i++){
	while (list_hadrons[3*index_ev+2]==0){
		index_Pip+=blist_hadrons[3*index_ev];
		index_Dp+=blist_hadrons[3*index_ev+1];
		index_ev++;
	}
	t_K_2->GetEntry(i);
	dep1=abs(Km - origZ_m_3);
	mum_jpsi_pt = pt_m_3;
	mum_jpsi_eta = eta_m_3;
	mum_jpsi_ip = ip_m_3;
	mum_jpsi_sig = sig_m_3;
	mum_jpsi_vtxX = origX_m_3;
	mum_jpsi_vtxY = origY_m_3;
	mum_jpsi_vtxZ = origZ_m_3;

	mum_jpsi_px = px_m_3;
	mum_jpsi_py = py_m_3;
	mum_jpsi_pz = pz_m_3;

	p1.SetPxPyPzE(px_m_3,py_m_3,pz_m_3,e_m_3);
	double a = p1.P()*p1.P();

   	for (Int_t j=0; j<blist_hadrons[3*index_ev]; j++) {
  		t1->GetEntry(j+index_Pip);
		p2.SetPxPyPzE(px_p,py_p,pz_p,e_p);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Pip - origZ_p);
		mup_jpsi_pt = pt_p;
		mup_jpsi_eta = eta_p;
		mup_jpsi_ip = ip_p;
		mup_jpsi_sig = sig_p;
		mup_jpsi_vtxX = origX_p;
		mup_jpsi_vtxY = origY_p;
		mup_jpsi_vtxZ = origZ_p;

		mup_jpsi_px = px_p;
		mup_jpsi_py = py_p;
		mup_jpsi_pz = pz_p;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m_3 - origX_p)+p1.Py()*(origY_m_3 - origY_p)+p1.Pz()*(origZ_m_3 - origZ_p);
		double e = p2.Px()*(origX_m_3 - origX_p)+p2.Py()*(origY_m_3 - origY_p)+p2.Pz()*(origZ_m_3 - origZ_p);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m_3 + origX_p + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m_3 + origY_p + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m_3 + origZ_p + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_3, eta_m_3, pt_cut) && IsinAcc(pt_p, eta_p, pt_cut) && PrevDecay(origZ_p) && PrevDecay(origZ_m_3)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_3, eta_m_3, pt_cut) && IsinAcc(pt_p, eta_p, pt_cut)){
			nb_init++;
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+1]; j++) {
	  	t_Dp->GetEntry(j+index_Dp);
		p2.SetPxPyPzE(px_p_2,py_p_2,pz_p_2,e_p_2);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Dp - origZ_p_2);
		mup_jpsi_pt = pt_p_2;
		mup_jpsi_eta = eta_p_2;
		mup_jpsi_ip = ip_p_2;
		mup_jpsi_sig = sig_p_2;
		mup_jpsi_vtxX = origX_p_2;
		mup_jpsi_vtxY = origY_p_2;
		mup_jpsi_vtxZ = origZ_p_2;

		mup_jpsi_px = px_p_2;
		mup_jpsi_py = py_p_2;
		mup_jpsi_pz = pz_p_2;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m_3 - origX_p_2)+p1.Py()*(origY_m_3 - origY_p_2)+p1.Pz()*(origZ_m_3 - origZ_p_2);
		double e = p2.Px()*(origX_m_3 - origX_p_2)+p2.Py()*(origY_m_3 - origY_p_2)+p2.Pz()*(origZ_m_3 - origZ_p_2);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m_3 + origX_p_2 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m_3 + origY_p_2 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m_3 + origZ_p_2 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup  && IsinAcc(pt_m_3, eta_m_3, pt_cut) && IsinAcc(pt_p_2, eta_p_2, pt_cut) && PrevDecay(origZ_p_2) && PrevDecay(origZ_m_3)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup  && IsinAcc(pt_m_3, eta_m_3, pt_cut) && IsinAcc(pt_p_2, eta_p_2, pt_cut)){
			nb_init++;
		}
	}
	for (Int_t j=0; j<blist_hadrons[3*index_ev+2]; j++) {
	  	t_K->GetEntry(j+i-count);
		p2.SetPxPyPzE(px_p_3,py_p_3,pz_p_3,e_p_3);
		ptot=p1+p2;
		m_inv=ptot.M();
		px=ptot.Px();
		py=ptot.Py();
		pz=ptot.Pz();
		e=ptot.E();
		pt=ptot.Pt();
		y=log((e+pz)/(e-pz))*0.5;
		ev=index_ev;
		dep2=abs(Kp - origZ_p_3);
		mup_jpsi_pt = pt_p_3;
		mup_jpsi_eta = eta_p_3;
		mup_jpsi_ip = ip_p_3;
		mup_jpsi_sig = sig_p_3;
		mup_jpsi_vtxX = origX_p_3;
		mup_jpsi_vtxY = origY_p_3;
		mup_jpsi_vtxZ = origZ_p_3;

		mup_jpsi_px = px_p_3;
		mup_jpsi_py = py_p_3;
		mup_jpsi_pz = pz_p_3;

		double b = p1.Px()*p2.Px()+p1.Py()*p2.Py()+p1.Pz()*p2.Pz();
		double c = p2.P()*p2.P();
		double d = p1.Px()*(origX_m_3 - origX_p_3)+p1.Py()*(origY_m_3 - origY_p_3)+p1.Pz()*(origZ_m_3 - origZ_p_3);
		double e = p2.Px()*(origX_m_3 - origX_p_3)+p2.Py()*(origY_m_3 - origY_p_3)+p2.Pz()*(origZ_m_3 - origZ_p_3);

		double s = (b*e-c*d)/(a*c-b*b);
		double t = (a*e-b*d)/(a*c-b*b);

		origx = (origX_m_3 + origX_p_3 + s*p1.Px() + t*p2.Px())/2+gRandom->Gaus(0,0.48)*1000;
		origy = (origY_m_3 + origY_p_3 + s*p1.Py() + t*p2.Py())/2+gRandom->Gaus(0,0.48)*1000;
		dep_av = (origZ_m_3 + origZ_p_3 + s*p1.Pz() + t*p2.Pz())/2+gRandom->Gaus(0,0.48)*1000;

		double xip = origx - px*dep_av/pz;
		double yip = origy - py*dep_av/pz;
		ip = sqrt(xip*xip+yip*yip);
		sigip=1.2533*(30 + 60/pt);
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_3, eta_m_3, pt_cut) && IsinAcc(pt_p_3, eta_p_3, pt_cut) && PrevDecay(origZ_p_3) && PrevDecay(origZ_m_3)){
			t2.Fill();
			nb_final++;
		}
		if (m_inv > m_lim_inf && m_inv < m_lim_sup && IsinAcc(pt_m_3, eta_m_3, pt_cut) && IsinAcc(pt_p_3, eta_p_3, pt_cut)){
			nb_init++;
		}
	}
	if ((count+1) >= list_hadrons[3*index_ev+2]){
		count=0;
		index_Pip+=blist_hadrons[3*index_ev];
		index_Dp+=blist_hadrons[3*index_ev+1];
		index_ev++;
	}else{
		count++;
	}
  } 
   cout << "Efficiency : " << (double)nb_final / (double)nb_init << endl;

   muon_out->cd();
   t2.Write();
   muon_out->Close();

}

int main(int argc, char ** argv)
{
	if (argc > 10)
	{
		muon_pairs_bkg(stoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], stod(argv[9]), stod(argv[10]), stod(argv[11]));
	}
	return 0;
}
