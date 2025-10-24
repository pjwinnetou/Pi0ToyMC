#include <iostream>
#include <cmath>
#include <memory>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "utils.h"
using namespace std;

struct DecayResult {
  bool accepted = false;
  TLorentzVector pi0;
  TLorentzVector g1;
  TLorentzVector g2;
};

class Pi0ToyMC {
public:
  double aE_=0.13, bE_=0.12, cE_=0.07;
  double aPos_=0.0016, bPos_=0.0032;

  //Pi0ToyMC(TF1* fpt, double yMax, double etaMax, double mPi0, unsigned int seed=0, double ptMin=0., double ptMax=50.)
  Pi0ToyMC(TH1* h1, double yMax, double etaMax, double mPi0,  double aE, double bE, double cE, double aPos, double bPos, unsigned int seed=0)
  : h1_(h1), yMax_(yMax), etaMax_(etaMax), mPi0_(mPi0), aE_(aE), bE_(bE), cE_(cE), aPos_(aPos), bPos_(bPos), rng_(seed ? seed : 0)
  {
    if (!h1_) {
      throw std::runtime_error("Pi0ToyMC: h1 (TH1*) is null!");
      //throw std::runtime_error("Pi0ToyMC: fpt (TF1*) is null!");
    }
  }

  double radius = 93.0;
  double sigma_vz= 0.0;

  const double ptsinglecut_ = 0.5;
  DecayResult GenerateOne() {
    DecayResult out;

    const double pt  = h1_->GetRandom();
    const double y   = Uniform(-yMax_, +yMax_);
    const double phi = Uniform(-M_PI, M_PI);

    TLorentzVector pi0 = MakeMotherFromPtYPhiM(pt, y, phi, mPi0_);

    TLorentzVector g1_rf, g2_rf;
    MakePi0ToGG_Isotropic_RF(g1_rf, g2_rf);

    const TVector3 beta = pi0.BoostVector();
    TLorentzVector g1 = g1_rf; g1.Boost(beta);
    TLorentzVector g2 = g2_rf; g2.Boost(beta);

    TLorentzVector sum = g1 + g2;

    /*if (abs(g1.Eta()) > etaMax_ || abs(g2.Eta()) > etaMax_) {
      out.accepted = false;
      return out;
      }

      if (g1.Pt() < ptsinglecut_ || g2.Pt() < ptsinglecut_) {
      out.accepted = false;
      return out;
      }*/

    out.accepted = true;
    out.pi0 = pi0;
    out.g1  = g1;
    out.g2  = g2;
    return out;
  }


  double def_energy_scale(double E){
    double a  = 1.024;
    double b = 0.081;    
    return (a-0.001/E)*(1-b/sqrt(E));
    //return 1;
  }

  double energy_scale(double E){
    double c = 0.000;
    return 1+c;
  }

  double smear_energy(double E){
    return sqrt(aE_*aE_/E + bE_*bE_/(E*E) + cE_*cE_);
  }

  double smear_position_eta(double E){
    return sqrt(bPos_*bPos_ + aPos_*aPos_/E);
  }

  double smear_position_phi(double E){
    return sqrt(bPos_*bPos_ + aPos_*aPos_/E);
  }

  Double_t GetShiftedEta(float _vz, float _eta)
  {
    double theta = 2*atan(exp(-_eta));
    double z = radius / tan(theta);
    double zshifted = z - _vz;
    double thetashifted = atan2(radius,zshifted);
    double etashifted = -log(tan(thetashifted/2.0));
    return etashifted;
  }

  TLorentzVector GetSmearedPhoton(TLorentzVector truthvec, double vz)
  {
    const double energy = truthvec.E();
    const double phi_tr = truthvec.Phi();
    const double eta_tr = truthvec.Eta();

    
    double smeareta = GetShiftedEta(vz, eta_tr);

    smeareta = rng_.Gaus(smeareta, smear_position_eta(energy));
    double smearphi = rng_.Gaus(phi_tr,  smear_position_phi(energy));

    double smeare = rng_.Gaus(energy, energy*smear_energy(energy));
    smeare = smeare * def_energy_scale(energy) * energy_scale(smeare) * 1.02;

    double smearpt = smeare / cosh(smeareta);
    if (smearpt < 0) smearpt = 0;

    TLorentzVector vsm;
    vsm.SetPtEtaPhiM(smearpt, smeareta, smearphi, 0);
    return vsm;
  }

  TRandom3 rng_;

private:
  TF1* fpt_ = nullptr;
  TH1* h1_ = nullptr;
  double yMax_ = 1.0;
  double etaMax_ = 1.0;
  double mPi0_ = 0.1349768;

  double Uniform(double a, double b) { return a + (b - a) * rng_.Rndm(); }

  static TLorentzVector MakeMotherFromPtYPhiM(double pt, double y, double phi, double m) {
    double mT = sqrt(m*m + pt*pt);
    double pz = mT * sinh(y);
    double p = sqrt(pz*pz + pt*pt);
    double eta = atanh(pz/p);
    double energy = sqrt(m*m + p*p);
    TLorentzVector v;
    v.SetPtEtaPhiE(pt, eta, phi, energy);
    return v;
  }

  void MakePi0ToGG_Isotropic_RF(TLorentzVector& g1_rf, TLorentzVector& g2_rf) {
    const double E = mPi0_ / 2.0;
    const double p = E;

    const double cosTheta = Uniform(-1.0, 1.0);
    const double sinTheta = sqrt(1 - cosTheta*cosTheta);
    const double phi_d    = Uniform(-M_PI, M_PI);

    const double px = p * sinTheta * cos(phi_d);
    const double py = p * sinTheta * sin(phi_d);
    const double pz = p * cosTheta;

    g1_rf.SetPxPyPzE( px, py, pz, E );
    g2_rf.SetPxPyPzE( -px, -py, -pz, E );
  }
};

int main(int argc, char** argv) {
  long long nEvents = 100;
  double yMax = 0.4;
  double etaMax = 0.4;
  unsigned int seed = 0;
  double mPi0 = 0.1349;
  double ptMin = 0.0;
  double ptMax = 10.0;

  int param_index = 0;

  if (argc > 1) nEvents = std::stoll(argv[1]);
  if (argc > 2) seed    = static_cast<unsigned int>(std::stoul(argv[2]));
  if (argc > 3) yMax    = std::stod(argv[3]);
  if (argc > 4) etaMax  = std::stod(argv[4]);
  if (argc > 5) param_index = std::stoi(argv[5]);

  double aE, bE, cE, aPos, bPos;
  GetParameterSet(param_index, aE, bE, cE, aPos, bPos);

  std::cout << "Configuration setup " << std::endl;
  std::cout << "  nEvents = " << nEvents << std::endl;
  std::cout << "  seed    = " << seed    << "\n"
            << "  yMax    = " << yMax    << "\n"
            << "  etaMax  = " << etaMax  << "\n";
  
  std::cout << "------------ Parameter set -------------" << std::endl;
  std::cout << "Index " << param_index << ": aE=" << aE << " bE=" << bE << " cE=" << cE << " aPos=" << aPos << " bPos=" << bPos << std::endl;

  //TFile *fit_file = new TFile("fit_pythiamb_cross_section_pi0.root","read");
  //TF1 *fpt = (TF1*) fit_file->Get("fit_restricted");
  //fpt->SetRange(0,10);

  TFile *inputhistfile = new TFile("hist_pt0.0_all.root","read");
  TH1D *h1 = (TH1D*) inputhistfile->Get("hpt");

  auto fout = std::make_unique<TFile>(Form("pi0_toymc_index%d.root",param_index), "RECREATE");
  TH1D h_pt_pi0("h_pt_pi0", ";p_{T}^{#pi^{0}} [GeV];Events", 200, ptMin, ptMax);
  TH1D h_y_pi0 ("h_y_pi0",  ";y^{#pi^{0}};Events",          120, -yMax, yMax);

  TH1D h_alpha ("h_alpha", ";alpha;",100,0,1);
  TH1D h_eta_g ("h_eta_g",  ";#eta^{#gamma};Photons",       120, -etaMax, etaMax);
  TH1D h_pt_g  ("h_pt_g",   ";p_{T}^{#gamma} [GeV];Photons",200, 0, ptMax);

  TH1D h_dphi("h_dphi",";dphi;",100,-3*M_PI,3*M_PI);
  TH2D h_pt1pt2("h_pt1pt2",";pt1;pt2",100,0,5,100,0,5);
  TH2D h_dRpt("h_dRpt",";pt;dR",100,0,5,100,0,1);

  TH1D h_m_resid("h_m_resid",";(m_{#gamma1+#gamma2}-E_{#pi^{0}}) [GeV];Events", 200, -1e-9, 1e-9);
  TH1D h_mass("h_mass",";m_{#gamma1+#gamma2} [GeV];",100,0.05,0.25);
  TH1D h_mass_smeared_tot("h_mass_smeared_tot",";m_{#gamma1+#gamma2} [GeV];",200,0.0,0.35);
  TH1D *h_mass_smeared[nPtBins][nAlphaBins];
  for(int i=0; i< nPtBins; i++){
    for(int j=0; j<nAlphaBins ; j++){
      h_mass_smeared[i][j] = new TH1D(Form("h_mass_smeared_ptbin%d_alphabin%d",i,j),";m_{#gamma1+#gamma2} [GeV];",200,0.0,0.35);
    }
  }

  long long nTried = 0;
  long long nAccepted = 0;

  float singleptcut = 0.5;
  float singleetacut = 0.3;

  Pi0ToyMC gen(h1, yMax, etaMax, mPi0, aE, bE, cE, aPos, bPos, seed);

  TRandom3 rndm(seed);
  long long binCounts[nPtBins][nAlphaBins] = {0};
  long long totalAccepted = 0;
  const long long targetPerBin = nEvents;

  while(totalAccepted < nPtBins * nAlphaBins * targetPerBin)
  {
    ++nTried;
    double vz = gen.rng_.Gaus(0.0, gen.sigma_vz);
    DecayResult ev = gen.GenerateOne();
    if (!ev.accepted) continue;

    double pt = ev.pi0.Pt();
    double y = ev.pi0.Rapidity();
    double pt1 = ev.g1.Pt();
    double pt2 = ev.g2.Pt();
    double eta1 = ev.g1.Eta();
    double eta2 = ev.g2.Eta();
    double phi1 = ev.g1.Phi();
    double phi2 = ev.g2.Phi();
    double dphi = phi2 - phi1;
    double dR = sqrt(dphi*dphi + (eta2-eta1)*(eta2-eta1));

    float alpha = abs(pt1-pt2) / (pt1+pt2);
    h_alpha.Fill(alpha);
        
    h_pt_pi0.Fill(pt);
    h_y_pi0 .Fill(y);

    h_eta_g.Fill(eta1);
    h_eta_g.Fill(eta2);
    h_pt_g.Fill(pt1);
    h_pt_g.Fill(pt2);
    h_dphi.Fill(dphi);

    TLorentzVector sum = ev.g1 + ev.g2;
    h_m_resid.Fill(sum.M()  - ev.pi0.M());
    h_pt1pt2.Fill(ev.g1.Pt(), ev.g2.Pt());
    h_dRpt.Fill(ev.pi0.Pt(), dR); 
    h_mass.Fill(sum.M());

    TLorentzVector vsm1 = gen.GetSmearedPhoton(ev.g1,vz);
    TLorentzVector vsm2 = gen.GetSmearedPhoton(ev.g2,vz);

    if(vsm1.Pt() < singleptcut || vsm2.Pt() < singleptcut) continue;
    if(abs(vsm1.Eta()) > singleetacut || abs(vsm2.Eta()) > singleetacut) continue;
    TLorentzVector vsm = vsm1 + vsm2;
    //if(vsm.Pt()<2*1.2*singleptcut) continue;
    float smpt = vsm.Pt();
    if(smpt< 1) continue;
    
    float sm_mass = vsm.M();
    h_mass_smeared_tot.Fill(sm_mass);

    int thisptbin = GetPtBin(smpt);
    alpha = (vsm1.E()-vsm2.E())/(vsm1.E() +vsm2.E());
    int thisalphabin = GetAlphaBin(thisptbin,alpha);
    if(thisptbin == -1 || thisalphabin == -1) continue;

    if(binCounts[thisptbin][thisalphabin] < targetPerBin) {
      h_mass_smeared[thisptbin][thisalphabin]->Fill(sm_mass);
      binCounts[thisptbin][thisalphabin]++;
      totalAccepted++;
    }
  }

  std::cout << "\nTried: " << nTried
            << "  Accepted (both |eta_g| <= " << etaMax << "): " << nAccepted
            << "  Efficiency: " << (nTried ? (double)nAccepted / nTried : 0.0) << std::endl;

  for(int i=0; i< nPtBins; i++){
    for(int j=0; j<nAlphaBins ; j++){
      h_mass_smeared[i][j]->Write();
    }
  }
  fout->Write();
  fout->Close();
  std::cout << "Done.. " << std::endl;

  return 0;
}
