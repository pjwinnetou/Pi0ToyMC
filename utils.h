#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>
#include <stdexcept>

const int nPtBins = 6;
float ptbin[nPtBins+1] = {1, 1.5, 2, 2.5, 3.5, 5, 8};
const int nAlphaBins = 4;
float alphabin[nPtBins][nAlphaBins+1] = { {0,0.045, 0.105, 0.195, 1}, {0, 0.105, 0.215, 0.325, 1}, {0, 0.145, 0.295, 0.435, 1}, {0, 0.175, 0.345, 0.505, 1}, {0, 0.205, 0.415, 0.605, 1}, {0, 0.265, 0.505, 0.705, 1} };

int GetPtBin(float pt)
{
  int bin = -1;
  for(int ipt=0; ipt<nPtBins;++ipt){
    if(pt > ptbin[ipt] && pt < ptbin[ipt+1]){
      bin = ipt;
      break;
    }
  }
  return bin;
}

int GetAlphaBin(int _ptbin, float alpha)
{
  int bin = -1;
  for(int i=0; i<nAlphaBins;++i){
    if(alpha > alphabin[_ptbin][i] && alpha < alphabin[_ptbin][i+1]){
      bin = i;
      break;
    }
  }
  return bin;
}

static const std::vector<double> aE_grid = {0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20};
static const std::vector<double> bE_grid = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16};
static const std::vector<double> cE_grid = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
static const std::vector<double> aPos_grid = {0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008};
static const std::vector<double> bPos_grid = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006};

const int total_combinations = aE_grid.size() * bE_grid.size() * cE_grid.size() * aPos_grid.size() * bPos_grid.size();

inline void GetParameterSet(int index, double& aE, double& bE, double& cE, double& aPos, double& bPos)
{
    int nA = aE_grid.size();
    int nB = bE_grid.size();
    int nC = cE_grid.size();
    int nAp = aPos_grid.size();
    int nBp = bPos_grid.size();

    int total = nA * nB * nC * nAp * nBp;
    if (index < 0 || index >= total)
        throw std::runtime_error("Invalid parameter index!");

    int iA  =  index % nA; index /= nA;
    int iB  =  index % nB; index /= nB;
    int iC  =  index % nC; index /= nC;
    int iAp =  index % nAp; index /= nAp;
    int iBp =  index % nBp;

    aE   = aE_grid[iA];
    bE   = bE_grid[iB];
    cE   = cE_grid[iC];
    aPos = aPos_grid[iAp];
    bPos = bPos_grid[iBp];
}

void getrelativewidth(double a, double ae, double &b, double &be)
{
  double ratio = b/a;
  double sigmaratio = sqrt(be*be/(a*a) + ae*ae*b*b/(a*a));
  b = ratio;
  be = sigmaratio;
}


#endif
