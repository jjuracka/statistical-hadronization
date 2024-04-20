// ROOT macro to implement the statistical hadronization model for the determination of particle yields in LHC collisions
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TMath.h>

// global variables
const Double_t V = 5280.0, Verr = 410.0; // volume of the fireball in fm^3
const Double_t muB = 0.7, muBerr = 3.8; // baryon chemical potential in MeV
const Double_t T = 156.5, Terr = 1.5; // temperature in MeV

const Int_t nIter = 1e6; // number of iterations for the infinite sum approximation

struct Particle {
    std::string name;
    Int_t stable;
    Int_t PDGcode;
    Int_t degeneracy;
    Int_t statistics;
    Double_t mass;
    Double_t yield;
    Double_t yieldErr;

    void calcYield() {
        Double_t sum = 0.0;
        for (Int_t k = 1; k <= nIter; k++) {
            sum += degeneracy*V/(2*TMath::Power(TMath::Pi(), 2)) * TMath::Power(-1*statistics, k+1) * TMath::Power(mass, 2)*T/k * TMath::BesselK(2, k*mass/T) * TMath::Exp(k*muB/T);
        }
        yield = sum;
    }

    void calcYieldErr() {
        Double_t dNdV = 0.0, dNdmuB = 0.0, dNdT = 0.0;
        for (Int_t k = 1; k <= nIter; k++) {
            dNdV += degeneracy/(2*TMath::Power(TMath::Pi(), 2)) * TMath::Power(-1*statistics, k+1) * TMath::Power(mass, 2)*T/k * TMath::BesselK(2, k*mass/T) * TMath::Exp(k*muB/T);
            dNdmuB += degeneracy*V/(2*TMath::Power(TMath::Pi(), 2)) * TMath::Power(-1*statistics, k+1) * TMath::Power(mass, 2) * TMath::BesselK(2, k*mass/T) * TMath::Exp(k*muB/T);
            dNdT += V*degeneracy*TMath::Power(mass, 2)*TMath::Power(-1*statistics, k+1)*TMath::Exp(k*muB/T)*TMath::BesselK(2, k*mass/T)/(2*k*TMath::Power(TMath::Pi(), 2)) 
                    - V*degeneracy*TMath::Power(mass, 3)*TMath::Power(-1*statistics, k+1)*(-TMath::BesselK(1, k*mass/T)/2 - TMath::BesselK(3, k*mass/T)/2)*TMath::Exp(k*muB/T)/(2*T*TMath::Power(TMath::Pi(), 2)) 
                    - V*degeneracy*TMath::Power(mass, 2)*muB*TMath::Power(-1*statistics, k+1)*TMath::Exp(k*muB/T)*TMath::BesselK(2, k*mass/T)/(2*T*TMath::Power(TMath::Pi(), 2)); // this monstrosity was obtained by using symbolic differentiation in Python SymPy
        }
        yieldErr = TMath::Sqrt(TMath::Power(Verr*dNdV, 2) + TMath::Power(muBerr*dNdmuB, 2) + TMath::Power(Terr*dNdT, 2));
    }
};

std::vector<Particle> particles;

// load data from external file
void loadData() {
    std::ifstream file("PartList_PPB2021_CBHN.txt");
    if (!file.is_open()) {
        std::cout << "ERROR: Unable to open file." << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string col1, col2; // replace with actual column names
        std::getline(ss, col1, '\t');
        std::getline(ss, col2, '\t');
        // Push columns into corresponding vectors
        column1.push_back(col1);
        column2.push_back(col2);
        // Add more columns as needed
    }
    file.close();
}



void calcYields() {

}