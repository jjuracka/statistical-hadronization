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
    Int_t stable;
    std::string name;
    Int_t PDGcode;
    Int_t degeneracy;
    Int_t statistics;
    Double_t mass;
    Double_t yield;
    Double_t yieldErr;

    Particle(Int_t stable, std::string name, Int_t PDGcode, Int_t degeneracy, Int_t statistics, Double_t mass) : stable(stable), name(name), PDGcode(PDGcode), degeneracy(degeneracy), statistics(statistics), mass(mass), yield(0.0), yieldErr(0.0) {}

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

std::vector<std::pair<std::string, Particle>> particles;

// load data from external file
void loadData() {
    std::ifstream file("PartList_PPB2021_CBHN.txt");
    if (!file.is_open()) {
        std::cout << "ERROR: Unable to open file." << std::endl;
        return;
    }
    std::string line;
    Int_t counter = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines
        std::stringstream ss(line);
        Int_t stable;
        std::string name;
        Int_t PDGcode;
        Int_t degeneracy;
        Int_t statistics;
        Double_t mass;
        ss >> stable >> name >> PDGcode >> degeneracy >> statistics >> mass;
        Particle p(stable, name, PDGcode, degeneracy, statistics, mass);
        p.calcYield();
        p.calcYieldErr();
        particles.push_back(std::make_pair(name, p));
        counter++;
        std::cout << Form("particle %d: %s, %d, %d, %d, %d, %f; Yield = %d +- %d", counter, name.c_str(), stable, PDGcode, degeneracy, statistics, mass, p.yield, p.yieldErr) << std::endl;
    }
    file.close();
}

void calcYields() {
    std::cout << "LOADING DATA" << std::endl;
    loadData();
    std::cout << "ALL DATA LOADED" << std::endl;

    // print results
    // std::cout << "RESULTS" << std::endl;
    // for (auto p : particles) {
    //     std::cout << "Particle: " << p.first << std::endl;
    //     std::cout << "Yield: " << p.second.yield << " +/- " << p.second.yieldErr << std::endl;
    // }
}