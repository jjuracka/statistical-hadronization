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
    bool stable;
    std::string name;
    Int_t PDGcode;
    Int_t degeneracy;
    Int_t statistics;
    Double_t mass;
    Double_t yield;
    Double_t yieldErr;
    std::vector<std::pair<Double_t, std::vector<Int_t>>> decayChannels;

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

    void loadDecayChannels() {
        if (stable) return;
        std::ifstream file(Form("particles/%s_decay.txt", name.c_str()));
        if (!file.is_open()) {
            std::cout << "\tERROR: Unable to open file with decay information." << std::endl;
            return;
        }
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            Double_t branchingRatio;
            ss >> branchingRatio;
            std::vector<Int_t> daughters;
            Int_t daughter;
            while (ss >> daughter) {
                daughters.push_back(daughter);
            }
            decayChannels.push_back(std::make_pair(branchingRatio, daughters));
        }
    }

    Particle(bool stable, std::string name, Int_t PDGcode, Int_t degeneracy, Int_t statistics, Double_t mass) : stable(stable), name(name), PDGcode(PDGcode), degeneracy(degeneracy), statistics(statistics), mass(mass) {calcYield(); calcYieldErr(); loadDecayChannels();}
};

std::vector<std::pair<Int_t, Particle>> particles;

// load data from external file
void loadData() {
    std::ifstream file("particles/PartList_PPB2014_CBHN.txt");
    if (!file.is_open()) {
        std::cout << "ERROR: Unable to open file." << std::endl;
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines
        std::stringstream ss(line);
        bool stable;
        std::string name;
        Int_t PDGcode;
        Int_t degeneracy;
        Int_t statistics;
        Double_t mass;
        ss >> stable >> name >> PDGcode >> degeneracy >> statistics >> mass;
        Particle p(stable, name, PDGcode, degeneracy, statistics, mass);
        particles.push_back(std::make_pair(PDGcode, p));
        std:cout << Form("loaded %s (PDG %i): stable %i, mass %f, degeneracy %d, statistics %d", p.name.c_str(), p.PDGcode, p.stable, p.mass, p.degeneracy, p.statistics) << std::endl;
        std::cout << "\t" << Form("calculated yield: %f +- %f", p.yield, p.yieldErr) << std::endl;
        if (!p.stable) {
            for (const auto& decayChannel : p.decayChannels) {
                Double_t branchingRatio = decayChannel.first;
                const auto& daughters = decayChannel.second;
                std::cout << "\t" << Form("decay channel with BR %f:", branchingRatio/100);
                for (const auto& daughter : daughters) {
                    std::cout << " " << daughter;
                }
                std::cout << std::endl;
            }
        }
    }
    file.close();
}

void calcYields() {
    loadData(); // loads data and calculates primary yields

    // Adjust the yield for each particle based on decay channels
    // this breaks the code currently
    // for (auto& pair : particles) {
    //     auto& particle = pair.second;
    //     if (!particle.stable) {
    //         for (const auto& decayChannel : particle.decayChannels) {
    //             Double_t branchingRatio = decayChannel.first;
    //             const auto& daughters = decayChannel.second;
    //             for (const auto& daughter : daughters) {
    //                 particles[std::abs(daughter)].second.yield += particle.yield * branchingRatio/100/daughters.size();
    //             }
    //         }
    //     }
    // }
}