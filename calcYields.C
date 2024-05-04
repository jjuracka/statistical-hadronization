// ROOT macro to implement the statistical hadronization model for the determination of particle yields in LHC collisions

// includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <TFile.h>
#include <TMath.h>

// global variables
const Double_t V = 5280.0, Verr = 410.0; // volume of the fireball in fm^3
const Double_t muB = 0.7, muBerr = 3.8; // baryon chemical potential in MeV
const Double_t T = 156.5; // temperature in MeV, not assuming uncertainty
const Int_t nIter = 100; // number of iterations for the infinite sum approximation

// structure to store particle information loaded from file
struct Particle {
    // store relevant info
    bool stable;
    std::string name;
    Int_t PDGcode;
    Int_t degeneracy;
    Int_t statistics;
    Double_t mass;
    Double_t yield;
    Double_t yieldErr;
    std::vector<std::pair<Double_t, std::vector<Int_t>>> decayChannels;

    // calculate primary yield and uncertainty
    void calcYield() {
        Double_t sum = 0.0;
        for (Int_t k = 1; k <= nIter; k++) {
            sum += degeneracy*V/(2*TMath::Power(TMath::Pi(), 2)) * TMath::Power(-1*statistics, k+1) * TMath::Power(mass, 2)*T/k * TMath::BesselK(2, k*mass/T) * TMath::Exp(k*muB/T);
        }
        yield = sum;
    }

    void calcYieldErr() {
        Double_t dNdV = 0.0, dNdmuB = 0.0;
        for (Int_t k = 1; k <= nIter; k++) {
            dNdV += degeneracy/(2*TMath::Power(TMath::Pi(), 2)) * TMath::Power(-1*statistics, k+1) * TMath::Power(mass, 2)*T/k * TMath::BesselK(2, k*mass/T) * TMath::Exp(k*muB/T);
            dNdmuB += degeneracy*V/(2*TMath::Power(TMath::Pi(), 2)) * TMath::Power(-1*statistics, k+1) * TMath::Power(mass, 2) * TMath::BesselK(2, k*mass/T) * TMath::Exp(k*muB/T);
        }
        yieldErr = TMath::Sqrt(TMath::Power(Verr*dNdV, 2) + TMath::Power(muBerr*dNdmuB, 2));
    }

    // load decay channels from file
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
            while (ss >> daughter) daughters.push_back(daughter);
            decayChannels.push_back(std::make_pair(branchingRatio, daughters));
        }
    }

    // constructor
    Particle(bool stable, std::string name, Int_t PDGcode, Int_t degeneracy, Int_t statistics, Double_t mass) : stable(stable), name(name), PDGcode(PDGcode), degeneracy(degeneracy), statistics(statistics), mass(mass) {calcYield(); calcYieldErr(); loadDecayChannels();}
};

std::vector<std::pair<Int_t, Particle>> particles; // for storage of particles, with PDG code as key

// load data from external file
void loadData() {
    std::ifstream file("particles/PartList_PPB2014_CBHN.txt"); // input file
    if (!file.is_open()) {
        std::cout << "ERROR: Unable to open file." << std::endl;
        return;
    }
    // read file line by line and store particle information
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
        std::cout << "\t" << Form("calculated primary yield: %f +- %f", p.yield, p.yieldErr) << std::endl;
        // load decay channels if particle is unstable
        if (!p.stable) {
            for (const auto& decayChannel : p.decayChannels) {
                Double_t branchingRatio = decayChannel.first;
                const auto& daughters = decayChannel.second;
                std::cout << "\t" << Form("decay channel with BR %f:", branchingRatio/100);
                for (const auto& daughter : daughters) std::cout << " " << daughter;
                std::cout << std::endl;
            }
        }
    }
    file.close();
}

// find particle by PDG code
Particle* getParticleByPDG(int PDGcode) {
    for (auto& pair : particles) {
        if (pair.first == PDGcode) return &pair.second;
    }
    return nullptr;
}

// maps to store secondary yields and uncertainties
std::map<int, double> yields;
std::map<int, double> uncertainties;
// calculate secondary yields
void calculateSecondaryYields(const Particle& particle, double yield, double uncertainty) {
    // if the particle is stable, just add the yield and uncertainty to the map
    if (particle.stable) {
        yields[particle.PDGcode] += yield;
        uncertainties[particle.PDGcode] += uncertainty;
        return;
    }

    // if the particle is unstable, loop over all decay channels and calculate the secondary yields
    for (const auto& decayChannel : particle.decayChannels) {
        double branchingRatio = decayChannel.first;
        const auto& daughters = decayChannel.second;

        for (const auto& daughterPDG : daughters) {
            Particle* daughter = getParticleByPDG(daughterPDG);
            if (daughter == nullptr) continue; // skip this daughter if it's not found
            double daughterYield = yield * branchingRatio / (100 * daughters.size()); // calculate secondary yield...
            double daughterUncertainty = uncertainty * branchingRatio / (100 * daughters.size()); // ... and its uncertainty
            calculateSecondaryYields(*daughter, daughterYield, daughterUncertainty); // recursion
        }
    }

    yields[particle.PDGcode] += yield;
    uncertainties[particle.PDGcode] += uncertainty;
}

void calcYields() {
    loadData(); // loads data and calculates primary yields

    for (const auto& pair : particles) {
        const Particle& particle = pair.second;
        calculateSecondaryYields(particle, particle.yield, particle.yieldErr);
    }

    // now loop over all stable particles again and print the total yield
    std::cout << "\nTotal yields of stable particles:" << std::endl;
    std::vector<std::string> namesStable;
    std::vector<Double_t> yieldsStable;
    std::vector<Double_t> uncertaintiesStable;
    for (const auto& pair : particles) {
        const Particle& particle = pair.second;
        if (particle.stable) {
            std::cout << Form("%s (PDG %i): %f +- %f", particle.name.c_str(), particle.PDGcode, yields[particle.PDGcode], uncertainties[particle.PDGcode]) << std::endl;
            namesStable.push_back(particle.name);
            yieldsStable.push_back(yields[particle.PDGcode]);
            uncertaintiesStable.push_back(uncertainties[particle.PDGcode]);
        }
    }
    // create a plot with the results for yields of stable particles
    TCanvas* cYields = new TCanvas("cYields", "Yields of stable particles", 1600, 900);
    cYields->SetLogy();
    cYields->SetBottomMargin(0.2);
    TH1D* hYields = new TH1D("hYields", "yields of stable particles;;yield", namesStable.size(), 0, namesStable.size());
    for (Int_t i = 0; i < namesStable.size(); i++) {
        hYields->GetXaxis()->SetBinLabel(i+1, namesStable[i].c_str());
        hYields->SetBinContent(i+1, yieldsStable[i]);
        hYields->SetBinError(i+1, uncertaintiesStable[i]);
    }
    hYields->SetStats(0);
    hYields->Draw("E");
    cYields->SaveAs("yields.pdf");
}