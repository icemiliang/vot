// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Oct 16th 2018

#include "otx.hh"

namespace votx {
    void OTX::setup(const int pMaxIterD, const int pMaxIterH, const double pThres,
        const double pLearnRate, const bool pFlagVerb, const double vPlotScale,
        const std::string pPrefix, const bool pFlagDebug) {
        vFlagVerb = pFlagVerb;
        vFlagDebug = pFlagDebug;
        vPrefix = pPrefix;
        vMaxIterD = pMaxIterD;
        vMaxIterH = pMaxIterH;
        vThres = pThres;
        vLearnRate = pLearnRate;

        vNumDiracXs = vDiracXs.size();
        vNumEmpiricalXs = vEmpiricalXs.size();

        // Set h to 1 and gradient to 0
        vHs = std::vector<double>(vNumDiracXs, 1.0);
        vGradient = std::vector<double>(vNumDiracXs, 0.0);
    }

    void OTX::import_data(const std::string pDiracFile, const std::string pEmpiricalFile) {
        std::cout << "--> Importing data... " << std::endl;
        read_diracX_from_file(pDiracFile);
        read_empiricalX_from_file(pEmpiricalFile);
        std::cout << std::endl;
        check_total_mass();
    }

    void OTX::read_diracX_from_file(std::string pFilename) {
        std::string line, tmp, sNum, sDim, sDirac;
        std::ifstream file(pFilename);
        do {getline(file, line);} while (line[0] == '#');
        std::stringstream pline(line);
        pline >> tmp >> sNum >> sDim;
        int num = stoi(sNum);
        int dim = stoi(sDim);
        vDim = dim;

        // Read samples line by line
        for (int i = 0; i < num; i++) {
            do {getline(file, line);} while (line[0] == '#');
            std::stringstream pline(line);
            pline >> sDirac;
            double* coors = new double[dim];
            for (int j = 0; j < dim; j++) {
                ASSERT_VOT(pline >> tmp, "# coor elements != dimension");
                coors[j] = stod(tmp);
            }
            DiracX *tmpdiracx = new DiracX(dim, *coors, stod(sDirac));
            vDiracXs.push_back(tmpdiracx);
        }
        vNumDiracXs = vDiracXs.size();
    }

    void OTX::read_empiricalX_from_file(std::string pFilename) {
        std::string line, tmp, sNum, sDim, sMass;
        std::ifstream file(pFilename);
        do {getline(file, line);} while (line[0] == '#'); 
        std::stringstream pLine(line);
        pLine >> tmp >> sNum >> sDim;
        int num = stoi(sNum);
        int dim = stoi(sDim);

        // Read samples line by line
        for (int i = 0; i < num; i++) {
            do {getline(file, line);} while (line[0] == '#');
            std::stringstream pline(line);
            pline >> sMass;
            double* coors = new double[dim];
            for (int j = 0; j < dim; j++) {
                ASSERT_VOT(pline >> tmp, "# coor elements != dimension");
                coors[j] = stod(tmp);
            }
            EmpiricalX *tmpempiricalx = new EmpiricalX(dim, *coors, stod(sMass));
            vEmpiricalXs.push_back(tmpempiricalx);
        }
        vNumEmpiricalXs = vEmpiricalXs.size();
    }

    // Total mass of empiricals = total dirac of diracs
    void OTX::check_total_mass() {
        double totaldirac = 0, totalmass = 0;
        for (auto d : vDiracXs) totaldirac += d->dirac();
        for (auto e : vEmpiricalXs) totalmass += e->mass();

        if (vFlagVerb) {
            std::cout << "--> Checking total measures... " << std::endl;
            std::cout << "       Total Dirac measure: " << totaldirac << std::endl;
            std::cout << "       Total Empirical measure: " << totalmass << std::endl;
            std::cout << "       Diff: " << fabs(totaldirac - totalmass) << std::endl;
            std::cout << std::endl;
        }
        ASSERT_VOT(fabs(totaldirac - totalmass) / totaldirac < OTX_TOLERANCE*vNumEmpiricalXs,
               "Total Dirac measure: " << totaldirac << 
               " not equal to total Empirical measure: " << totalmass);
    }

    // Main function
    void OTX::cluster() {
        int iterD = 0;
        while (iterD < vMaxIterD) {
            int iterH = 0;
            while (iterH < vMaxIterH) {
                // break if map converges 
                if (update_map(iterD, iterH)) {
                    if (vFlagVerb) write_results(iterD,iterH);
                    break;
                } 
                if (vFlagVerb) write_results(iterD,iterH);
                iterH++;
            }

            // break if dirac converges
            if (update_dirac()) { write_results(iterD,iterH); break; }
            write_results(iterD,iterH);

            if (vFlagDebug) getchar();
            
            iterD++;
        }
    }

    bool OTX::update_map(int iterD, int iterH) {
        if (vFlagVerb || vFlagDebug) {
            std::cout << "updating map... iterD: " << iterD 
                                    << ", iterH: " << iterH << std::endl;    
        }
        // If gradient norm less than threshold
        if (compute_gradient() < vThres) { return true; }
        update_h();
        return false;
    }

    bool OTX::update_dirac() {
        if (vFlagVerb || vFlagDebug) { std::cout << "updating dirac... " << std::endl; }
        int converge = true;
        // Create a new double array to replace the old one.
        for (auto d : vDiracXs) {
            // If index set does not change, no need to update.
            if (!d->update_flag()) { if (vFlagDebug) d->print(); continue; }

            double* coor = new double[d->dim()];
            std::fill(&coor[0], &coor[d->dim()], 0.0);
            std::unordered_set<int>& tmpSet = d->set();
            for (auto idx : tmpSet) for (int i = 0; i < vDim; i++) 
                coor[i] += (*vEmpiricalXs[idx])[i]*vEmpiricalXs[idx]->mass();

            d->update_coor(*coor);
            *d /= d->mass();

            if (vFlagDebug) d->print();
            converge = false;
        }
        // 0 means no dirac points have been updated.
        return converge;
    }

    double OTX::compute_gradient() {
        // Reset mass
        for (auto d : vDiracXs) { d->set_mass(0.0); d->clear_idx(); d->reset_upate_flag(); }

        // Loop over empiricals, find nearest dirac by brute-force add mass
        for (int i = 0; i < vNumEmpiricalXs; i++) {
            double dist = OTX_MAX_DOUBLE;
            int newIdx = -1;

            // Brute-force search for nearest dirac
            for (int j = 0; j < vNumDiracXs; j++) {
                double tmpDist = vEmpiricalXs[i]->distL22(*vDiracXs[j]) - vHs[j];
                if (tmpDist < dist) {
                    dist = tmpDist;
                    newIdx = j;
                }
            }
            ASSERT_VOT(newIdx >= 0, "index < 0, no nearby dirac found for a empirical.");
            // Assign index 

            int oldIdx = vEmpiricalXs[i]->di();
            if (oldIdx != newIdx) {
                vEmpiricalXs[i]->set_dirac_idx(newIdx);
                vDiracXs[newIdx]->set_to_update();
                if (oldIdx != -1) vDiracXs[oldIdx]->set_to_update();
            }
            vDiracXs[newIdx]->add_empirical_idx(i);
            vDiracXs[newIdx]->add_mass(vEmpiricalXs[i]->mass());
        }

        if (vFlagDebug) {
            for (auto e : vEmpiricalXs) e->print();
            for (auto d : vDiracXs) {
                d->print();
                std::unordered_set<int> tmpSet = d->set();
                std::cout << "enclosed indices: ";
                for (auto j : tmpSet) std::cout << j << ", ";
                std::cout << std::endl;
            }
        }

        // update gradient 
        double maxDer = 0.0;
        for (int i = 0; i < vNumDiracXs; i++) {
            vGradient[i] = vDiracXs[i]->mass() - vDiracXs[i]->dirac();
            if (vGradient[i] > maxDer) maxDer = vGradient[i];
        }

        if (vFlagVerb || vFlagDebug) {
            std::cout << "Largest derivative: " << maxDer << std::endl;
        }
        
        return maxDer;
    }

    void OTX::update_h() {
        for (int i = 0; i < vNumDiracXs; i++) {
            vHs[i] -= vLearnRate * vGradient[i];
            if((vFlagVerb || vFlagDebug) && vHs[i] < 0) {
                std::cout << "[warning] H: " << vHs[i] 
                          << ", voro++ cannot plot if h is negative. " << std::endl;
            }
        }
    }

    void OTX::write_results(const int iterD, const int iterH) {
        std::cout << "writing results...\n    File prefix: " << vPrefix << std::endl;

        std::string dFilename = vPrefix + "_" + std::to_string(iterD) 
                                        + "_" + std::to_string(iterH) + ".vot";
        std::string iFilename = vPrefix + "_" + std::to_string(iterD) 
                                        + "_" + std::to_string(iterH) + ".index";        
        write_dirac(dFilename);
        write_empirical(iFilename);
    }

    void OTX::write_dirac(std::string filename) {
        std::ofstream of(filename);
        of << "dirac " << vNumDiracXs << " " << vDim << std::endl;
        for (int i = 0; i < vNumDiracXs; i++) {
            of << "{ mass: " << vDiracXs[i]->mass() << " dirac: " << vDiracXs[i]->dirac() 
               << " h: " << vHs[i] << " } ";
            for (int j = 0; j < vDim; j++) {
                of << (*vDiracXs[i])[j] << " ";
            }
            of << std::endl;
        }
        of.close();
    }

    void OTX::write_empirical(std::string filename) {
        std::ofstream of(filename);
        of << "empirical " << vNumEmpiricalXs << " " << vDim << std::endl;
        for (int i = 0; i < vNumEmpiricalXs; i++) {
            of << "{ mass: " << vEmpiricalXs[i]->mass() << " dirac index: " << vEmpiricalXs[i]->di() << " } ";
            for (int j = 0; j < vDim; j++) {
                of << (*vEmpiricalXs[i])[j] << " ";
            }
            of << std::endl;
        }
        of.close();
    }
}
