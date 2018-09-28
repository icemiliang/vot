// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#include "ot.hh"

namespace vot {

    void OT::import_data(std::string pDiracFile, std::string pEmpiricalFile) {
        std::cout << "--> Importing data... " << std::endl;
        read_dirac_from_file(pDiracFile);
        read_empirical_from_file(pEmpiricalFile);
        std::cout << std::endl;
        check_total_mass();
    }

    void OT::setup(const int pMaxIterD, const int pMaxIterH, const double pThres, 
                   const double pLearnRate, const bool pFlagVerbose, const double pPlotScale,
                   const std::string pFilePrefix) {
        mFlagVerbose = pFlagVerbose;
        mOutFilePrefix = pFilePrefix;
        mMaxIterD = pMaxIterD;
        mMaxIterH = pMaxIterH;
        mThres = pThres;
        mLearnRate = pLearnRate;

        mDiagram->setup(mFlagVerbose);
    }

    void OT::read_dirac_from_file(std::string pFilename) {
        std::string line, temp, tempNumP;
        std::ifstream measureFile(pFilename);
        getline(measureFile, line);
        std::stringstream pLine(line);
        pLine >> temp >> tempNumP;
        int numP = stoi(tempNumP);
        // x, y, z, Dirac, mass, h, fix
        std::string x, y, z, d, m, h, f;
        for (int i = 0; i < numP; i++) {
            getline(measureFile, line);
            std::stringstream pLine(line);
            pLine >> x >> y >> z >> d >> m >> h >> f;
            // Dirac measure must be non-negative
            assert(stod(d) >= 0.0);
            
            // Add the dirac to the mDiagram
            if (m.empty())
                mDiagram->add_dirac(stod(x), stod(y), stod(z), stod(d), 0.0, 1.0, 0.0); 
            else if (h.empty())
                mDiagram->add_dirac(stod(x), stod(y), stod(z), stod(d), stod(m), 1.0, 0.0);   
            else if (f.empty()) 
                mDiagram->add_dirac(stod(x), stod(y), stod(z), stod(d), stod(m), stod(h), 0.0);
            else
                mDiagram->add_dirac(stod(x), stod(y), stod(z), stod(d), stod(m), stod(h), stod(f));    
        }
    }

    void OT::read_empirical_from_file(std::string pFilename) {
        std::string line, temp, tempNumM;
        std::ifstream measureFile(pFilename);
        getline(measureFile, line);
        std::stringstream pLine(line);
        pLine >> temp >> tempNumM;
        int numM = stoi(tempNumM);
        // x, y, z, mass, diracIndex
        std::string x, y, z, m, dIdx;
        for (int i = 0; i < numM; i++) {
            getline(measureFile, line);
            std::stringstream pLine(line);
            pLine >> x >> y >> z >> m >> dIdx;
            if (dIdx.empty()) dIdx = "-1";
            // All measures are non-negative
            assert(stod(m) >= 0.0);

            // Add an empirical measure to mDiagram
            mDiagram->add_empirical(stod(x), stod(y), stod(z), stod(m), stoi(dIdx));
        }
    }

    void OT::check_total_mass() {
        double totalDirac = mDiagram->total_dirac_mass();
        double totalEmpirical = mDiagram->total_empirical_mass();
        if (mFlagVerbose) {
            std::cout << "--> Checking total measures... " << std::endl;
            std::cout << "       Total Dirac measure: " << totalDirac << std::endl;
            std::cout << "       Total Empirical measure: " << totalEmpirical << std::endl;
            std::cout << "       Diff: " << fabs(totalDirac - totalEmpirical) << std::endl;
            std::cout << std::endl;
        }
        ASSERT_VOT(fabs(totalDirac - totalEmpirical) / totalDirac < otTOLERANCE,
               "Total Dirac: " << totalDirac << " not equal to total Empirical: " << totalEmpirical);
    }

    void OT::cluster() {
        int iterD = 0;
        mDiagram->write_results(0, 0, mOutFilePrefix);
        while (iterD < mMaxIterD) {
            // OT 
            int iterH = 0;
            while (iterH < mMaxIterH) {
                if (mDiagram->update(METHOD_NEWTON, mThres, mLearnRate, iterD, iterH)) {
                    break;
                }
                iterH++;
                // mDiagram->write_results(iterD, iterH, mOutFilePrefix);
            }
            // Update dirac
            if (mDiagram->update_dirac()) {
                break;
            }
            iterD++;
            mDiagram->write_results(iterD, 0, mOutFilePrefix);
        }
        // mDiagram->write_results(iterD, iterH, mOutFilePrefix);
    }

}
