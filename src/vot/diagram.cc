// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#include "diagram.hh"

namespace vot {
    
    Diagram::Diagram(const otBBox bBox, const int numBlocks) {
        mDiagram = new voro::container_poly(bBox.xMin, bBox.xMax,
                                             bBox.yMin, bBox.yMax,
                                             bBox.zMin, bBox.zMax,
                                             numBlocks, numBlocks, numBlocks,
                                             false, false, false, 8);
        mNumEmpiricals = 0;
        mNumDiracs = 0;
    }

    void Diagram::setup(const bool verbose, const bool debug) {
        mFlagVerbose = verbose;
        mFlagDebug = debug;
    }

    // x, y, z, dirac, mass, h, fix
    void Diagram::add_dirac(double x, double y, double z, double d, double m, double h, bool fix) {
        mDiagram->put(mNumDiracs, x, y, z, sqrt(h)); // Add a centroid/particle in container
        mDiracs.push_back(Dirac(x,y,z,d,m,fix));
        mHs.push_back(h);
        mNumDiracs++;
    }

    void Diagram::add_empirical(double x, double y, double z, double m, int dIndex) {
        mEmpiricals.push_back(Empirical(x,y,z,m,dIndex));
        mNumEmpiricals++;
    }

    void Diagram::draw_diagram(std::string pFilename) {
        std::cout << "--> Writing power diagram to " << pFilename << std::endl;
        mDiagram->draw_diagram_gnuplot(pFilename.c_str());
    }

    void Diagram::write_dirac_as_input(std::string pFilename) {
        std::cout << "--> Writing resulting Dirac measures to " << pFilename << std::endl;
        std::ofstream outP(pFilename);
        outP << "Dirac " << mNumDiracs << std::endl;
        for (int i = 0; i < mNumDiracs; i++) {
            // Index starts from 0
            outP << i << " "
                 << mDiracs[i].x() << " " << mDiracs[i].y() << " "
                 << mDiracs[i].z() << " " << mDiracs[i].mass() << " "
                 << "0 " << mHs[i] << std::endl;
        }
        outP.close();
    }

    void Diagram::write_dirac(std::string pFilename) {
        std::cout << "--> Writing resulting Diract measures to a " << pFilename << std::endl;
        std::ofstream outP(pFilename);
        for (int i = 0; i < mNumDiracs; i++) {
            // Index starts from 0
            outP << mDiracs[i].x() << " " << mDiracs[i].y() << " "
                 << mDiracs[i].z() << " "  << mDiracs[i].dirac() << " "
                 << mDiracs[i].mass() << " "  << mHs[i] << std::endl;
        }
        outP.close();
    }

    void Diagram::write_results(int const iterD, int const iterH, std::string pFilePrefix) {
        std::cout << "File prefix: " << pFilePrefix << std::endl;
        std::string diagramFileName = pFilePrefix + "_diagram_" + std::to_string(iterD) 
                                                  + "_" + std::to_string(iterH)+ ".gnu";
        std::string pFileName = pFilePrefix + "_" + std::to_string(iterD) 
                                            + "_" + std::to_string(iterH) + ".vot";
        std::string pFileNameInput = pFilePrefix + "_" + std::to_string(iterD) 
                                                 + "_" + std::to_string(iterH) + "_input.vot";

        // Reconstruct the diagram in case the data stored in the voro++ object is not updated.
        mDiagram->clear();
        for (int i = 0; i < mNumDiracs; i++) {
            mDiagram->put(i, mDiracs[i].x(), mDiracs[i].y(), mDiracs[i].z(), sqrt(mHs[i]));
        }
        
        draw_diagram(diagramFileName);
        write_dirac(pFileName);
        write_dirac_as_input(pFileNameInput);
    }

    double Diagram::total_dirac_mass() {
        double total = 0;
        for (int i = 0; i < mNumDiracs; i++) {
            total += mDiracs[i].dirac();
        }
        return total;
    }

    double Diagram::total_empirical_mass() {
        double total = 0;
        for (int i = 0; i < mNumEmpiricals; i++) {
            total += mEmpiricals[i].mass();
        }
        return total;
    }

    // method: gradient descent or newton
    // Return true if gradient doesn't change too much.
    bool Diagram::update(int method, double thres, double step, int iterP, int iterH) {
        // Hessian is computed during updating the diagram
        Eigen::SparseMatrix<double> hessian(mNumDiracs, mNumDiracs);
        if (method == METHOD_NEWTON)
            mDiagram->update_diagram(mDiracs, hessian);
        else 
            mDiagram->update_diagram();

        // Return true if gradient doesn't change too much.
        if (mFlagVerbose) 
            std::cout << "iterP: " << iterP << ", iterH: " << iterH << ". ";
        
        if (compute_gradient() < thres) {
            return true;
        }

        if (method == METHOD_NEWTON)
            update_h(hessian);
        else
            update_h(step);

        return false;
    }

    bool Diagram::update_bf(int method, double thres, double step, int iterP, int iterH) {
        // Hessian is computed during updating the diagram

        // Return true if gradient doesn't change too much.
        if (mFlagVerbose) 
            std::cout << "iterP: " << iterP << ", iterH: " << iterH << ". ";
        
        if (compute_gradient_bf() < thres) {
            return true;
        }

        update_h(step);

        return false;
    }

    bool Diagram::update_dirac() {
        std::vector<Dirac> newDiracs(mNumDiracs, Dirac(0, 0, 0, 0, 0, false));
        // Sum up all empirical measures
        for (int iterE = 0; iterE < mNumEmpiricals; iterE++) {
            newDiracs[mEmpiricals[iterE].cellIndex()] += 
                                    Dirac(mEmpiricals[iterE].x() * mEmpiricals[iterE].mass(),
                                          mEmpiricals[iterE].y() * mEmpiricals[iterE].mass(),
                                          mEmpiricals[iterE].z() * mEmpiricals[iterE].mass(),
                                          0,
                                          mEmpiricals[iterE].mass(),
                                          false);
        }

        // Compute new Dirac measures
        for (int iterD = 0; iterD < mNumDiracs; iterD++) {
            // Use old dirac measure
            newDiracs[iterD].reset_dirac(mDiracs[iterD].dirac());

            // Divide x,y,z by mass
            double mass = newDiracs[iterD].mass();
            // If cell is empty, do not update the dirac
            if (mass < vot::otTOLERANCE) {
                newDiracs[iterD] = mDiracs[iterD];
                continue;
            } 
            newDiracs[iterD] /= mass;
            newDiracs[iterD].set_mass(mass);
        }

        // Check if old and new Dirac measures have similar locations
        for (int iterD = 0; iterD < mNumDiracs; iterD++) {
            if (fabs(mDiracs[iterD].x() - newDiracs[iterD].x()) > otTOLERANCE ||
                fabs(mDiracs[iterD].y() - newDiracs[iterD].y()) > otTOLERANCE ||
                fabs(mDiracs[iterD].z() - newDiracs[iterD].z()) > otTOLERANCE) {
                // Update Dirac measures
                mDiracs = newDiracs;
                
                // Reset h to 1
                std::fill(mHs.begin(), mHs.end(), 1.0);
                return false;
            }
        }

        // Update Dirac measures
        mDiracs = newDiracs;
        return true;
    }

    double Diagram::compute_gradient() {
        // Reset cell mass
        for (int i = 0; i < mNumDiracs; i++) {
            mDiracs[i].reset_mass();
        }
        
        // Assign cellIndex to every measure and add mass to that cell
        int cellIndex = -1;
        for (int i = 0; i < mNumEmpiricals; i++) {
            cellIndex = mDiagram->cell_index(mEmpiricals[i].x(),
                                             mEmpiricals[i].y(), 
                                             mEmpiricals[i].z());
            mEmpiricals[i].set_cellIndex(cellIndex);
            mDiracs[cellIndex].add_mass(mEmpiricals[i].mass());
        }

        // Multiple-point check
        if (mFlagDebug) {
            double sumMass = 0;
            double sumH = 0;
            for (int i = 0; i < mNumDiracs; i++) {
                sumMass += mDiracs[i].mass();
                sumH += mHs[i];
            }
            std::cout << "[Debug] Total mass of Diracs: " << sumMass << std::endl;
            std::cout << "[Debug] Total h of Diracs: " << sumH << std::endl;
            getchar();
        }


        // Compute gradient and their norm
        double gradNorm = 0;
        mGradient = Eigen::VectorXd::Zero(mNumDiracs);
        for (int i = 0; i < mNumDiracs; i++) {
            mGradient(i) = mDiracs[i].mass() - mDiracs[i].dirac();
        }
        gradNorm = mGradient.norm();
        if (mFlagVerbose) {
            std::cout << "Gradient norm: " << gradNorm << std::endl;
        }

        return gradNorm;
    }

    // brute force
    double Diagram::compute_gradient_bf() {
        // Reset cell mass
        for (int i = 0; i < mNumDiracs; i++) {
            mDiracs[i].reset_mass();
        }

        // Assign cellIndex to every measure and add mass to that cell
        for (int i = 0; i < mNumEmpiricals; i++) {

            // Find assignment by bruteforce
            double dist = otMaxDouble;
            int idx = -1;
            for (int j = 0; j < mNumDiracs; j++) {
                double tmpDist1 = mEmpiricals[i].costL2_tmp(mDiracs[j]);
                double tmpDist2 = tmpDist1 - mHs[j];
                if (tmpDist2 < dist) {
                    dist = tmpDist2;
                    idx = j;
                }
            }

            mEmpiricals[i].set_cellIndex(idx);
            mDiracs[idx].add_mass(mEmpiricals[i].mass());            
        }

        // Multiple-point check
        if (mFlagDebug) {
            double sumMass = 0;
            double sumH = 0;
            for (int i = 0; i < mNumDiracs; i++) {
                sumMass += mDiracs[i].mass();
                sumH += mHs[i];
            }
            std::cout << "[Debug] Total mass of Diracs: " << sumMass << std::endl;
            std::cout << "[Debug] Total h of Diracs: " << sumH << std::endl;
            getchar();
        }

        // Compute gradient and its norm
        double gradNorm = 0;
        mGradient = Eigen::VectorXd::Zero(mNumDiracs);
        for (int i = 0; i < mNumDiracs; i++) {
            mGradient(i) = mDiracs[i].mass() - mDiracs[i].dirac();
        }
        gradNorm = mGradient.norm();
        if (mFlagVerbose) {
            std::cout << "Gradient norm: " << gradNorm << std::endl;
        }

        return gradNorm;
    }

    // method: gradient descent or newton
    void Diagram::update_h(double step) {
        for (int i = 0; i < mNumDiracs; i++) {
            mHs[i] -= step * mGradient[i]; // delta h
            ASSERT_VOT(mHs[i] > 0, "Negative h not allowed");
        }
        // I haven't found a way to directly update h
        // The alternative is to clear the diagram and define a new one with new h
        mDiagram->clear();
        for (int i = 0; i < mNumDiracs; i++) {
            mDiagram->put(i, mDiracs[i].x(), mDiracs[i].y(), mDiracs[i].z(), sqrt(mHs[i]));
        }
    }

    // method: gradient descent or newton
    void Diagram::update_h_bf(double step) {
        for (int i = 0; i < mNumDiracs; i++) {
            mHs[i] -= step * mGradient[i]; // delta h
            if(mHs[i] <= 0) {
                std::cout << "[warning] H: " << mHs[i] << std::endl;
            }
        }
    }

    void Diagram::update_h(Eigen::SparseMatrix<double>& pHessian) {
        Eigen::VectorXd deltaH;
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
        deltaH = solver.compute(pHessian).solve(mGradient);

        for (int i = 0; i < mNumDiracs; i++) {
            mHs[i] -= deltaH(i);
        }

        // I haven't found a way to directly update h
        // The alternative is to clear the diagram and define a new one with new h
        // Liang Mi
        mDiagram->clear();
        for (int i = 0; i < mNumDiracs; i++) {
            mDiagram->put(i, mDiracs[i].x(), mDiracs[i].y(), mDiracs[i].z(), sqrt(mHs[i]));
        }
    }

    double Diagram::compute_wasserstein() {
        double wd = 0;
        for (int i = 0; i < mNumEmpiricals; i++) {
            wd += mEmpiricals[i].costL2(mDiracs[mEmpiricals[i].cellIndex()]);
        }
        return wd;
    }

}
