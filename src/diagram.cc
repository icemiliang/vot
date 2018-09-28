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

    void Diagram::setup(const bool verbose) {
        mFlagVerbose = verbose;
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

    bool Diagram::update_dirac_beta() {
        std::vector<Dirac> newDiracs(mNumDiracs, Dirac(0, 0, 0, 0, 0, false));
        for (int iterE = 0; iterE < mNumEmpiricals; iterE++) {
            // Index multiplied by 2
            newDiracs[mEmpiricals[iterE].cellIndex()] += 
                        Dirac(mEmpiricals[iterE].x() * mEmpiricals[iterE].mass(),
                              mEmpiricals[iterE].y() * mEmpiricals[iterE].mass(),
                              mEmpiricals[iterE].z() * mEmpiricals[iterE].mass(),
                              0,
                              mEmpiricals[iterE].mass());
        }

        // Update existing centroids
        for (int i = 0; i < mNumDiracs; i++) {
            newDiracs[i].reset_dirac(mDiracs[i].dirac());
            newDiracs[i].set_fix(mDiracs[i].fix());
            std::cout << "Mass #" << i << mDiracs[i].mass() << std::endl;

            // If centroid fixed, then only update mass
            if (newDiracs[i].fix()) {
                newDiracs[i] = Dirac(mDiracs[i].x(),
                                               mDiracs[i].y(),
                                               mDiracs[i].z(),
                                               mDiracs[i].dirac(),
                                               newDiracs[i].mass());
            }
            else if (newDiracs[i].mass() == 0) {
                newDiracs[i] = mDiracs[i];
            }
            else {
                newDiracs[i] /= newDiracs[i].mass();
            }
        }

        mDiracs = newDiracs;
        std::cout << "Mass #" << 0 << mDiracs[0].mass() << std::endl;
        // Reset h to 1
        std::fill(mHs.begin(), mHs.end(), 10.0);
        mHs.push_back(10.0);
        return true;
    }

    double Diagram::interpolate() {
        double maxLength = otMinDouble;
        // double totalLength = 0.0;
        int maxIndex = -1;

        // Find the longest segment and insert a Dirac in the middle
        for (int iterD = 0; iterD < mNumDiracs-1; iterD++) {
            Point tmp = Point(mDiracs[iterD+1].x()-mDiracs[iterD].x(),
                              mDiracs[iterD+1].y()-mDiracs[iterD].y(),
                              mDiracs[iterD+1].z()-mDiracs[iterD].z());
            // No need to do L-1 norm, save time.
            double len = tmp.norm_squared();
            std::cout << "Segment length: " << len << std::endl;
            // totalLength += len;
            if (len > maxLength) {
                maxLength = len;
                maxIndex = iterD;
            }
        }

        ASSERT_VOT(maxIndex >= 0 && maxLength > otMinDouble,
                   "Something is wrong during interpolation." << " maxIndex: " << maxIndex
                   << ", maxLength: " << maxLength);
        // Insert a new Dirac
        Dirac tmpDirac = Dirac((mDiracs[maxIndex].x() + mDiracs[maxIndex+1].x())/2,
                               (mDiracs[maxIndex].y() + mDiracs[maxIndex+1].y())/2,
                               (mDiracs[maxIndex].z() + mDiracs[maxIndex+1].z())/2,
                               (mDiracs[maxIndex].dirac() + mDiracs[maxIndex+1].dirac())/3,
                               0,
                               false);
        mDiracs.insert(mDiracs.begin()+maxIndex+1, tmpDirac);

        // Reduce dirac values of adjacent Dirac
        mDiracs[maxIndex].reset_dirac(mDiracs[maxIndex].dirac()/3*2);
        mDiracs[maxIndex+2].reset_dirac(mDiracs[maxIndex+2].dirac()/3*2);
        mNumDiracs += 1;

        return maxLength;
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

        // Compute gradient and their norm
        double gradNorm = 0;
        mGradient = Eigen::VectorXd::Zero(mNumDiracs);
        for (int i = 0; i < mNumDiracs; i++) {
            // This is not the real gradients but it's equivalent to the real ones, why? Liang Mi
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
            assert(mHs[i] > 0);
        }
        // I haven't found a way to directly update h
        // The alternative is to clear the diagram and define a new one with new h
        mDiagram->clear();
        for (int i = 0; i < mNumDiracs; i++) {
            mDiagram->put(i, mDiracs[i].x(), mDiracs[i].y(), mDiracs[i].z(), sqrt(mHs[i]));
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
