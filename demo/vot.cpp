// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#include "ot.hh"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

using namespace vot;
namespace po = boost::program_options;

void print_welcome_message() {
    std::cout << "Variational Wasserstein clustering" << std::endl << std::endl;
}

int main(int argc, char* argv[]) {

    try {
        std::string appName = boost::filesystem::basename(argv[0]); 
        std::string dmFile, emFile;
        int maxIterP, maxIterH; // Position of centroids & Minimizier H
        double thres, plotScale, learnRate;
        std::string outDir;
        bool verb;

        // Define and parse the program options //
        po::options_description desc("Options");
        desc.add_options()
            ("dirac,d", po::value<std::string>(&dmFile)->required(), 
                "Particles associated with Dirac measures")
            ("empirical,e", po::value<std::string>(&emFile)->required(), 
                "Test points associated with Empirical measures")
            ("iterD,p", po::value<int>(&maxIterP)->default_value(5), 
                "Max iteration for updating Dirac positions. Default is 5.")
            ("iterH,h", po::value<int>(&maxIterH)->default_value(1000), 
                "Max iteration for updating diagram/'H', or size of cells. Default is 1000.")
            ("threshold,t", po::value<double>(&thres)->default_value(1e-06), 
                "Threshold for terminating the program. Default is 1e-06")
            ("scale,s", po::value<double>(&plotScale)->default_value(1), 
                "Scale of the output diagram. Default is 1.")
            ("rate,r", po::value<double>(&learnRate)->default_value(0.2), 
                "Learning rate. Default is 0.2.")
            ("outdir,o", po::value<std::string>(&outDir), 
                "Output directory")
            ("verbose,v", po::value<bool>(&verb)->default_value(false), 
                "Verbose output (both console and files). Default is False.")
            ("help", "Print help messages");
        po::variables_map vm;
        
        try {
            po::store(po::parse_command_line(argc, argv, desc), vm); // can throw 

            //  --help option  //
            if (vm.count("help") || argc < 2) {
                print_welcome_message();
                std::cout << desc << std::endl;
                return SUCCESS;
            }

            po::notify(vm); // TODO: What does this line do?
        }
        catch(po::required_option& e) {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
            return ERROR_IN_COMMAND_LINE; 
        }
        catch(po::error& e){
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        // Displaying paramters // 
        int slashD = dmFile.find_last_of("/\\");
        int slashE = emFile.find_last_of("/\\");
        std::string outPrefix = "D" 
                              + dmFile.substr(slashD + 1, dmFile.length() - 5 - slashD) 
                              + "_E" 
                              + emFile.substr(slashE + 1, emFile.length() - 5 - slashE);
                              
        if (!outDir.empty()) outPrefix = outDir + "/" + outPrefix;

        std::cout << "Dirac measure file: " << dmFile << std::endl;
        std::cout << "Empirical measure file: " << emFile << std::endl;
        std::cout << "Max iteration for updating centroid positions: " << maxIterP << std::endl;
        std::cout << "Max iteration for updating 'H': " << maxIterH << std::endl;
        std::cout << "Learning rate: " << learnRate << std::endl;
        std::cout << "Threshold for terminating the program: " << thres << std::endl;
        std::cout << "Scale of the output diagram: " << plotScale << std::endl;
        std::cout << "Output prefix: " << outPrefix << std::endl;

        // ----------------------------//
        // ------- Start program ------//
        // ----------------------------//

        OT *ot = new OT();
        std::cout << std::endl;

        ot->import_data(dmFile, emFile);

        std::cout << "--> Setting up parameters... " << std::endl;
        ot->setup(maxIterP, maxIterH, thres, learnRate, verb, plotScale, outPrefix);
        std::cout << std::endl;

        std::cout << "--> Running variational Wasserstein clustering..." << std::endl << std::endl;
        ot->cluster();

    }
    catch(std::exception& e){
        std::cerr << "Unhandled Exception reached the top of main: " 
                  << e.what() << ", application will now exit" << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;
    }

    return SUCCESS;
}
