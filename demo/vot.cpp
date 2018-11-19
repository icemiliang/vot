// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Nov 18th 2018

#include "ot.hh"
#include <sys/stat.h>

using namespace vot;
using namespace std;

string appName, dmFile, emFile, outPrefix;
int maxIterP = 10, maxIterH = 1000;
double thres = 0.00000001, learnRate = 0.1, plotscale = 1;
bool verb = false, debug = false;

void print_usage() {
    cerr << "Usage:\n\n"
         << appName.substr(appName.find_last_of("/")+1)
         << " -e emfile -d dmfile [options]\n\n"
         << " -d [ --dirac ]       Dirac measure file\n"
         << " -e [ --empirical ]   empirical measure file\n"
         << " -p [ --iterD ]       max iteration for updating Dirac positions, default is 10.\n"
         << " -h [ --iterH ]       max iteration for updating 'H', default is 1000.\n"
         << " -t [ --threshold ]   threshold for terminating the program, default is 1e-06\n"
         << " -r [ --rate ]        learning rate, default is 0.1\n"
         << " -o [ --outdir ]      output directory\n"
         << " -v [ --verbose ]     verbose (both console and files), default is false.\n"
         << " -b [ --debug ]       multiple-point check, default is false.\n"
         << "    [ --help ]        print usage\n"
         << "\n"
         << " To run a this demo:\n"
         << "    votx -e data/empirical.votx -d data/dirac.votx\n\n";
}

void get_args(int argc, char **argv) {   
    appName =argv[0];
    
    if (argc <= 1) {
        cerr << "---- Variational Wasserstein clustering ---- \n\n";
        print_usage();
        exit(SUCCESS);
    }
    int i = 1;
    while (i < argc) {
        if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--maxiterp")) {
            maxIterP = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--maxiterh")) {
            maxIterH = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--dmfile")) {
            dmFile = argv[++i];
        }
        else if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--emfile")) {
            emFile =argv[++i];
        }
        else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threshold")) {
            thres = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--rate")) {
            learnRate = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verb = true;
        }
        else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--debug")) {
            debug = true;
        }
        else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "-outdir")) {
            struct stat st;
            if(stat(argv[++i],&st) != 0) {
                if (-1 == mkdir(argv[i], 0775)) {
                    printf("Error creating directory!\n");
                    exit(1);
                }
            }
            outPrefix += string(argv[i]) + "/";
        }
        else {                                  // illegal syntax
            cerr << "Unrecognized option.\n";
            exit(ERROR_IN_COMMAND_LINE);
        }
        i++;
    }
    if (dmFile.empty() || emFile.empty()) {
        cerr << "-d and -e options must be specified\n";
        exit(ERROR_IN_COMMAND_LINE);
    }
}

int main(int argc, char* argv[]) {
    get_args(argc, argv);

    // Show paramters // 
    int slashD = dmFile.find_last_of("/\\");
    int slashE = emFile.find_last_of("/\\");
    outPrefix +=  "D"  + dmFile.substr(slashD + 1, dmFile.length() - 6 - slashD) 
                     + "_E" + emFile.substr(slashE + 1, emFile.length() - 6 - slashE);
    
    cerr << "\n Dirac measure file: " << dmFile
         << "\n Empirical measure file: " << emFile
         << "\n Max iteration for updating centroid positions: " << maxIterP
         << "\n Max iteration for updating 'h': " << maxIterH
         << "\n Learning rate: " << learnRate
         << "\n Threshold for terminating the program: " << thres
         << "\n Output prefix: " << outPrefix
         << "\n Verbose flag: " << verb
         << "\n Debug flag: " << debug << std::endl;

    // ----------------------------//
    // ------- Start program ------//
    // ----------------------------//

    OT *ot = new OT();
    std::cout << "\n";

    ot->import_data(dmFile, emFile);

    std::cout << "--> Setting up parameters... \n";
    ot->setup(maxIterP, maxIterH, thres, learnRate, verb, plotscale, outPrefix, debug);
    std::cout << std::endl;

    std::cout << "--> Running clustering... \n\n";
    ot->cluster();

    return SUCCESS;
}

