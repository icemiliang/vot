// vot: Variational Optimal Transportation
//
// Author   : Mi,Liang (Arizona State University)
// Email    : icemiliang@gmail.com
// Date     : Sept 27th 2018

#ifndef OTX_HH
#define OTX_HH

#include "pointX.hh"

namespace vot {

class OTX {
public:
	OTX(){};
	~OTX();

	void setup(const int pMaxIterD, const int pMaxIterH, const double pThres,
			   const double pLearnRate, const bool pFlagVerb, const double vPlotScale,
			   const std::string pFilePrefix, const bool pFlagDebug);

	void import_data(const std::string pDiracFile, const std::string pEmpiricalFile);
	void cluster();

protected:
	void read_empiricalX_from_file(const std::string filename);
	void read_diracX_from_file(const std::string filename);
	void check_total_mass();

	void write_results(const int iterD, const int iterH);
	void write_dirac(const std::string filename); 
	void write_empirical(const std::string filename);

	bool update_map(int iterD, int iterH);
	void update_h();
	bool update_dirac();
	
	double compute_gradient();

	std::vector<DiracX*> vDiracXs;
	std::vector<EmpiricalX*> vEmpiricalXs;

	std::vector<double> vHs;
	std::vector<double> vGradient;

	double vLearnRate;
	double vThres;
	bool vFlagVerb, vFlagDebug;
	int vMaxIterD, vMaxIterH;
	int vNumDiracXs, vNumEmpiricalXs;

	int vDim;

	std::string vPrefix;
};

}

#endif
