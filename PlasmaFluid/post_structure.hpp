#pragma once
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <fstream>
#include "domain_structure.hpp"
#include "variable_structure.hpp"
#include "config_structure.hpp"

#include "PFM.hpp"

using namespace std;
class CPost 
{
	public:

		CPost();
		void Init() ;
		void OutputFlow( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int, int ) ;
		void PlotDecomposition( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var ) ;
		void OutputAverageFlow( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int ) ;
};
