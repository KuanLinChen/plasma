#pragma once
//#include "PFM.hpp"
//#include "config_structure.hpp"
//#include "domain_structure.hpp"
#include "variable_structure.hpp"

using namespace std;
class CPost 
{
	public:

		CPost();
		void Init() ;
		void OutputFlow( boost::shared_ptr<CDomain> &m, boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int, int ) ;
		void OutputAverageFlow( boost::shared_ptr<CConfig> &config, boost::shared_ptr<CVariable> &var, int ) ;
};
