#include "Superflow/WeightComponents.h"

#include <sstream>      // std::ostringstream

using namespace std;

namespace sflow {
	std::string WeightComponents::str() const
	{
		ostringstream oss;
		oss << " susynt: " << susynt
			<< " lepSf: " << lepSf
			<< " btag: " << btag
			<< " trigger: " << trigger
			<< " qflip: " << qflip
			<< " fake: " << fake;
		return oss.str();
	}
};