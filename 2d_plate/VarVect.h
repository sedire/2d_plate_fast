#ifndef _PLATE_2D_VARVECT_
#define _PLATE_2D_VARVECT_ 1

#include <vector>
#include "plate_var_types.h"
using std::vector;

class VarVect
{
public:
	VarVect();
	VarVect( int _varNum );
	~VarVect();
	void setup( int _varNum );

	vector<PL_NUM> Nk;
	vector<PL_NUM> Nk1;
	vector<PL_NUM> d1N;
	vector<PL_NUM> d2N;

	vector<PL_NUM> Nk0;			//don't really know why we need these. some computational tricks, probably
	vector<PL_NUM> d1N0;
	vector<PL_NUM> d2N0;
};

#endif