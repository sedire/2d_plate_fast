#include "OrthoBuilder.h"

SolInfo::SolInfo()
{

}

void SolInfo::setup( int _varNum )
{
	o.resize( ( _varNum / 2 + 1 ) * ( _varNum / 2 + 2 ) / 2, 0.0 );		//the matrix is diagonal. such size is to store less  TODO: maybe I can make the size smaller by 1
	//zi.resize( _varNum / 2, vector< PL_NUM>( _varNum, 0.0 ) );
	//z5.resize( _varNum, 0.0 );

	C.resize( _varNum / 2, 0.0 );				//FIXME may be another size here
}

void SolInfo::flushO()
{
	for( int i = 0; i < o.size(); ++i )
	{
		o[i] = 0.0;
	}
}

SolInfo::~SolInfo()
{

}

OrthoBuilder::OrthoBuilder()
{
}

OrthoBuilder::~OrthoBuilder()
{

}

OrthoBuilderGSh::OrthoBuilderGSh( int _varNum)
{
	eq_num = _varNum;
	LL.resize( eq_num / 2, vector<PL_NUM>( eq_num / 2, 0.0 ) );
	UU.resize( eq_num / 2, vector<PL_NUM>( eq_num / 2, 0.0 ) );
}

OrthoBuilderGodunov::OrthoBuilderGodunov( int _varNum)
{
	eq_num = _varNum;
}

void OrthoBuilder::flushO( int x )
{
	solInfoMap[x].flushO();
}

void OrthoBuilder::setParams( int _Km )
{
	Km = _Km;

	try 
	{
		solInfoMap.resize( Km );
		for( int i = 0; i < solInfoMap.size(); ++i )
		{
			solInfoMap[i].setup( eq_num );
		}
	}
	catch( bad_alloc &ba )
	{
		cout << ba.what() << endl;
	}

	for( int i = 0; i < NODES_ON_Y; ++i )
	{
		for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
		{
			z5[i][j] = 0.0;
		}
	}
}

void OrthoBuilder::setInitVects( const vector<PL_NUM>& N1, const vector<PL_NUM>& N2, const vector<PL_NUM>& N3, const vector<PL_NUM>& N4, const vector<PL_NUM>& N5 )
{
	//for( int i = 0; i < eq_num; ++i )
	//{
	//	solInfoMap[0].z1[i] = N1[i];
	//	solInfoMap[0].z2[i] = N2[i];
	//	solInfoMap[0].z3[i] = N3[i];
	//	solInfoMap[0].z4[i] = N4[i];
	//	solInfoMap[0].z5[i] = N5[i];
	//}
}

void OrthoBuilder::LUsolve( vector<vector<PL_NUM>>& AA, vector<PL_NUM>& ff, vector<PL_NUM>* xx )
{
	if( AA.size() < 1 )
	{
		cout << "ERROR in LUsolve: AA.size() is less than 1\n";
		return;
	}
	int AAsize = AA.size();

	vector<PL_NUM> storage( AAsize, 0.0 );
	PL_NUM tmpNum = 0;

	int nzRow = 0;
	if( fabs( AA[0][0] ) < ALMOST_ZERO )		//TODO may be I can join this and the next block. I mean the checks on AA[i][j] != 0
	{
		//make AA[0][0] nonzero by changing the order of rows 
		int nzFound = false;
		for( nzRow; nzRow < AAsize; ++nzRow )
		{
			if( fabs( AA[nzRow][0] ) > ALMOST_ZERO )			//CHECK, may be it is better to put != 0 here
			{
				nzFound = true;
				break;
			}
		}
		if( nzFound == false )
		{
			cout << "ERROR in LUsolve: no nonzero elements found\n";
			return;
		}

		for( int col = 0; col < AAsize; ++col )
		{
			storage[col] = AA[nzRow][col];
			AA[nzRow][col] = AA[0][col];
			AA[0][col] = storage[col];
		}
		tmpNum = ff[nzRow];
		ff[nzRow] = ff[0];
		ff[0] = tmpNum;
	}
	for( int i = 0; i < AAsize; ++i )
	{
		if( fabs( AA[i][i] ) < ALMOST_ZERO )				//TODO may be there should be fabs( ) < DELTA?
		{
			for( int row = i + 1; row < AAsize; ++row )
			{
				if( fabs( AA[row][i] ) > ALMOST_ZERO )			//TODO may be there should be fabs( ) > DELTA?
				{
					for( int col = 0; col < AAsize; ++col )
					{
						storage[col] = AA[row][col];			//TODO may be we don't need a whole vector here, just single number. (I mean storage)
						AA[row][col] = AA[i][col];
						AA[i][col] = storage[col];
					}
					tmpNum = ff[row];
					ff[row] = ff[i];
					ff[i] = tmpNum;
					break;
				}
			}
		}
	}

	for( int i = 0; i < AAsize; ++i )
	{
		if( fabs( AA[i][i] ) < ALMOST_ZERO )
		{
			cout << "alert! UFO detected!\n";
		}
	}

	//vector<vector<PL_NUM>> LL( AAsize, vector<PL_NUM>( AAsize, 0.0 ) );			//TODO here we can use less memory, UU and LL are trialgular
	//vector<vector<PL_NUM>> UU( AAsize, vector<PL_NUM>( AAsize, 0.0 ) );			//TODO initialization of arrays is slow

	//Crout's algorithm, theory is in book: Kincaid, Cheney - Numerical analysis: mathematics of scientific computing
	for( int k = 0; k < AAsize; ++k )
	{
		UU[k][k] = 1;
		for( int i = k; i < AAsize; ++i )
		{
			LL[i][k] = AA[i][k];
			for( int s = 0; s <= k - 1; ++s )
			{
				LL[i][k] -= LL[i][s] * UU[s][k];
			}
		}
		for( int j = k + 1; j < AAsize; ++j )
		{
			UU[k][j] = AA[k][j];
			for( int s = 0; s <= k - 1; ++s )
			{
				UU[k][j] -= LL[k][s] * UU[s][j];
			}
			UU[k][j] /= LL[k][k];
		}
	}

	//now we can find xx, by solving LL * zz = ff and then UU * x = zz; for theory look at the same book as for LU decomposition
	for( int i = 0; i < AAsize; ++i )
	{
		(*xx)[i] = ff[i];
		for( int  j = 0; j <= i - 1; ++j )
		{
			(*xx)[i] -= LL[i][j] * (*xx)[j];
		}
		(*xx)[i] /= LL[i][i];
	}
	for( int i = AAsize - 1; i >= 0; --i )
	{
		for( int  j = i + 1; j < AAsize; ++j )
		{
			(*xx)[i] -= UU[i][j] * (*xx)[j];
		}
	}
}


void OrthoBuilderGodunov::orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt )
{

}

void OrthoBuilderGodunov::buildSolution( vector<VarVect>* _mesh )
{
	cout << "WARNING: OrthoBuilderGodunov is void\n";
}

void OrthoBuilderGSh::orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt )		//baseV are from 0 to eq_num / 2 (including), eq_num = 10 * nx. 10 is the number of basic variables, nx is the number of lines
{
	//PL_NUM k11 = 1.414213562373095;
	//PL_NUM norm = 0;
	//vector<PL_NUM> omega2;
	//omega2.resize( eq_num * eq_num, 0.0 );

	//theory is on p.45-46 of Scott, Watts article
	if( baseV < eq_num / 2 )
	{
		for( int bvIt = 0; bvIt < baseV; ++bvIt )
		{
			for( int k = 0; k < eq_num; ++k )
			{
				solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] += (*NtoOrt)[k] * zi[n + 1][bvIt][k];
			}
			for( int k = 0; k < eq_num; ++k )
			{
				(*NtoOrt)[k] -= solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] * zi[n + 1][bvIt][k];
			}
		}
		for( int k = 0; k < eq_num; ++k )
		{
			solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] += (*NtoOrt)[k] * (*NtoOrt)[k];
		}
		solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] = sqrtl( fabs( solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] ) );

		for( int k = 0; k < eq_num; ++k )
		{
			zi[n + 1][baseV][k] = (*NtoOrt)[k] / solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV];
//			(*NtoOrt)[k] = solInfoMap[n + 1].zi[baseV][k];
		}
	}
	else
	{
		//look at the oroginal. m.b. there is a misake here
		for( int bvIt = 0; bvIt < eq_num / 2; ++bvIt )
		{
			for( int k = 0; k < eq_num; ++k )
			{
				solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] += (*NtoOrt)[k] * zi[n + 1][bvIt][k];
			}
			for( int k = 0; k < eq_num; ++k )
			{
				(*NtoOrt)[k] -= solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] * zi[n + 1][bvIt][k];
			}
		}
		for( int k = 0; k < eq_num; ++k )
		{
			//solInfoMap[n + 1].z5[k] = (*NtoOrt)[k];
			z5[n + 1][k] = (*NtoOrt)[k];
		}
	}
}


void OrthoBuilderGSh::buildSolution( vector<VarVect>* _mesh )
{
	vector<vector<PL_NUM>> M;
	vector<PL_NUM> f11;

	vector<PL_NUM> x1;
	vector<PL_NUM> res;
	vector<PL_NUM> res2;
	vector<PL_NUM> dx1;

	int msize = eq_num / 2;						//caution here!
	M.resize( msize, vector<PL_NUM>( msize, 0.0) );
	f11.resize( msize, 0.0 );
	x1.resize( msize, 0.0 );
	res.resize( msize, 0.0 );
	res2.resize( msize, 0.0 );
	dx1.resize( msize, 0.0 );
	
	//simply supported plate NO CURRENT PASSING THROUGH THE BOUNDARY

	int totLines = eq_num / EQ_NUM;
	int _a = EQ_NUM / 2;
	for( int line = 0; line < totLines; ++line )
	{
		for( int vNum = 0; vNum < eq_num / 2; ++vNum )
		{
			M[line * _a + 0][vNum] = zi[Km - 1][vNum][line * EQ_NUM + 0];		//TODO potential lags here!
			M[line * _a + 1][vNum] = zi[Km - 1][vNum][line * EQ_NUM + 1];
			M[line * _a + 2][vNum] = zi[Km - 1][vNum][line * EQ_NUM + 4];
			M[line * _a + 3][vNum] = zi[Km - 1][vNum][line * EQ_NUM + 6];
			M[line * _a + 4][vNum] = zi[Km - 1][vNum][line * EQ_NUM + 8];
		}
		/*f11[line * _a + 0] = -solInfoMap[Km - 1].z5[line * EQ_NUM + 0];
		f11[line * _a + 1] = -solInfoMap[Km - 1].z5[line * EQ_NUM + 1];
		f11[line * _a + 2] = -solInfoMap[Km - 1].z5[line * EQ_NUM + 4];
		f11[line * _a + 3] = -solInfoMap[Km - 1].z5[line * EQ_NUM + 6];
		f11[line * _a + 4] = -solInfoMap[Km - 1].z5[line * EQ_NUM + 8];*/
		f11[line * _a + 0] = -z5[Km - 1][line * EQ_NUM + 0];
		f11[line * _a + 1] = -z5[Km - 1][line * EQ_NUM + 1];
		f11[line * _a + 2] = -z5[Km - 1][line * EQ_NUM + 4];
		f11[line * _a + 3] = -z5[Km - 1][line * EQ_NUM + 6];
		f11[line * _a + 4] = -z5[Km - 1][line * EQ_NUM + 8];
	}

	LUsolve( M, f11, &x1 );

	//refinement. I do not know the theoretical source of this procedure yet. just rewrote it
	//TODO test this
	for( int i = 0; i < eq_num / 2; ++i )
	{
		res[i] = f11[i];
		for( int j = 0; j < eq_num / 2; ++j )
		{
			res[i] -= M[i][j] * x1[j];
		}
	}

	LUsolve( M, res, &dx1 );

	for( int i = 0; i < eq_num / 2; ++i )
	{
		x1[i] += dx1[i];
	}

	PL_NUM ndx = 0.0;
	PL_NUM ndx2 = 0.0;
	PL_NUM temp;

	for( int i = 0; i < eq_num / 2; ++i )
	{
		ndx += dx1[i] * dx1[i];
	}
	temp = ndx;							//FIXME may be we don't need temp

	while( fabs( ndx - ndx2 ) > 0.0 )		//FIXME may be > DELTA ??
	{
		ndx = temp;						//may be just ndx = ndx2
		for( int i = 0; i < eq_num / 2; ++i )
		{
			res2[i] = res[i];				//FIXME may be we do not need res2 here. use just res
			for( int j = 0; j < eq_num / 2; ++j )
			{
				res2[i] -= M[i][j] * dx1[j];
			}
			res[i] = res2[i];
		}
		LUsolve( M, res2, &dx1 );
		for( int i = 0; i < eq_num / 2; ++i )
		{
			x1[i] += dx1[i];
		}

		for( int i = 0; i < eq_num / 2; ++i )
		{
			res2[i] = res[i];
			for( int j = 0; j < eq_num / 2; ++j )
			{
				res2[i] -= M[i][j] * dx1[j];
			}
			res[i] = res2[i];
		}
		LUsolve( M, res2, &dx1 );
		for( int i = 0; i < eq_num / 2; ++i )
		{
			x1[i] += dx1[i];
		}


		ndx2 = 0.0;
		for( int i = 0; i < eq_num / 2; ++i )
		{
			ndx2 += dx1[i] * dx1[i];
		}
		temp = ndx2;
	}
	//refinement is over

	//now we determine coefficients for base solutions
	//the right-hand side:
	for( int i = 0; i < eq_num / 2; ++i )
	{
		solInfoMap[Km - 1].C[i] = x1[i];
	}
	//all the other points:
	for( int _x = Km - 2; _x >= 0; --_x )
	{
		for( int i = eq_num / 2 - 1; i >= 0; --i )
		{
			solInfoMap[_x].C[i] = solInfoMap[_x + 1].C[i] - solInfoMap[_x + 1].o[eq_num / 2 * ( eq_num / 2 + 1 ) / 2 + i];
			for( int j = eq_num / 2 - 1; j > i; --j )
			{
				solInfoMap[_x].C[i] -= solInfoMap[_x + 1].o[j * ( j + 1 ) / 2 + i] * solInfoMap[_x].C[j];
			}
			solInfoMap[_x].C[i] /= solInfoMap[_x + 1].o[i * ( i + 1 ) / 2 + i];
		}
	}

	//now using the coefficients we write down the solution
	for( int _x = 0; _x < Km; ++_x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			(*_mesh)[_x].Nk1[i] = 0.0;
			for( int vNum = 0; vNum < eq_num / 2; ++vNum )
			{
				(*_mesh)[_x].Nk1[i] += solInfoMap[_x].C[vNum] * zi[_x][vNum][i];			//FIXME lags may happen here
			}
			/*(*_mesh)[_x].Nk1[i] += solInfoMap[_x].z5[i];*/
			(*_mesh)[_x].Nk1[i] += z5[_x][i];
		}
	}

	//force the BCs to be zero at y == a/2
	//TODO why do we need this??
	for( int line = 0; line < eq_num / EQ_NUM; ++line )
	{
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 0] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 1] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 4] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 6] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 8] = 0.0;
	}
}