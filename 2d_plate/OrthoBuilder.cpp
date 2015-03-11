#include "OrthoBuilder.h"

SolInfo::SolInfo() {}
SolInfo::~SolInfo() {}

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

OrthoBuilder::~OrthoBuilder() {}
OrthoBuilder::OrthoBuilder( const int _varNum, const int _Km ) :
	varNum( _varNum ),
	Km( _Km ),
	zi( NODES_ON_Y, vector<vector<PL_NUM> >( EQ_NUM * NUMBER_OF_LINES / 2, vector<PL_NUM>( EQ_NUM * NUMBER_OF_LINES, 0.0 ) ) ),
	z5( NODES_ON_Y, vector<PL_NUM>( EQ_NUM * NUMBER_OF_LINES, 0.0 ) )
{

}

OrthoBuilderGSh::OrthoBuilderGSh( int _varNum, int _Km ) :
	OrthoBuilder( _varNum, _Km )
{
	LL.resize( varNum / 2, vector<PL_NUM>( varNum / 2, 0.0 ) );
	UU.resize( varNum / 2, vector<PL_NUM>( varNum / 2, 0.0 ) );
}

OrthoBuilderGodunov::OrthoBuilderGodunov( int _varNum, int _Km ):
	OrthoBuilder( _varNum, _Km )
{

}

void OrthoBuilder::flushO( int x )
{
	solInfoMap[x].flushO();
}

void OrthoBuilder::setParams()
{
	try 
	{
		solInfoMap.resize( Km );
		for( int i = 0; i < solInfoMap.size(); ++i )
		{
			solInfoMap[i].setup( varNum );
		}
	}
	catch( bad_alloc &ba )
	{
		cout << ba.what() << endl;
		std::abort();
	}

	orthoDone.resize( Km - 1, 0 );

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
	cout << "Warning! I am void!\n";
}

void OrthoBuilder::setOrthoDoneInfo( int y )
{
	orthoDone[y] = true;
}

void OrthoBuilder::resetOrthoDoneInfo()
{
	for( int y = 0; y < orthoDone.size(); ++y )
	{
		orthoDone[y] = false;
	}
}

void OrthoBuilder::LUsolve( vector<vector<PL_NUM>>& AA, vector<PL_NUM>& ff, vector<PL_NUM>* xx )
{
	if( AA.size() < 1 )
	{
		cout << "ERROR in LUsolve: AA.size() is less than 1\n";
		return;
	}
	if( ff.size() < 1 )
	{
		cout << "ERROR in LUsolve: ff.size() is less than 1\n";
		return;
	}
	if( xx == 0 )
	{
		cout << "ERROR in LUsolve: xx is null\n";
		return;
	}
	size_t AAsize = AA.size();

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
			cout << "WARNING in LUsolve: AA has less than ALMOST_ZERO on diagonal and may be singular\n";
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


void OrthoBuilderGodunov::orthonorm( int baseV, int n, PL_NUM* NtoOrt )
{
	cout << "Warning: I am void\n";
}

void OrthoBuilderGodunov::buildSolution( vector<VarVect>* _mesh )
{
	cout << "WARNING: OrthoBuilderGodunov is void\n";
}

//void OrthoBuilderGSh::calcScalarProdsPar( int baseV, int n, vector<PL_NUM>* NtoOrt )
//{
//	int begIt = omp_get_thread_num() * ( varNum / NUM_OF_THREADS + 1 );
//	int endIt = begIt + varNum / NUM_OF_THREADS + 1;
//	if( omp_get_thread_num() == NUM_OF_THREADS - 1 )
//	{
//		endIt = varNum;
//	}
//
//	omegaPar[omp_get_thread_num()] = (*NtoOrt)[begIt] * solInfoMap[n + 1].zi[baseV][begIt];
//	for( int k = begIt + 1; k < endIt; ++k )
//	{
//		omegaPar[omp_get_thread_num()] += (*NtoOrt)[k] * solInfoMap[n + 1].zi[baseV][k];
//	}
//}

//void OrthoBuilderGSh::calcScalarProdsPar2( int baseV, int n, vector<PL_NUM>* NtoOrt )
//{
//	int begIt = omp_get_thread_num() * ( varNum / NUM_OF_THREADS + 1 );
//	int endIt = begIt + varNum / NUM_OF_THREADS + 1;
//	if( omp_get_thread_num() == NUM_OF_THREADS - 1 )
//	{
//		endIt = varNum;
//	}
//
//	omegaPar[omp_get_thread_num()] = solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] * solInfoMap[n + 1].zi[baseV][begIt];
//	for( int k = begIt + 1; k < endIt; ++k )
//	{
//		omegaPar[omp_get_thread_num()] = solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] * solInfoMap[n + 1].zi[baseV][begIt];
//	}
//}

inline void OrthoBuilderGSh::setNextSolVects( int n, const PL_NUM decompVect[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] )
{
	for( int vNum = 0; vNum < EQ_NUM * NUMBER_OF_LINES / 2; ++vNum )
	{
		for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
		{
			zi[n + 1][vNum][j] = decompVect[vNum][j];
		}
	}
	for( int j = 0; j < EQ_NUM * NUMBER_OF_LINES; ++j )
	{
		z5[n + 1][j] = decompVect[EQ_NUM * NUMBER_OF_LINES / 2][j];
	}
}

inline PL_NUM OrthoBuilder::getInfNorm( PL_NUM* vect, int vectSize )
{
	if( vect != 0 )
	{
		PL_NUM ret = fabs( vect[0] );
		for( int i = 1; i < vectSize; ++i )
		{
			if( fabs( vect[i] ) > ret )
			{
				ret = fabs( vect[i] );
			}
		}
		return ret;
	}
	else
	{
		return -1;
	}
}

inline int OrthoBuilderGSh::checkOrtho( int n, 
									PL_NUM vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES], 
									PL_NUM vectSetOrig[EQ_NUM * NUMBER_OF_LINES / 2 + 1][EQ_NUM * NUMBER_OF_LINES] )
{
	int ret = 0;
	PL_NUM eps = ORTHONORM_CHECK_EPS;

	//check for the main basis vectors
	for( int vNum = 1; vNum < EQ_NUM * NUMBER_OF_LINES / 2; ++vNum )
	{
		if( getInfNorm( vectSetOrtho[vNum], EQ_NUM * NUMBER_OF_LINES ) * solInfoMap[n + 1].o[vNum * ( vNum + 1 ) / 2 + vNum] <
			eps * getInfNorm( vectSetOrig[vNum], EQ_NUM * NUMBER_OF_LINES ) )
		{
			ret = 1;
			break;
		}
	}

	if( ret != 1 )	//check for the vector of particular solution
	{
		if( getInfNorm( vectSetOrtho[EQ_NUM * NUMBER_OF_LINES / 2], EQ_NUM * NUMBER_OF_LINES ) <
			eps * getInfNorm( vectSetOrig[EQ_NUM * NUMBER_OF_LINES / 2], EQ_NUM * NUMBER_OF_LINES ) )
		{
			ret = 1;
		}
	}

	return ret;
}

void OrthoBuilderGSh::orthonorm( int baseV, int n, PL_NUM* NtoOrt )		//baseV are from 0 to varNum / 2 (including), varNum = 10 * nx. 10 is the number of basic variables, nx is the number of lines
{
	PL_NUM k11 = 1.414213562373095;
	PL_NUM norm( 0.0 );

	//theory is on p.45-46 of Scott, Watts article
	if( baseV < varNum / 2 )
	{
		norm = 0.0;
		for( int i = 0; i < varNum; ++i )
		{
			norm += NtoOrt[i] * NtoOrt[i];
			omega2[i] = 0.0;
		}
		norm = sqrtf( norm );

		for( int bvIt = 0; bvIt < baseV; ++bvIt )
		{
			for( int k = 0; k < varNum; ++k )
			{
				//solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] += solInfoMap[n + 1].zi[baseV][k] * solInfoMap[n + 1].zi[bvIt][k];			//problems here
				solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] += NtoOrt[k] * zi[n + 1][bvIt][k];			//problems here
			}
			for( int k = 0; k < varNum; ++k )
			{
				NtoOrt[k] -= solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] * zi[n + 1][bvIt][k];
			}
		}

		for( int k = 0; k < varNum; ++k )
		{
			solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] += NtoOrt[k] * NtoOrt[k];			//problems here
		}
		solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] = sqrtl( fabs( solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] ) );

		if( norm / solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] <= k11 )
		{
			for( int k = 0; k < varNum; ++k )
			{
				zi[n + 1][baseV][k] = NtoOrt[k] / solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV];
			}
		}
		else
		{
			//cout << " === no ortho!!\n";
			for( int bvIt = 0; bvIt < baseV; ++bvIt )
			{
				for( int k = 0; k < varNum; ++k )
				{
					omega2[bvIt] += NtoOrt[k] * zi[n + 1][bvIt][k];
				}
				for( int k = 0; k < varNum; ++k )
				{
					NtoOrt[k] -= omega2[bvIt] * zi[n + 1][bvIt][k];
				}
			}
			for( int k = 0; k < varNum; ++k )
			{
				omega2[baseV] += NtoOrt[k] * NtoOrt[k];
			}
			omega2[baseV] = sqrtl( fabs( omega2[baseV] ) );
			for( int k = 0; k < varNum; ++k )
			{
				zi[n + 1][baseV][k] = NtoOrt[k] / omega2[baseV];
				//	(*NtoOrt)[k] = solInfoMap[n + 1].zi[baseV][k];
			}
			for( int bvIt = 0; bvIt < baseV; ++bvIt )
			{
				solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] += omega2[bvIt];
			}
			solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + baseV] = omega2[baseV];
		}
		//caution: previously there was a k11-procedure possibly from the modified gram-schmidt
	}
	else
	{
		//look at the original. m.b. there is a misake here
		for( int bvIt = 0; bvIt < varNum / 2; ++bvIt )
		{
			for( int k = 0; k < varNum; ++k )
			{
				solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] += NtoOrt[k] * zi[n + 1][bvIt][k];
			}
			for( int k = 0; k < varNum; ++k )
			{
				NtoOrt[k] -= solInfoMap[n + 1].o[baseV * ( baseV + 1 ) / 2 + bvIt] * zi[n + 1][bvIt][k];
			}
		}
		for( int k = 0; k < varNum; ++k )
		{
			z5[n + 1][k] = NtoOrt[k];
		}
	}
}

void OrthoBuilderGSh::buildSolution( vector<VarVect>* _mesh )
{
	static const int msize = NUMBER_OF_LINES * EQ_NUM / 2;
	Matrix<PL_NUM, msize, msize, RowMajor> M;
	Matrix<PL_NUM, msize, 1> f11;
	Matrix<PL_NUM, msize, 1> x1;
	Matrix<PL_NUM, msize, 1> res;
	Matrix<PL_NUM, msize, 1> res2;
	Matrix<PL_NUM, msize, 1> dx1;

	for( int i = 0; i < msize; ++i )
	{
		for( int j = 0; j < msize; ++j )
		{
			M( j, i ) = 0.0;
		}
		f11( i ) = 0.0;
		x1( i ) = 0.0;
		res( i ) = 0.0;
		res2( i ) = 0.0;
		dx1( i ) = 0.0;
	}
	
	//simply supported plate NO CURRENT PASSING THROUGH THE BOUNDARY

	int totLines = varNum / EQ_NUM;
	int _a = EQ_NUM / 2;
	for( int line = 0; line < totLines; ++line )
	{
		for( int vNum = 0; vNum < varNum / 2; ++vNum )
		{
			M( line * _a + 0, vNum ) = zi[Km - 1][vNum][line * EQ_NUM + 0];		//TODO potential lags here!
			M( line * _a + 1, vNum ) = zi[Km - 1][vNum][line * EQ_NUM + 1];
			M( line * _a + 2, vNum ) = zi[Km - 1][vNum][line * EQ_NUM + 4];
			M( line * _a + 3, vNum ) = zi[Km - 1][vNum][line * EQ_NUM + 6];
			M( line * _a + 4, vNum ) = zi[Km - 1][vNum][line * EQ_NUM + 8];
		}
		f11( line * _a + 0 ) = -z5[Km - 1][line * EQ_NUM + 0];
		f11( line * _a + 1 ) = -z5[Km - 1][line * EQ_NUM + 1];
		f11( line * _a + 2 ) = -z5[Km - 1][line * EQ_NUM + 4];
		f11( line * _a + 3 ) = -z5[Km - 1][line * EQ_NUM + 6];
		f11( line * _a + 4 ) = -z5[Km - 1][line * EQ_NUM + 8];
	}

	//EigenSolver<Matrix<PL_NUM, msize, msize, RowMajor>> es( M );
	//if( es.info() == Success )
	//{
	//	Matrix< complex<PL_NUM>, msize, 1> eigv = es.eigenvalues();

	//	PL_NUM minL = sqrt( eigv( 0 ).imag() * eigv( 0 ).imag() + eigv( 0 ).real() * eigv( 0 ).real() );
	//	PL_NUM maxL = sqrt( eigv( 0 ).imag() * eigv( 0 ).imag() + eigv( 0 ).real() * eigv( 0 ).real() );
	//	for( int ii = 1; ii < eigv.size(); ++ii )
	//	{
	//		if( sqrt( eigv( ii ).imag() * eigv( ii ).imag() + eigv( ii ).real() * eigv( ii ).real() ) > maxL )
	//		{
	//			maxL = sqrt( eigv( ii ).imag() * eigv( ii ).imag() + eigv( ii ).real() * eigv( ii ).real() );
	//		}
	//		if( sqrt( eigv( ii ).imag() * eigv( ii ).imag() + eigv( ii ).real() * eigv( ii ).real() ) < minL )
	//		{
	//			minL = sqrt( eigv( ii ).imag() * eigv( ii ).imag() + eigv( ii ).real() * eigv( ii ).real() );
	//		}
	//	}
	//	cout << "cond number is " << maxL / minL << endl;
	//}
	x1 = M.fullPivLu().solve( f11 );

	//refinement. I do not know the theoretical source of this procedure yet. just rewrote it
	//TODO test this
	res = f11 - M * x1;
	dx1 = M.fullPivLu().solve( res );
	x1 = x1 + dx1;

	PL_NUM nx = 0.0;
	PL_NUM ndx = 0.0;
	PL_NUM ndx2 = 0.0;
	PL_NUM temp;

	ndx = dx1.lpNorm<Infinity>();
	temp = ndx;							//FIXME may be we don't need temp
	int iterCount = 0;

	do
	{
		ndx = temp;
		++iterCount;

		res2 = res - M * dx1;
		res = res2;
		dx1 = M.fullPivLu().solve( res2 );
		x1 = x1 + dx1;


		res2 = res - M * dx1;
		res = res2;
		dx1 = M.fullPivLu().solve( res2 );
		x1 = x1 + dx1;

		nx = x1.lpNorm<Infinity>();
		ndx2 = dx1.lpNorm<Infinity>();
		temp = ndx2;
	} while( ndx2 < 0.9 * ndx && ndx2 / nx >= 2 * EPS_W );
	cout << " " << iterCount << " refnmt iterations\n";
	cout << " refnmt relative error is " << (PL_NUM)( ( M * x1 - f11 ).norm() / f11.norm() ) << endl;
	//refinement is over

	//now we determine coefficients for base solutions
	//the right-hand side:
	for( int i = 0; i < varNum / 2; ++i )
	{
		solInfoMap[Km - 1].C[i] = x1( i );
	}
	//all the other points:
	for( int _x = Km - 2; _x >= 0; --_x )
	{
		for( int i = varNum / 2 - 1; i >= 0; --i )
		{
			solInfoMap[_x].C[i] = solInfoMap[_x + 1].C[i] - solInfoMap[_x + 1].o[varNum / 2 * ( varNum / 2 + 1 ) / 2 + i];
			for( int j = varNum / 2 - 1; j > i; --j )
			{
				solInfoMap[_x].C[i] -= solInfoMap[_x + 1].o[j * ( j + 1 ) / 2 + i] * solInfoMap[_x].C[j];
			}
			solInfoMap[_x].C[i] /= solInfoMap[_x + 1].o[i * ( i + 1 ) / 2 + i];
		}
	}

	//now using the coefficients we write down the solution
	for( int _x = 0; _x < Km; ++_x )
	{
		for( int i = 0; i < varNum; ++i )
		{
			(*_mesh)[_x].Nk1[i] = 0.0;
			for( int vNum = 0; vNum < varNum / 2; ++vNum )
			{
				(*_mesh)[_x].Nk1[i] += solInfoMap[_x].C[vNum] * zi[_x][vNum][i];			//FIXME lags may happen here
			}
			/*(*_mesh)[_x].Nk1[i] += solInfoMap[_x].z5[i];*/
			(*_mesh)[_x].Nk1[i] += z5[_x][i];
		}
	}

	//force the BCs to be zero at y == a/2
	//TODO why do we need this??
	for( int line = 0; line < varNum / EQ_NUM; ++line )
	{
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 0] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 1] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 4] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 6] = 0.0;
		(*_mesh)[Km - 1].Nk1[line * EQ_NUM + 8] = 0.0;
	}
}