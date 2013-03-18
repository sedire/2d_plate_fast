#include "RungeKutta.h"

RungeKutta::RungeKutta( int _eq_num )
{
	eq_num = _eq_num;
	rgk_u = 0.3;
	rgk_v = 0.6;

//calculating Runge-Kutta method's parameters
	rgk_C2 = ( 2.0 * rgk_v - 1.0 ) / 12.0 / rgk_u / ( rgk_v - rgk_u ) / ( 1.0 - rgk_u );
	rgk_C3 = ( 1.0 - 2.0 * rgk_u ) / 12.0 / rgk_v / ( rgk_v - rgk_u ) / ( 1.0 - rgk_v );
	rgk_C4 = ( 6.0 * rgk_u * rgk_v - 4.0 * rgk_u - 4.0 * rgk_v + 3.0 ) / 12.0 / ( 1.0 - rgk_u ) / ( 1.0 - rgk_v );
	rgk_C1 = 1.0 - rgk_C2 - rgk_C3 - rgk_C4;
	rgk_d21 = rgk_u;
	rgk_d32 = 1.0 / 24.0 / rgk_C3 / rgk_u / ( 1.0 - rgk_v );
	rgk_d31 = rgk_v - rgk_d32;
	rgk_d43 = ( 1.0 - 2.0 * rgk_u ) / 12.0 / rgk_C4 / rgk_v / ( rgk_v - rgk_u );
	rgk_d42 = -( rgk_v * ( 4.0 * rgk_v - 5.0 ) - rgk_u + 2.0 ) / 24.0 / rgk_C4 / rgk_u / ( rgk_v - rgk_u ) / ( 1.0 - rgk_v );
	rgk_d41 = 1.0 - rgk_d42 - rgk_d43;
}

void RungeKutta::calc( /*const vector<vector<PL_NUM>>& A*/PL_NUM A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES], PL_NUM *f, PL_NUM dx, int thrNum, int hom, vector<PL_NUM>* x )
{
	//if( A.size() != eq_num || 
	//	A[0].size() != eq_num ||
	//	f.size() != eq_num ||
	//	x == 0 ||
	//	x->size() != eq_num ) {
	//	cout << "Error in Runge-Kutta calc: bad input\n";
	//	return;
	//}

	int sizeOfF = EQ_NUM * NUMBER_OF_LINES;

	for( int i = 0; i < sizeOfF; ++i ) {
		f1[thrNum][i] = 0.0;
		f2[thrNum][i] = 0.0;
		f3[thrNum][i] = 0.0;
		f4[thrNum][i] = 0.0;
	}

	for( int i = 0; i < sizeOfF; ++i ) {				//f1 = dx * Fhi( x )
		for( int j = 0; j < eq_num; ++j ) {
			f1[thrNum][i] += dx * A[i][j] * (*x)[j];
		}
	}
	for( int i = 0; i < sizeOfF; ++i ) {				//f2 = dx * Fhi( x + d21 * f1 )
		for( int j = 0; j < eq_num; ++j ) {
			f2[thrNum][i] += dx * A[i][j] * ( (*x)[j] + rgk_d21 * f1[thrNum][j] );
		}
	}
	for( int i = 0; i < sizeOfF; ++i ) {				//f3 = dx * Fhi( x + d31 * f1 + d32 * f2 )
		for( int j = 0; j < eq_num; ++j ) {
			f3[thrNum][i] += dx * A[i][j] * ( (*x)[j] + rgk_d31 * f1[thrNum][j] + rgk_d32 * f2[thrNum][j] );
		}
	}
	for( int i = 0; i < sizeOfF; ++i ) {				//f2 = dx * Fhi( x + d41 * f1 + d42 * f2 + d43 * f3 )
		for( int j = 0; j < eq_num; ++j ) {
			f4[thrNum][i] += dx * A[i][j] * ( (*x)[j] + rgk_d41 * f1[thrNum][j] + rgk_d42 * f2[thrNum][j] + rgk_d43 * f3[thrNum][j] );
		}
	}
	if( hom != 0 ) {
		for( int i = 0; i < sizeOfF; ++i ) {
			f1[thrNum][i] += dx * f[i];
			f2[thrNum][i] += dx * f[i];
			f3[thrNum][i] += dx * f[i];
			f4[thrNum][i] += dx * f[i];
		}
	}

	for( int i = 0; i < sizeOfF; ++i ) {
		(*x)[i] = (*x)[i] + rgk_C1 * f1[thrNum][i] + rgk_C2 * f2[thrNum][i] + rgk_C3 * f3[thrNum][i] + rgk_C4 * f4[thrNum][i];
	}
}