#include "RungeKutta.h"

RungeKutta::RungeKutta( int _varNum )
{
	varNum = _varNum;
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

	buildBorders();
}

RungeKutta::~RungeKutta() {}

void RungeKutta::buildBorders()
{
	for( int i = 0; i < EQ_NUM * NUMBER_OF_LINES; ++i )
	{
		if( i < 2 * EQ_NUM )
		{
			lb[i] = 0;
			rb[i] = EQ_NUM * 4;
		}
		else if( i >= 2 * EQ_NUM && i < EQ_NUM * ( NUMBER_OF_LINES - 2 ) )
		{
			lb[i] = ( i / EQ_NUM - 2 ) * EQ_NUM;
			rb[i] = lb[i] + EQ_NUM * 5;
		}
		else if( i >= EQ_NUM * ( NUMBER_OF_LINES - 2 ) )
		{
			lb[i] = EQ_NUM * NUMBER_OF_LINES - 4 * EQ_NUM;
			rb[i] = EQ_NUM * NUMBER_OF_LINES;
		}
	}
}

//void RungeKutta::calc( PL_NUM A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES], PL_NUM *f, PL_NUM dx, int thrNum, int hom, const vector<PL_NUM>& x, PL_NUM *x1 )	//this one is optimized to take into accout matrix block structure
//{
//	//if( A.size() != varNum || 
//	//	A[0].size() != varNum ||
//	//	f.size() != varNum ||
//	//	x == 0 ||
//	//	x->size() != varNum ) {
//	//	cout << "Error in Runge-Kutta calc: bad input\n";
//	//	return;
//	//}
//
//	int sizeOfF = EQ_NUM * NUMBER_OF_LINES;
//
//	for( int i = 0; i < sizeOfF; ++i ) {
//		f1[thrNum][i] = 0.0;
//	}
//
//	//this is just rolled out RgK
//	PL_NUM tmpProd = 0.0;
//	for( int i = 0; i < sizeOfF; ++i ) {				//f1 = dx * Fhi( x )
//		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
//			f1[thrNum][i] += A[i][j] * (x)[j];
//			f1[thrNum][i] += A[i][j + 1] * (x)[j + 1];
//			f1[thrNum][i] += A[i][j + 2] * (x)[j + 2];
//			f1[thrNum][i] += A[i][j + 3] * (x)[j + 3];
//			f1[thrNum][i] += A[i][j + 4] * (x)[j + 4];
//		}
//		f1[thrNum][i] *= dx;
//		f2[thrNum][i] = f1[thrNum][i];
//		f3[thrNum][i] = f1[thrNum][i];
//		f4[thrNum][i] = f1[thrNum][i];
//	}
//	for( int i = 0; i < sizeOfF; ++i ) {				//f2 = dx * Fhi( x + d21 * f1 )
//		tmpProd = 0.0;
//		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
//			tmpProd += A[i][j] * f1[thrNum][j];
//			tmpProd += A[i][j + 1] * f1[thrNum][j + 1];
//			tmpProd += A[i][j + 2] * f1[thrNum][j + 2];
//			tmpProd += A[i][j + 3] * f1[thrNum][j + 3];
//			tmpProd += A[i][j + 4] * f1[thrNum][j + 4];
//		}
//		f2[thrNum][i] += dx * rgk_d21 * tmpProd;
//		f3[thrNum][i] += dx * rgk_d31 * tmpProd;
//		f4[thrNum][i] += dx * rgk_d41 * tmpProd;
//	}
//	for( int i = 0; i < sizeOfF; ++i ) {				//f3 = dx * Fhi( x + d31 * f1 + d32 * f2 )
//		tmpProd = 0.0;
//		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
//			tmpProd += A[i][j] * f2[thrNum][j];
//			tmpProd += A[i][j + 1] * f2[thrNum][j + 1];
//			tmpProd += A[i][j + 2] * f2[thrNum][j + 2];
//			tmpProd += A[i][j + 3] * f2[thrNum][j + 3];
//			tmpProd += A[i][j + 4] * f2[thrNum][j + 4];
//		}
//		f3[thrNum][i] += dx * rgk_d32 * tmpProd;
//		f4[thrNum][i] += dx * rgk_d42 * tmpProd;
//	}
//	for( int i = 0; i < sizeOfF; ++i ) {				//f2 = dx * Fhi( x + d41 * f1 + d42 * f2 + d43 * f3 )
//		tmpProd = 0.0;
//		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
//			tmpProd += A[i][j] * f3[thrNum][j];
//			tmpProd += A[i][j + 1] * f3[thrNum][j + 1];
//			tmpProd += A[i][j + 2] * f3[thrNum][j + 2];
//			tmpProd += A[i][j + 3] * f3[thrNum][j + 3];
//			tmpProd += A[i][j + 4] * f3[thrNum][j + 4];
//		}
//		f4[thrNum][i] += dx * rgk_d43 * tmpProd;
//	}
//	if( hom != 0 ) {
//		for( int i = 0; i < sizeOfF; ++i ) {
//			f1[thrNum][i] += dx * f[i];
//			f2[thrNum][i] += dx * f[i];
//			f3[thrNum][i] += dx * f[i];
//			f4[thrNum][i] += dx * f[i];
//		}
//	}
//
//	for( int i = 0; i < sizeOfF; ++i ) {
//		(x1)[i] = (x)[i] + rgk_C1 * f1[thrNum][i] + rgk_C2 * f2[thrNum][i] + rgk_C3 * f3[thrNum][i] + rgk_C4 * f4[thrNum][i];
//	}
//}

void RungeKutta::calc( PL_NUM A[EQ_NUM * NUMBER_OF_LINES][EQ_NUM * NUMBER_OF_LINES], PL_NUM *f, PL_NUM dx, int thrNum, int hom, const vector<PL_NUM>& x, PL_NUM *x1 )	//this one is optimized to take into accout matrix block structure
{
	//if( A.size() != varNum || 
	//	A[0].size() != varNum ||
	//	f.size() != varNum ||
	//	x == 0 ||
	//	x->size() != varNum ) {
	//	cout << "Error in Runge-Kutta calc: bad input\n";
	//	return;
	//}

	int sizeOfF = EQ_NUM * NUMBER_OF_LINES;

	for( int i = 0; i < sizeOfF; ++i ) {
		f1[thrNum][i] = 0.0;
	}

	//this is just rolled out RgK
	PL_NUM tmpProd = 0.0;
	for( int i = 0; i < sizeOfF; ++i ) {				//f1 = dx * Fhi( x )
		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
			f1[thrNum][i] += A[i][j] * (x)[j];
			f1[thrNum][i] += A[i][j + 1] * (x)[j + 1];
			f1[thrNum][i] += A[i][j + 2] * (x)[j + 2];
			f1[thrNum][i] += A[i][j + 3] * (x)[j + 3];
			f1[thrNum][i] += A[i][j + 4] * (x)[j + 4];
		}
		f1[thrNum][i] *= dx;
		f2[thrNum][i] = f1[thrNum][i];
		f3[thrNum][i] = f1[thrNum][i];
		f4[thrNum][i] = f1[thrNum][i];
	}
	for( int i = 0; i < sizeOfF; ++i ) {				//f2 = dx * Fhi( x + d21 * f1 )
		tmpProd = 0.0;
		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
			tmpProd += A[i][j] * f1[thrNum][j];
			tmpProd += A[i][j + 1] * f1[thrNum][j + 1];
			tmpProd += A[i][j + 2] * f1[thrNum][j + 2];
			tmpProd += A[i][j + 3] * f1[thrNum][j + 3];
			tmpProd += A[i][j + 4] * f1[thrNum][j + 4];
		}
		f2[thrNum][i] += dx * rgk_d21 * tmpProd;
		f3[thrNum][i] += dx * rgk_d31 * tmpProd;
		f4[thrNum][i] += dx * rgk_d41 * tmpProd;
	}
	for( int i = 0; i < sizeOfF; ++i ) {				//f3 = dx * Fhi( x + d31 * f1 + d32 * f2 )
		tmpProd = 0.0;
		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
			tmpProd += A[i][j] * f2[thrNum][j];
			tmpProd += A[i][j + 1] * f2[thrNum][j + 1];
			tmpProd += A[i][j + 2] * f2[thrNum][j + 2];
			tmpProd += A[i][j + 3] * f2[thrNum][j + 3];
			tmpProd += A[i][j + 4] * f2[thrNum][j + 4];
		}
		f3[thrNum][i] += dx * rgk_d32 * tmpProd;
		f4[thrNum][i] += dx * rgk_d42 * tmpProd;
	}
	for( int i = 0; i < sizeOfF; ++i ) {				//f2 = dx * Fhi( x + d41 * f1 + d42 * f2 + d43 * f3 )
		tmpProd = 0.0;
		for( int j = lb[i]; j < rb[i]; j = j + 5 ) {
			tmpProd += A[i][j] * f3[thrNum][j];
			tmpProd += A[i][j + 1] * f3[thrNum][j + 1];
			tmpProd += A[i][j + 2] * f3[thrNum][j + 2];
			tmpProd += A[i][j + 3] * f3[thrNum][j + 3];
			tmpProd += A[i][j + 4] * f3[thrNum][j + 4];
		}
		f4[thrNum][i] += dx * rgk_d43 * tmpProd;
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
		(x1)[i] = (x)[i] + rgk_C1 * f1[thrNum][i] + rgk_C2 * f2[thrNum][i] + rgk_C3 * f3[thrNum][i] + rgk_C4 * f4[thrNum][i];
	}
}