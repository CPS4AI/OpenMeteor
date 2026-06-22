
#ifndef PRECOMPUTE_H
#define PRECOMPUTE_H

#pragma once
#include "globals.h"


class Precompute
{
private:
	void initialize();

public:
	Precompute();
	~Precompute();

	void getRandomShares(RSSVectorMyType &r, size_t size);
	void getRandomShares(RSSVectorSmallType &r, size_t size);
	void getDividedShares(RSSVectorMyType &r, RSSVectorMyType &rPrime, int d, size_t size);
	void getZeroShares(RSSVectorMyType &z, size_t size);
	void getRandomBitShares(RSSVectorSmallType &a, size_t size);
	void getProductShares(const RSSVectorMyType &a, const RSSVectorMyType &b, RSSVectorMyType &c, size_t size);
	void getMatrixProductShares(const RSSVectorMyType &a, const RSSVectorMyType &b, RSSVectorMyType &c,
								size_t rows, size_t common_dim, size_t columns,
								size_t transpose_a, size_t transpose_b);
	void getProductShares(const RSSVectorSmallType &a, const RSSVectorSmallType &b, RSSVectorSmallType &c, size_t size);
	void getBitProductShares(const RSSVectorSmallType &a, const RSSVectorSmallType &b, RSSVectorSmallType &c, size_t size);
	void getBitToArithmeticShares(const RSSVectorSmallType &bits, RSSVectorMyType &arith, size_t size);
	void getBitToPrimeShares(const RSSVectorSmallType &bits, RSSVectorSmallType &arith, size_t size);
	void getTruncationMasks(RSSVectorMyType &r, RSSVectorMyType &rTrunc, size_t power, size_t size);
	void getNonZeroPrimeJShares(MEVectorSmallType &gamma, size_t size);
	void getPreMSBObjects(RSSVectorMyType &lambda, MEVectorSmallType &lambdaBits, MEVectorSmallType &lambdaMSB, MEVectorSmallType &sBits,
						  MEVectorSmallType &betaBit, MEVectorSmallType &betaPrime,
						  MEVectorSmallType &gamma, size_t size);
	void getSelectorBitShares(RSSVectorSmallType &c, RSSVectorMyType &m_c, size_t size);
	void getShareConvertObjects(RSSVectorMyType &r, RSSVectorSmallType &shares_r, RSSVectorSmallType &alpha, size_t size);
	//void getTriplets(RSSVectorMyType &a, RSSVectorMyType &b, RSSVectorMyType &c, 
	//				size_t rows, size_t common_dim, size_t columns);
	//void getTriplets(RSSVectorMyType &a, RSSVectorMyType &b, RSSVectorMyType &c, size_t size);
	//void getTriplets(RSSVectorSmallType &a, RSSVectorSmallType &b, RSSVectorSmallType &c, size_t size);

};


#endif
