
#pragma once
#include "Precompute.h"
#include "AESObject.h"
#include "connect.h"
#include "tools.h"
#include <thread>

extern AESObject* aes_next;
extern AESObject* aes_prev;
extern AESObject* aes_indep;
extern void push_preprocessing_cost();
extern void pop_preprocessing_cost();

namespace
{
	struct PreprocessingCostScope
	{
		PreprocessingCostScope() { push_preprocessing_cost(); }
		~PreprocessingCostScope() { pop_preprocessing_cost(); }
	};

	void reconstructBitSharesLocal(const RSSVectorSmallType &a, vector<smallType> &b, size_t size)
	{
		assert(a.size() == size && "size mismatch for reconstructBitSharesLocal");
		assert(b.size() == size && "size mismatch for reconstructBitSharesLocal");

		vector<smallType> a_next(size), a_prev(size);
		for (size_t i = 0; i < size; ++i)
		{
			a_next[i] = a[i].first;
			b[i] = a[i].first ^ a[i].second;
		}

		thread *threads = new thread[2];
		threads[0] = thread(sendVector<smallType>, ref(a_next), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<smallType>, ref(a_prev), prevParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (size_t i = 0; i < size; ++i)
			b[i] = (b[i] ^ a_prev[i]) & 1;
	}

	void reconstructPrimeSharesLocal(const RSSVectorSmallType &a, vector<smallType> &b, size_t size)
	{
		assert(a.size() == size && "size mismatch for reconstructPrimeSharesLocal");
		assert(b.size() == size && "size mismatch for reconstructPrimeSharesLocal");

		vector<smallType> a_next(size), a_prev(size);
		for (size_t i = 0; i < size; ++i)
		{
			a_next[i] = a[i].first;
			b[i] = additionModPrime[a[i].first][a[i].second];
		}

		thread *threads = new thread[2];
		threads[0] = thread(sendVector<smallType>, ref(a_next), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<smallType>, ref(a_prev), prevParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (size_t i = 0; i < size; ++i)
			b[i] = additionModPrime[b[i]][a_prev[i]];
	}

	void reshareAdditiveRing(const vector<myType> &additive, RSSVectorMyType &out, size_t size)
	{
		assert(additive.size() == size && "size mismatch for reshareAdditiveRing");
		assert(out.size() == size && "size mismatch for reshareAdditiveRing");

		vector<myType> localRandom(size), prevRandom(size), localShare(size), nextShare(size);
		for (size_t i = 0; i < size; ++i)
			localRandom[i] = aes_indep->get64Bits();

		thread *threads = new thread[2];
		threads[0] = thread(sendVector<myType>, ref(localRandom), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<myType>, ref(prevRandom), prevParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (size_t i = 0; i < size; ++i)
			localShare[i] = additive[i] + localRandom[i] - prevRandom[i];

		threads = new thread[2];
		threads[0] = thread(sendVector<myType>, ref(localShare), prevParty(partyNum), size);
		threads[1] = thread(receiveVector<myType>, ref(nextShare), nextParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (size_t i = 0; i < size; ++i)
			out[i] = make_pair(localShare[i], nextShare[i]);
	}

	void reshareAdditivePrime(const vector<smallType> &additive, RSSVectorSmallType &out, size_t size)
	{
		assert(additive.size() == size && "size mismatch for reshareAdditivePrime");
		assert(out.size() == size && "size mismatch for reshareAdditivePrime");

		vector<smallType> localRandom(size), prevRandom(size), localShare(size), nextShare(size);
		for (size_t i = 0; i < size; ++i)
			localRandom[i] = aes_indep->randModPrime();

		thread *threads = new thread[2];
		threads[0] = thread(sendVector<smallType>, ref(localRandom), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<smallType>, ref(prevRandom), prevParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (size_t i = 0; i < size; ++i)
			localShare[i] = subtractModPrime[additionModPrime[additive[i]][localRandom[i]]][prevRandom[i]];

		threads = new thread[2];
		threads[0] = thread(sendVector<smallType>, ref(localShare), prevParty(partyNum), size);
		threads[1] = thread(receiveVector<smallType>, ref(nextShare), nextParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (size_t i = 0; i < size; ++i)
			out[i] = make_pair(localShare[i], nextShare[i]);
	}

	void maskBooleanJFromRSS(Precompute &precompute, const RSSVectorSmallType &secret,
							 MEVectorSmallType &out, size_t size)
	{
		assert(secret.size() == size && "size mismatch for maskBooleanJFromRSS");
		assert(out.size() == size && "size mismatch for maskBooleanJFromRSS");

		RSSVectorSmallType outputMask(size), deltaShares(size);
		vector<smallType> delta(size);
		precompute.getRandomBitShares(outputMask, size);
		for (size_t i = 0; i < size; ++i)
		{
			deltaShares[i].first = (secret[i].first ^ outputMask[i].first) & 1;
			deltaShares[i].second = (secret[i].second ^ outputMask[i].second) & 1;
		}

		reconstructBitSharesLocal(deltaShares, delta, size);
		for (size_t i = 0; i < size; ++i)
			out[i] = make_pair(delta[i], outputMask[i]);
	}

	void maskPrimeJFromRSS(Precompute &precompute, const RSSVectorSmallType &secret,
						   MEVectorSmallType &out, size_t size)
	{
		assert(secret.size() == size && "size mismatch for maskPrimeJFromRSS");
		assert(out.size() == size && "size mismatch for maskPrimeJFromRSS");

		RSSVectorSmallType outputMask(size), deltaShares(size);
		vector<smallType> delta(size);
		precompute.getRandomShares(outputMask, size);
		for (size_t i = 0; i < size; ++i)
		{
			deltaShares[i].first = subtractModPrime[secret[i].first][outputMask[i].first];
			deltaShares[i].second = subtractModPrime[secret[i].second][outputMask[i].second];
		}

		reconstructPrimeSharesLocal(deltaShares, delta, size);
		for (size_t i = 0; i < size; ++i)
			out[i] = make_pair(delta[i], outputMask[i]);
	}
}

Precompute::Precompute(){initialize();}
Precompute::~Precompute(){}
void Precompute::initialize(){}

void Precompute::getRandomShares(RSSVectorMyType &r, size_t size)
{
	PreprocessingCostScope scope;
	assert(r.size() == size && "size mismatch for getRandomShares");
	for (size_t i = 0; i < size; ++i)
		r[i] = std::make_pair(aes_prev->get64Bits(), aes_next->get64Bits());
}

void Precompute::getRandomShares(RSSVectorSmallType &r, size_t size)
{
	PreprocessingCostScope scope;
	assert(r.size() == size && "size mismatch for getRandomShares");
	for (size_t i = 0; i < size; ++i)
		r[i] = std::make_pair(aes_prev->randModPrime(), aes_next->randModPrime());
}

// Currently, r = 3 and rPrime = 3 * 2^d
// TODO: ReLU produces a bug with this. Why? funcRELU does not even call getDividedShares()
void Precompute::getDividedShares(RSSVectorMyType &r, RSSVectorMyType &rPrime, int d, size_t size)
{
	PreprocessingCostScope scope;
	assert(r.size() == size && "r.size is incorrect");
	assert(rPrime.size() == size && "rPrime.size is incorrect");

	for (int i = 0; i < size; ++i)
	{
		r[i].first = 0;
		r[i].second = 0;
		rPrime[i].first = 0;
		rPrime[i].second = 0;
		// r[i].first = 1;
		// r[i].second = 1;
		// rPrime[i].first = d;
		// rPrime[i].second = d;		
	}
}

void Precompute::getZeroShares(RSSVectorMyType &z, size_t size)
{
	PreprocessingCostScope scope;
	assert(z.size() == size && "r.size is incorrect");

	vector<myType> local(size), remote(size);
	for (size_t i = 0; i < size; ++i)
	{
		myType prev = aes_prev->get64Bits();
		myType next = aes_next->get64Bits();
		local[i] = prev - next;
	}

	thread *threads = new thread[2];
	threads[0] = thread(sendVector<myType>, ref(local), prevParty(partyNum), size);
	threads[1] = thread(receiveVector<myType>, ref(remote), nextParty(partyNum), size);
	for (int i = 0; i < 2; i++)
		threads[i].join();
	delete[] threads;

	for (size_t i = 0; i < size; ++i)
		z[i] = std::make_pair(local[i], remote[i]);
}

void Precompute::getRandomBitShares(RSSVectorSmallType &a, size_t size)
{
	PreprocessingCostScope scope;
	assert(a.size() == size && "size mismatch for getRandomBitShares");
	for (size_t i = 0; i < size; ++i)
		a[i] = std::make_pair(aes_prev->get8Bits() & 1, aes_next->get8Bits() & 1);
}

void Precompute::getProductShares(const RSSVectorMyType &a, const RSSVectorMyType &b, RSSVectorMyType &c, size_t size)
{
	PreprocessingCostScope scope;
	assert(a.size() == size && "size mismatch for getProductShares");
	assert(b.size() == size && "size mismatch for getProductShares");
	assert(c.size() == size && "size mismatch for getProductShares");

	vector<myType> local(size), remote(size);
	for (size_t i = 0; i < size; ++i)
	{
		myType productShare = a[i].first * b[i].first
							+ a[i].first * b[i].second
							+ a[i].second * b[i].first;
		myType zeroShare = aes_prev->get64Bits() - aes_next->get64Bits();
		local[i] = productShare + zeroShare;
	}

	thread *threads = new thread[2];
	threads[0] = thread(sendVector<myType>, ref(local), prevParty(partyNum), size);
	threads[1] = thread(receiveVector<myType>, ref(remote), nextParty(partyNum), size);
	for (int i = 0; i < 2; i++)
		threads[i].join();
	delete[] threads;

	for (size_t i = 0; i < size; ++i)
		c[i] = std::make_pair(local[i], remote[i]);
}

void Precompute::getMatrixProductShares(const RSSVectorMyType &a, const RSSVectorMyType &b, RSSVectorMyType &c,
										size_t rows, size_t common_dim, size_t columns,
										size_t transpose_a, size_t transpose_b)
{
	PreprocessingCostScope scope;
	assert(a.size() == rows*common_dim && "size mismatch for getMatrixProductShares a");
	assert(b.size() == common_dim*columns && "size mismatch for getMatrixProductShares b");
	assert(c.size() == rows*columns && "size mismatch for getMatrixProductShares c");

	vector<myType> local(rows*columns), remote(rows*columns);
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			myType productShare = 0;
			for (size_t k = 0; k < common_dim; ++k)
			{
				const RSSMyType &left = transpose_a ? a[k*rows + i] : a[i*common_dim + k];
				const RSSMyType &right = transpose_b ? b[j*common_dim + k] : b[k*columns + j];
				productShare += left.first * right.first
							  + left.first * right.second
							  + left.second * right.first;
			}
			myType zeroShare = aes_prev->get64Bits() - aes_next->get64Bits();
			local[i*columns + j] = productShare + zeroShare;
		}
	}

	thread *threads = new thread[2];
	threads[0] = thread(sendVector<myType>, ref(local), prevParty(partyNum), rows*columns);
	threads[1] = thread(receiveVector<myType>, ref(remote), nextParty(partyNum), rows*columns);
	for (int i = 0; i < 2; i++)
		threads[i].join();
	delete[] threads;

	for (size_t i = 0; i < rows*columns; ++i)
		c[i] = std::make_pair(local[i], remote[i]);
}

void Precompute::getProductShares(const RSSVectorSmallType &a, const RSSVectorSmallType &b, RSSVectorSmallType &c, size_t size)
{
	PreprocessingCostScope scope;
	assert(a.size() == size && "size mismatch for getProductShares");
	assert(b.size() == size && "size mismatch for getProductShares");
	assert(c.size() == size && "size mismatch for getProductShares");

	vector<smallType> local(size), remote(size);
	for (size_t i = 0; i < size; ++i)
	{
		smallType productShare = additionModPrime[multiplicationModPrime[a[i].first][b[i].first]]
								[multiplicationModPrime[a[i].first][b[i].second]];
		productShare = additionModPrime[productShare][multiplicationModPrime[a[i].second][b[i].first]];
		smallType zeroShare = subtractModPrime[aes_prev->randModPrime()][aes_next->randModPrime()];
		local[i] = additionModPrime[productShare][zeroShare];
	}

	thread *threads = new thread[2];
	threads[0] = thread(sendVector<smallType>, ref(local), prevParty(partyNum), size);
	threads[1] = thread(receiveVector<smallType>, ref(remote), nextParty(partyNum), size);
	for (int i = 0; i < 2; i++)
		threads[i].join();
	delete[] threads;

	for (size_t i = 0; i < size; ++i)
		c[i] = std::make_pair(local[i], remote[i]);
}

void Precompute::getBitProductShares(const RSSVectorSmallType &a, const RSSVectorSmallType &b, RSSVectorSmallType &c, size_t size)
{
	PreprocessingCostScope scope;
	assert(a.size() == size && "size mismatch for getBitProductShares");
	assert(b.size() == size && "size mismatch for getBitProductShares");
	assert(c.size() == size && "size mismatch for getBitProductShares");

	vector<smallType> local(size), remote(size);
	for (size_t i = 0; i < size; ++i)
	{
		smallType productShare = (a[i].first & b[i].first)
								^ (a[i].first & b[i].second)
								^ (a[i].second & b[i].first);
		smallType zeroShare = (aes_prev->get8Bits() ^ aes_next->get8Bits()) & 1;
		local[i] = productShare ^ zeroShare;
	}

	thread *threads = new thread[2];
	threads[0] = thread(sendVector<smallType>, ref(local), prevParty(partyNum), size);
	threads[1] = thread(receiveVector<smallType>, ref(remote), nextParty(partyNum), size);
	for (int i = 0; i < 2; i++)
		threads[i].join();
	delete[] threads;

	for (size_t i = 0; i < size; ++i)
		c[i] = std::make_pair(local[i], remote[i]);
}

void Precompute::getBitToArithmeticShares(const RSSVectorSmallType &bits, RSSVectorMyType &arith, size_t size)
{
	PreprocessingCostScope scope;
	assert(bits.size() == size && "size mismatch for getBitToArithmeticShares");
	assert(arith.size() == size && "size mismatch for getBitToArithmeticShares");

	vector<myType> tripleShare(size, 0), additiveShare(size, 0);

	if (partyNum == PARTY_A)
	{
		vector<myType> maskForB(size), maskForC(size);
		for (size_t i = 0; i < size; ++i)
		{
			maskForB[i] = aes_indep->get64Bits();
			maskForC[i] = aes_indep->get64Bits();
			tripleShare[i] = maskForB[i] * maskForC[i];
		}
		sendVector<myType>(maskForB, PARTY_B, size);
		sendVector<myType>(maskForC, PARTY_C, size);
	}
	else if (partyNum == PARTY_B)
	{
		vector<myType> mask(size), xPlusMask(size), yPlusMask(size);
		receiveVector<myType>(mask, PARTY_A, size);
		for (size_t i = 0; i < size; ++i)
		{
			myType x = myType(bits[i].first & 1) * myType(bits[i].second & 1);
			xPlusMask[i] = x + mask[i];
		}
		sendVector<myType>(xPlusMask, PARTY_C, size);
		receiveVector<myType>(yPlusMask, PARTY_C, size);
		for (size_t i = 0; i < size; ++i)
		{
			myType x = myType(bits[i].first & 1) * myType(bits[i].second & 1);
			tripleShare[i] = x * yPlusMask[i];
		}
	}
	else
	{
		vector<myType> mask(size), xPlusMask(size), yPlusMask(size);
		receiveVector<myType>(mask, PARTY_A, size);
		receiveVector<myType>(xPlusMask, PARTY_B, size);
		for (size_t i = 0; i < size; ++i)
			yPlusMask[i] = myType(bits[i].second & 1) + mask[i];
		sendVector<myType>(yPlusMask, PARTY_B, size);
		for (size_t i = 0; i < size; ++i)
			tripleShare[i] = myType(0) - mask[i] * xPlusMask[i];
	}

	for (size_t i = 0; i < size; ++i)
	{
		myType first = myType(bits[i].first & 1);
		myType second = myType(bits[i].second & 1);
		additiveShare[i] = first - 2 * first * second + 4 * tripleShare[i];
	}

	reshareAdditiveRing(additiveShare, arith, size);
}

void Precompute::getBitToPrimeShares(const RSSVectorSmallType &bits, RSSVectorSmallType &arith, size_t size)
{
	PreprocessingCostScope scope;
	assert(bits.size() == size && "size mismatch for getBitToPrimeShares");
	assert(arith.size() == size && "size mismatch for getBitToPrimeShares");

	const smallType two = 2 % PRIME_NUMBER;
	const smallType four = 4 % PRIME_NUMBER;
	vector<smallType> tripleShare(size, 0), additiveShare(size, 0);

	if (partyNum == PARTY_A)
	{
		vector<smallType> maskForB(size), maskForC(size);
		for (size_t i = 0; i < size; ++i)
		{
			maskForB[i] = aes_indep->randModPrime();
			maskForC[i] = aes_indep->randModPrime();
			tripleShare[i] = multiplicationModPrime[maskForB[i]][maskForC[i]];
		}
		sendVector<smallType>(maskForB, PARTY_B, size);
		sendVector<smallType>(maskForC, PARTY_C, size);
	}
	else if (partyNum == PARTY_B)
	{
		vector<smallType> mask(size), xPlusMask(size), yPlusMask(size);
		receiveVector<smallType>(mask, PARTY_A, size);
		for (size_t i = 0; i < size; ++i)
		{
			smallType x = multiplicationModPrime[bits[i].first & 1][bits[i].second & 1];
			xPlusMask[i] = additionModPrime[x][mask[i]];
		}
		sendVector<smallType>(xPlusMask, PARTY_C, size);
		receiveVector<smallType>(yPlusMask, PARTY_C, size);
		for (size_t i = 0; i < size; ++i)
		{
			smallType x = multiplicationModPrime[bits[i].first & 1][bits[i].second & 1];
			tripleShare[i] = multiplicationModPrime[x][yPlusMask[i]];
		}
	}
	else
	{
		vector<smallType> mask(size), xPlusMask(size), yPlusMask(size);
		receiveVector<smallType>(mask, PARTY_A, size);
		receiveVector<smallType>(xPlusMask, PARTY_B, size);
		for (size_t i = 0; i < size; ++i)
			yPlusMask[i] = additionModPrime[bits[i].second & 1][mask[i]];
		sendVector<smallType>(yPlusMask, PARTY_B, size);
		for (size_t i = 0; i < size; ++i)
			tripleShare[i] = subtractModPrime[0][multiplicationModPrime[mask[i]][xPlusMask[i]]];
	}

	for (size_t i = 0; i < size; ++i)
	{
		smallType first = bits[i].first & 1;
		smallType second = bits[i].second & 1;
		smallType pairProduct = multiplicationModPrime[first][second];
		smallType negativePairs = multiplicationModPrime[two][pairProduct];
		smallType triple = multiplicationModPrime[four][tripleShare[i]];
		additiveShare[i] = additionModPrime[subtractModPrime[first][negativePairs]][triple];
	}

	reshareAdditivePrime(additiveShare, arith, size);
}

void Precompute::getTruncationMasks(RSSVectorMyType &r, RSSVectorMyType &rTrunc, size_t power, size_t size)
{
	PreprocessingCostScope scope;
	assert(r.size() == size && "size mismatch for getTruncationMasks r");
	assert(rTrunc.size() == size && "size mismatch for getTruncationMasks rTrunc");
	assert(power < BIT_SIZE && "getTruncationMasks power must be smaller than ring size");

	const size_t randomBitsPerValue = BIT_SIZE - 2;
	size_t bitSize = size * randomBitsPerValue;
	RSSVectorSmallType bits(bitSize);
	RSSVectorMyType bitShares(bitSize);
	getRandomBitShares(bits, bitSize);
	getBitToArithmeticShares(bits, bitShares, bitSize);

	for (size_t i = 0; i < size; ++i)
	{
		r[i] = make_pair(myType(0), myType(0));
		rTrunc[i] = make_pair(myType(0), myType(0));
		for (size_t localBit = 0; localBit < randomBitsPerValue; ++localBit)
		{
			size_t sourceBit = localBit + 2;
			size_t index = i*randomBitsPerValue + localBit;
			myType weight = myType(1) << (BIT_SIZE - 1 - sourceBit);
			myType truncatedWeight = weight >> power;

			r[i].first += bitShares[index].first * weight;
			r[i].second += bitShares[index].second * weight;
			if (truncatedWeight != 0)
			{
				rTrunc[i].first += bitShares[index].first * truncatedWeight;
				rTrunc[i].second += bitShares[index].second * truncatedWeight;
			}
		}
	}
}

void Precompute::getNonZeroPrimeJShares(MEVectorSmallType &gamma, size_t size)
{
	PreprocessingCostScope scope;
	assert(gamma.size() == size && "size mismatch for getNonZeroPrimeJShares");

	RSSVectorSmallType gammaSeed(size), gammaShares(size);
	getRandomShares(gammaSeed, size);
	getProductShares(gammaSeed, gammaSeed, gammaShares, size);
	for (size_t i = 0; i < size; ++i)
	{
		// PRIME_NUMBER is 3 mod 4, so delta^2 + 1 is always non-zero in F_p.
		if (partyNum == PARTY_A)
			gammaShares[i].first = additionModPrime[gammaShares[i].first][1];
		if (partyNum == PARTY_C)
			gammaShares[i].second = additionModPrime[gammaShares[i].second][1];
	}
	maskPrimeJFromRSS(*this, gammaShares, gamma, size);
}

void Precompute::getPreMSBObjects(RSSVectorMyType &lambda, MEVectorSmallType &lambdaBitsJ, MEVectorSmallType &lambdaMSB, MEVectorSmallType &sBits,
								  MEVectorSmallType &betaBit, MEVectorSmallType &betaPrime,
								  MEVectorSmallType &gamma, size_t size)
{
	PreprocessingCostScope scope;
	assert(lambda.size() == size && "size mismatch for getPreMSBObjects lambda");
	assert(lambdaBitsJ.size() == size*BIT_SIZE && "size mismatch for getPreMSBObjects lambdaBits");
	assert(lambdaMSB.size() == size && "size mismatch for getPreMSBObjects lambdaMSB");
	assert(sBits.size() == size*BIT_SIZE && "size mismatch for getPreMSBObjects sBits");
	assert(betaBit.size() == size && "size mismatch for getPreMSBObjects betaBit");
	assert(betaPrime.size() == size && "size mismatch for getPreMSBObjects betaPrime");
	assert(gamma.size() == size && "size mismatch for getPreMSBObjects gamma");

	size_t bitSize = size * BIT_SIZE;
	RSSVectorSmallType lambdaBits(bitSize), lambdaMSBBits(size), shiftedLambdaBits(bitSize), betaBits(size);
	RSSVectorMyType lambdaBitShares(bitSize);
	getRandomBitShares(lambdaBits, bitSize);
	getBitToArithmeticShares(lambdaBits, lambdaBitShares, bitSize);
	maskBooleanJFromRSS(*this, lambdaBits, lambdaBitsJ, bitSize);

	for (size_t i = 0; i < size; ++i)
	{
		lambda[i] = make_pair(myType(0), myType(0));
		for (size_t bit = 0; bit < BIT_SIZE; ++bit)
		{
			size_t index = i*BIT_SIZE + bit;
			myType weight = myType(1) << (BIT_SIZE - 1 - bit);
			lambda[i].first += lambdaBitShares[index].first * weight;
			lambda[i].second += lambdaBitShares[index].second * weight;
		}
		lambdaMSBBits[i] = lambdaBits[i*BIT_SIZE];
	}
	maskBooleanJFromRSS(*this, lambdaMSBBits, lambdaMSB, size);

	for (size_t i = 0; i < size; ++i)
	{
		for (size_t bit = 0; bit + 1 < BIT_SIZE; ++bit)
			shiftedLambdaBits[i*BIT_SIZE + bit] = lambdaBits[i*BIT_SIZE + bit + 1];
		shiftedLambdaBits[i*BIT_SIZE + BIT_SIZE - 1] = make_pair(smallType(0), smallType(0));
	}

	RSSVectorSmallType shiftedLambdaPrime(bitSize);
	getBitToPrimeShares(shiftedLambdaBits, shiftedLambdaPrime, bitSize);
	maskPrimeJFromRSS(*this, shiftedLambdaPrime, sBits, bitSize);

	getRandomBitShares(betaBits, size);
	maskBooleanJFromRSS(*this, betaBits, betaBit, size);
	RSSVectorSmallType betaPrimeShares(size);
	getBitToPrimeShares(betaBits, betaPrimeShares, size);
	maskPrimeJFromRSS(*this, betaPrimeShares, betaPrime, size);

	getNonZeroPrimeJShares(gamma, size);
}


//m_0 is random shares of 0, m_1 is random shares of 1 in RSSMyType. 
//This function generates random bits c and corresponding RSSMyType values m_c
void Precompute::getSelectorBitShares(RSSVectorSmallType &c, RSSVectorMyType &m_c, size_t size)
{
	PreprocessingCostScope scope;
	assert(c.size() == size && "size mismatch for getSelectorBitShares");
	assert(m_c.size() == size && "size mismatch for getSelectorBitShares");
	for(auto &it : c)
		it = std::make_pair(0,0);

	for(auto &it : m_c)
		it = std::make_pair(0,0);
}

//Shares of random r, shares of bits of that, and shares of wrap3 of that.
void Precompute::getShareConvertObjects(RSSVectorMyType &r, RSSVectorSmallType &shares_r, 
										RSSVectorSmallType &alpha, size_t size)
{
	PreprocessingCostScope scope;
	assert(shares_r.size() == size*BIT_SIZE && "getShareConvertObjects size mismatch");
	for(auto &it : r)
		it = std::make_pair(0,0);

	for(auto &it : shares_r)
		it = std::make_pair(0,0);

	for(auto &it : alpha)
		it = std::make_pair(0,0);
}
