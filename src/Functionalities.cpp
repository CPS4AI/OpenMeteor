
#pragma once
#include "Functionalities.h"
#include "Precompute.h"
#include <thread>


using namespace std;
extern Precompute PrecomputeObject;

/******************************** Functionalities 2PC ********************************/
// Share Truncation, truncate shares of a by power (in place) (power is logarithmic)

void funcReconstructBit(const RSSVectorSmallType &a, vector<smallType> &b, size_t size, string str, bool print)
{
	log_print("Reconst: RSSSmallType (bits), smallType (bit)");


		vector<smallType> a_next(size), a_prev(size);
		for (int i = 0; i < size; ++i)
		{
			a_prev[i] = 0;
			a_next[i] = a[i].first;
			b[i] = a[i].first;
			b[i] = b[i] ^ a[i].second;
		}

		thread *threads = new thread[2];
		threads[0] = thread(sendVector<smallType>, ref(a_next), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<smallType>, ref(a_prev), prevParty(partyNum), size);
		for (int i = 0; i < 2; i++)
			threads[i].join();
		delete[] threads;

		for (int i = 0; i < size; ++i)
			b[i] = b[i] ^ a_prev[i];

		if (print)
		{
			std::cout << str << ": \t\t";
			for (int i = 0; i < size; ++i)
				cout << (int)(b[i]) << " ";
			std::cout << std::endl;
		}

}

void funcReconstruct(const RSSVectorSmallType &a, vector<smallType> &b, size_t size, string str, bool print)
{
	log_print("Reconst: RSSSmallType, smallType");


		vector<smallType> a_next(size), a_prev(size);
		for (int i = 0; i < size; ++i)
		{
			a_prev[i] = 0;
			a_next[i] = a[i].first;
			b[i] = a[i].first;
			b[i] = additionModPrime[b[i]][a[i].second];
		}

		thread *threads = new thread[2];

		threads[0] = thread(sendVector<smallType>, ref(a_next), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<smallType>, ref(a_prev), prevParty(partyNum), size);

		for (int i = 0; i < 2; i++)
			threads[i].join();

		delete[] threads;

		for (int i = 0; i < size; ++i)
			b[i] = additionModPrime[b[i]][a_prev[i]];

		if (print)
		{
			std::cout << str << ": \t\t";
			for (int i = 0; i < size; ++i)
				cout << (int)(b[i]) << " ";
			std::cout << std::endl;
		}

}

void funcReconstruct(const RSSVectorMyType &a, vector<myType> &b, size_t size, string str, bool print)
{
	log_print("Reconst: RSSMyType, myType");
	assert(a.size() == size && "a.size mismatch for reconstruct function");

	vector<myType> a_next(size), a_prev(size, 0);
	for (int i = 0; i < size; ++i)
	{
		a_next[i] = a[i].first;
		b[i] = a[i].first + a[i].second;
	}

		thread *threads = new thread[2];

		threads[0] = thread(sendVector<myType>, ref(a_next), nextParty(partyNum), size);
		threads[1] = thread(receiveVector<myType>, ref(a_prev), prevParty(partyNum), size);

		for (int i = 0; i < 2; i++)
			threads[i].join();

		delete[] threads;

		for (int i = 0; i < size; ++i)
			b[i] = b[i] + a_prev[i];

		if (print)
		{
			std::cout << str << ": \t\t";
			for (int i = 0; i < size; ++i)
				print_linear(b[i], "SIGNED");
			std::cout << std::endl;
		}

}

namespace
{
	myType signedRightShift(myType value, size_t power)
	{
		assert(power < BIT_SIZE && "shift power must be smaller than ring size");
		if (power == 0)
			return value;

		myType shifted = value >> power;
		if (getMSB(value) == 0)
			return shifted;

		myType signMask = (~myType(0)) << (BIT_SIZE - power);
		return shifted | signMask;
	}

	RSSMyType negateRSS(const RSSMyType &share)
	{
		return make_pair(myType(0) - share.first, myType(0) - share.second);
	}
}

void Meteor_funcTruncate(const MEVectorType &a, MEVectorType &b, size_t size, size_t power)
{
	log_print("Meteor faithful truncation");
	assert(a.size() == size && "input a is size error for truncation");
	assert(b.size() == size && "output b is size error for truncation");
	assert(power < BIT_SIZE && "truncation power must be smaller than ring size");

	if (power == 0)
	{
		for (size_t i = 0; i < size; ++i)
			b[i] = a[i];
		return;
	}

	RSSVectorMyType r(size), rTrunc(size), openedShares(size);
	vector<myType> opened(size, 0);
	PrecomputeObject.getTruncationMasks(r, rTrunc, power, size);

	for (size_t i = 0; i < size; ++i)
	{
		openedShares[i] = a[i].second + r[i];
		if (partyNum == PARTY_A)
			openedShares[i].first += a[i].first;
		if (partyNum == PARTY_C)
			openedShares[i].second += a[i].first;
	}

	funcReconstruct(openedShares, opened, size, "faithful truncation masked value", false);

	for (size_t i = 0; i < size; ++i)
		b[i] = make_pair(signedRightShift(opened[i], power), negateRSS(rTrunc[i]));
}


/******************************** Functionalities MPC ********************************/
// Matrix Multiplication of a*b = c with transpose flags for a,b.
// Output is a share between PARTY_A and PARTY_B.
// a^transpose_a is rows*common_dim and b^transpose_b is common_dim*columns
// define Meteor MatMul
void Meteor_funcMatMul(const MEVectorType &a, const MEVectorType &b, MEVectorType &c,
					size_t rows, size_t common_dim, size_t columns,
				 	size_t transpose_a, size_t transpose_b, size_t truncation)
{
	log_print("Meteor_funcMatMul");
	cout << "WE ARE USING METEOR MATMUL!" << endl;
	assert(a.size() == rows*common_dim && "Matrix a incorrect for Mat-Mul");
	assert(b.size() == common_dim*columns && "Matrix b incorrect for Mat-Mul");
	assert(c.size() == rows*columns && "Matrix c incorrect for Mat-Mul");

	size_t final_size = rows*columns;
	RSSVectorMyType temp3(final_size), productMask(final_size), outputMask(final_size);
	RSSVectorMyType aMask(rows*common_dim), bMask(common_dim*columns);
	vector<myType> Delta(final_size, 0);

	for (size_t i = 0; i < rows*common_dim; ++i)
		aMask[i] = a[i].second;
	for (size_t i = 0; i < common_dim*columns; ++i)
		bMask[i] = b[i].second;
	PrecomputeObject.getMatrixProductShares(aMask, bMask, productMask, rows, common_dim, columns, transpose_a, transpose_b);
	PrecomputeObject.getRandomShares(outputMask, final_size);

	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			size_t out = i*columns + j;
			RSSMyType corr = productMask[out] - outputMask[out];
			temp3[out] = corr;
			for (size_t k = 0; k < common_dim; ++k)
			{
				const METype &left = transpose_a ? a[k*rows + i] : a[i*common_dim + k];
				const METype &right = transpose_b ? b[j*common_dim + k] : b[k*columns + j];
				temp3[out].first += left.first * right.second.first + left.second.first * right.first;
				temp3[out].second += left.first * right.second.second + left.second.second * right.first;

				if (partyNum == PARTY_A)
					temp3[out].first += left.first * right.first;
				if (partyNum == PARTY_C)
					temp3[out].second += left.first * right.first;
			}
		}
	}

	funcReconstruct(temp3, Delta, final_size, "Delta Reconst", false);

	for (int i = 0; i < final_size; ++i)
		c[i] = make_pair(Delta[i], outputMask[i]);

	if (truncation)
		Meteor_funcTruncate(c, c, final_size, truncation);
}

// 2-input dot product, used in private compare and relu.
void Meteor_funcDotProduct(const MEVectorType &a, const MEVectorType &b, MEVectorType &c, size_t size, bool truncation, size_t precision){

	log_print("Meteor Function for Dot Product");
	assert(a.size() == size && "a size is incorrect");
	assert(b.size() == size && "b size is incorrect");
	assert(c.size() == size && "c size is incorrect");

	RSSVectorMyType temp3(size), productMask(size), outputMask(size);
	vector<myType> Delta(size, 0), temp4(size, 0);
	RSSVectorMyType aMask(size), bMask(size);
	for (size_t i = 0; i < size; ++i)
	{
		aMask[i] = a[i].second;
		bMask[i] = b[i].second;
	}
	PrecomputeObject.getProductShares(aMask, bMask, productMask, size);
	PrecomputeObject.getRandomShares(outputMask, size);

	for (int i = 0; i < size; i++){
		RSSMyType corr = productMask[i] - outputMask[i];
		if(partyNum == PARTY_A){
			temp3[i].first = a[i].first * b[i].first + a[i].first * b[i].second.first + a[i].second.first * b[i].first + corr.first;
			temp3[i].second = a[i].first * b[i].second.second + a[i].second.second * b[i].first + corr.second;
		}
		else if(partyNum == PARTY_B){
			temp3[i].first = a[i].first * b[i].second.first + a[i].second.first * b[i].first + corr.first;
			temp3[i].second = a[i].first * b[i].second.second + a[i].second.second * b[i].first + corr.second;
		}
		else{
			temp3[i].first = a[i].first * b[i].second.first + a[i].second.first * b[i].first + corr.first;
			temp3[i].second = a[i].first * b[i].first + a[i].first * b[i].second.second + a[i].second.second * b[i].first + corr.second;
		}
		temp4[i] = temp3[i].second;

	}
	//funcReconstruct(temp3, Delta, size, "Delta Reconst", false);


	thread *threads = new thread[2];

	threads[0] = thread(sendVector<myType>, ref(temp4), prevParty(partyNum), size);
	threads[1] = thread(receiveVector<myType>, ref(Delta), nextParty(partyNum), size);

	for (int i = 0; i < 2; i++)
		threads[i].join();
	delete[] threads;

	for(int i = 0; i < size; i++){
		Delta[i] = Delta[i] + temp3[i].first + temp3[i].second;
	}

	for (int i = 0; i < size; ++i)
	{
		c[i] = make_pair(Delta[i], outputMask[i]);
	}

	if (truncation)
		Meteor_funcTruncate(c, c, size, precision);

}

void funcReconstruct3out3parallel(const vector<smallType> &a, vector<smallType> &b, size_t size, string str, bool print)
{
 	assert(a.size() == size && "a.size mismatch for reconstruct function");

	vector<smallType> a_next(size), a_prev(size);
 	for (int i = 0; i < size; ++i)
 	{
 		a_prev[i] = 0;
 		a_next[i] = 0;
 	}

 	thread *threads = new thread[4];

 	threads[0] = thread(sendVector<smallType>, ref(a), nextParty(partyNum), size);
 	threads[1] = thread(sendVector<smallType>, ref(a), prevParty(partyNum), size);
 	threads[2] = thread(receiveVector<smallType>, ref(a_next), nextParty(partyNum), size);
 	threads[3] = thread(receiveVector<smallType>, ref(a_prev), prevParty(partyNum), size);

 	for (int i = 0; i < 4; i++)
 		threads[i].join();

 	delete[] threads;

 	for (int i = 0; i < size; ++i)
 		b[i] = additionModPrime[additionModPrime[a[i]][a_prev[i]]][a_next[i]];
 }

void funcReconstruct3out3(const vector<smallType> &a, vector<smallType> &b, size_t size, string str, bool print)
{
	log_print("Reconst: smallType, smallType");
	assert(a.size() == size && "a.size mismatch for reconstruct function");

	vector<smallType> temp_A(size,0), temp_B(size, 0);

	if (partyNum == PARTY_A or partyNum == PARTY_B)
		sendVector<smallType>(a, PARTY_C, size);

	if (partyNum == PARTY_C)
	{
		receiveVector<smallType>(temp_A, PARTY_A, size);
		receiveVector<smallType>(temp_B, PARTY_B, size);
		addVectors<smallType>(temp_A, a, temp_A, size);
		addVectors<smallType>(temp_B, temp_A, b, size);
		sendVector<smallType>(b, PARTY_A, size);
		sendVector<smallType>(b, PARTY_B, size);
	}

	if (partyNum == PARTY_A or partyNum == PARTY_B)
		receiveVector<smallType>(b, PARTY_C, size);

	if (print)
	{
		std::cout << str << ": \t\t";
		for (int i = 0; i < size; ++i)
			print_linear(b[i], "SIGNED");
		std::cout << std::endl;
	}
}


void funcReconstruct3out3parallel(const vector<myType> &a, vector<myType> &b, size_t size, string str, bool print)
{
 	assert(a.size() == size && "a.size mismatch for reconstruct function");

	vector<myType> a_next(size), a_prev(size);
 	for (int i = 0; i < size; ++i)
 	{
 		a_prev[i] = 0;
 		a_next[i] = 0;
 	}

 	thread *threads = new thread[4];

 	threads[0] = thread(sendVector<myType>, ref(a), nextParty(partyNum), size);
 	threads[1] = thread(sendVector<myType>, ref(a), prevParty(partyNum), size);
 	threads[2] = thread(receiveVector<myType>, ref(a_next), nextParty(partyNum), size);
 	threads[3] = thread(receiveVector<myType>, ref(a_prev), prevParty(partyNum), size);

 	for (int i = 0; i < 4; i++)
 		threads[i].join();

 	delete[] threads;

 	for (int i = 0; i < size; ++i)
 		b[i] = a[i] + a_prev[i] + a_next[i];
 }

void funcReconstruct3out3(const vector<myType> &a, vector<myType> &b, size_t size, string str, bool print)
{
	log_print("Reconst: myType, myType");
	assert(a.size() == size && "a.size mismatch for reconstruct function");

	vector<myType> temp_A(size,0), temp_B(size, 0);

	if (partyNum == PARTY_A or partyNum == PARTY_B)
		sendVector<myType>(a, PARTY_C, size);

	if (partyNum == PARTY_C)
	{
		receiveVector<myType>(temp_A, PARTY_A, size);
		receiveVector<myType>(temp_B, PARTY_B, size);
		addVectors<myType>(temp_A, a, temp_A, size);
		addVectors<myType>(temp_B, temp_A, b, size);
		sendVector<myType>(b, PARTY_A, size);
		sendVector<myType>(b, PARTY_B, size);
	}

	if (partyNum == PARTY_A or partyNum == PARTY_B)
		receiveVector<myType>(b, PARTY_C, size);

	if (print)
	{
		std::cout << str << ": \t\t";
		for (int i = 0; i < size; ++i)
			print_linear(b[i], "SIGNED");
		std::cout << std::endl;
	}
}


void Meteor_funcDotProduct(const MEVectorSmallType &a, const MEVectorSmallType &b, MEVectorSmallType &c, size_t size)
{
	log_print("2-input meteor dot protudct in prime field");
	assert(a.size() == size && "Vector a incorrect for vec-prod");
	assert(b.size() == size && "vector b incorrect for vec-prod");
	assert(c.size() == size && "vector c incorrect for vec-prod");

	RSSVectorSmallType temp(size), productMask(size), outputMask(size);
	vector<smallType> temp4(size, 0);
	RSSVectorSmallType aMask(size), bMask(size);
	for (size_t i = 0; i < size; ++i)
	{
		aMask[i] = a[i].second;
		bMask[i] = b[i].second;
	}
	PrecomputeObject.getProductShares(aMask, bMask, productMask, size);
	PrecomputeObject.getRandomShares(outputMask, size);

	for(int i = 0; i < size; i++){
		RSSSmallType corr;
		corr.first = subtractModPrime[productMask[i].first][outputMask[i].first];
		corr.second = subtractModPrime[productMask[i].second][outputMask[i].second];

		temp[i].first = multiplicationModPrime[a[i].first][b[i].second.first];
		temp[i].first = additionModPrime[temp[i].first][multiplicationModPrime[a[i].second.first][b[i].first]];
		temp[i].first = additionModPrime[temp[i].first][corr.first];

		temp[i].second = multiplicationModPrime[a[i].first][b[i].second.second];
		temp[i].second = additionModPrime[temp[i].second][multiplicationModPrime[a[i].second.second][b[i].first]];
		temp[i].second = additionModPrime[temp[i].second][corr.second];

		if(partyNum == PARTY_A)
			temp[i].first = additionModPrime[temp[i].first][multiplicationModPrime[a[i].first][b[i].first]];
		if(partyNum == PARTY_C)
			temp[i].second = additionModPrime[temp[i].second][multiplicationModPrime[a[i].first][b[i].first]];
	}

	vector<smallType> Delta(size, 0);
	funcReconstruct(temp, Delta, size, "2-input small dot-prod Reconst", false);


	for(int i =0; i < size; i++)
	{
		c[i] = make_pair(Delta[i], outputMask[i]);

	}
}



// Term by term multiplication boolean shares
void Meteor_funcDotProductBits(const MEVectorSmallType &a, const MEVectorSmallType &b, MEVectorSmallType &c, size_t size)
{
	log_print("Meteor funcDotProductBits");
	assert(a.size() == size && "Vector a incorrect for DPB");
	assert(b.size() == size && "Vector b incorrect for DPB");
	assert(c.size() == size && "Vector c incorrect for DPB");

	RSSVectorSmallType temp(size), productMask(size), outputMask(size);
	RSSVectorSmallType aMask(size), bMask(size);
	for (size_t i = 0; i < size; ++i)
	{
		aMask[i] = a[i].second;
		bMask[i] = b[i].second;
	}
	PrecomputeObject.getBitProductShares(aMask, bMask, productMask, size);
	PrecomputeObject.getRandomBitShares(outputMask, size);

	for (int i = 0; i < size; i++)
	{
		RSSSmallType corr = productMask[i] ^ outputMask[i];
		temp[i].first = (a[i].first and b[i].second.first) ^ (a[i].second.first and b[i].first) ^ corr.first;
		temp[i].second = (a[i].first and b[i].second.second) ^ (a[i].second.second and b[i].first) ^ corr.second;

		if(partyNum == PARTY_A)
		{
			temp[i].first = temp[i].first ^ (a[i].first and b[i].first);
		}
		if(partyNum == PARTY_C)
		{
			temp[i].second = temp[i].second ^ (a[i].first and b[i].first);
		}
	}
	vector<smallType> Delta(size, 0);
	funcReconstructBit(temp, Delta, size, "2-input small bit-prod Reconst", false);

	for(int i =0; i < size; i++)
	{
		c[i] = make_pair(Delta[i], outputMask[i]);
	}
}

void Meteor_funcBit2A(const MEVectorSmallType &bits, MEVectorType &arith, size_t size)
{
	log_print("Meteor funcBit2A");
	assert(bits.size() == size && "Vector bits incorrect for Bit2A");
	assert(arith.size() == size && "Vector arith incorrect for Bit2A");

	RSSVectorSmallType bitMask(size);
	RSSVectorMyType lambdaArith(size), outputMask(size), deltaShares(size);
	vector<myType> Delta(size, 0);
	for (size_t i = 0; i < size; ++i)
	{
		assert((bits[i].first & 1) == bits[i].first && "Bit2A expects public deltas to be bits");
		bitMask[i] = bits[i].second;
	}

	PrecomputeObject.getBitToArithmeticShares(bitMask, lambdaArith, size);
	PrecomputeObject.getRandomShares(outputMask, size);

	for (size_t i = 0; i < size; ++i)
	{
		myType publicBit = myType(bits[i].first & 1);
		myType scale = myType(1) - 2 * publicBit;
		deltaShares[i].first = scale * lambdaArith[i].first - outputMask[i].first;
		deltaShares[i].second = scale * lambdaArith[i].second - outputMask[i].second;

		if (publicBit == 1)
		{
			if (partyNum == PARTY_A)
				deltaShares[i].first += 1;
			if (partyNum == PARTY_C)
				deltaShares[i].second += 1;
		}
	}

	funcReconstruct(deltaShares, Delta, size, "Meteor Bit2A Reconst", false);
	for (size_t i = 0; i < size; ++i)
		arith[i] = make_pair(Delta[i], outputMask[i]);
}

void Thunder_funcNMultGateOnline(const MEVectorType &c_1, MEVectorType &c_2, size_t size)
{
	log_print("only support 64 bit!");
	assert((size % 4) == 0 && "size is not in {64, 16, 4}, error!");
	size_t res_size = size / 4;
	RSSVectorMyType l0(res_size), l1(res_size), l2(res_size), l3(res_size);
	RSSVectorMyType p01(res_size), p02(res_size), p03(res_size), p12(res_size), p13(res_size), p23(res_size);
	RSSVectorMyType t012(res_size), t013(res_size), t023(res_size), t123(res_size);
	RSSVectorMyType q0123(res_size), outputMask(res_size), temp(res_size);
	vector<myType> Delta(res_size, 0);

	for (size_t i = 0; i < res_size; ++i)
	{
		l0[i] = c_1[4*i].second;
		l1[i] = c_1[4*i + 1].second;
		l2[i] = c_1[4*i + 2].second;
		l3[i] = c_1[4*i + 3].second;
	}

	PrecomputeObject.getProductShares(l0, l1, p01, res_size);
	PrecomputeObject.getProductShares(l0, l2, p02, res_size);
	PrecomputeObject.getProductShares(l0, l3, p03, res_size);
	PrecomputeObject.getProductShares(l1, l2, p12, res_size);
	PrecomputeObject.getProductShares(l1, l3, p13, res_size);
	PrecomputeObject.getProductShares(l2, l3, p23, res_size);
	PrecomputeObject.getProductShares(p01, l2, t012, res_size);
	PrecomputeObject.getProductShares(p01, l3, t013, res_size);
	PrecomputeObject.getProductShares(p02, l3, t023, res_size);
	PrecomputeObject.getProductShares(p12, l3, t123, res_size);
	PrecomputeObject.getProductShares(p01, p23, q0123, res_size);
	PrecomputeObject.getRandomShares(outputMask, res_size);

	for (size_t i = 0; i < res_size; ++i)
	{
		myType m0 = c_1[4*i].first;
		myType m1 = c_1[4*i + 1].first;
		myType m2 = c_1[4*i + 2].first;
		myType m3 = c_1[4*i + 3].first;

		auto addScaled = [](RSSMyType &acc, const RSSMyType &share, myType scalar) {
			acc.first += share.first * scalar;
			acc.second += share.second * scalar;
		};

		temp[i] = q0123[i] - outputMask[i];
		addScaled(temp[i], l0[i], m1*m2*m3);
		addScaled(temp[i], l1[i], m0*m2*m3);
		addScaled(temp[i], l2[i], m0*m1*m3);
		addScaled(temp[i], l3[i], m0*m1*m2);
		addScaled(temp[i], p01[i], m2*m3);
		addScaled(temp[i], p02[i], m1*m3);
		addScaled(temp[i], p03[i], m1*m2);
		addScaled(temp[i], p12[i], m0*m3);
		addScaled(temp[i], p13[i], m0*m2);
		addScaled(temp[i], p23[i], m0*m1);
		addScaled(temp[i], t012[i], m3);
		addScaled(temp[i], t013[i], m2);
		addScaled(temp[i], t023[i], m1);
		addScaled(temp[i], t123[i], m0);

		myType publicProduct = m0*m1*m2*m3;
		if (partyNum == PARTY_A)
			temp[i].first += publicProduct;
		if (partyNum == PARTY_C)
			temp[i].second += publicProduct;
	}

	funcReconstruct(temp, Delta, res_size, "4-input integer n-multiply Reconst", false);

	for(size_t i =0; i < res_size; i++)
		c_2[i] = make_pair(Delta[i], outputMask[i]);

}

void Thunder_funcNMultGate(const MEVectorType &c, size_t size, size_t N)
{
	size_t sizeLong = size * N;
	size_t rounds = size_t(log2(N)/2);

	//start_m();
	//start_communication();
	//start_time();
	// Setup for our improved 4-input mult gates over meteor, start
	size_t setup_size = N + N / 4 + N / 16;
	MEVectorType a(setup_size / 2,  make_pair(1, make_pair(0,0))),  b(setup_size / 2,  make_pair(2, make_pair(0,0))), ab(setup_size / 2);
	Meteor_funcDotProduct(a, b, ab, setup_size/2, false, 0);
	// end;

	//end_m("setup for 4-input gates for N fan-in");
	//end_time("setup for 4-mult with N fan-in");
	//pause_communication();
	//end_communication("setup for 4-mult with N fan-in");

   //start_communication();
   //start_time();
	vector<MEVectorType> c_tau(rounds);
	for(int i =0; i < rounds; i++){
		c_tau[i] = MEVectorType(sizeLong >> (2*(i+1)), make_pair(0, make_pair(0,0)));
	}
	Thunder_funcNMultGateOnline(c, c_tau[0], sizeLong);
	for(int i = 1; i < rounds; i++){
		Thunder_funcNMultGateOnline(c_tau[i-1], c_tau[i], sizeLong >> (2*i));
	}
	//end_time("online for 4-mult with N fan-in");
	//pause_communication();
	//end_communication("online for 4-mult with N fan-in");

}




void Meteor_funcMultiplyNeighbors(const MEVectorSmallType &c_1, MEVectorSmallType &c_2, size_t size)
{
	log_print("only support 64 bit!");
	assert((size % 4) == 0 && "size is not in {64, 16, 4}, error!");
	size_t res_size = size / 4;
	RSSVectorSmallType l0(res_size), l1(res_size), l2(res_size), l3(res_size);
	RSSVectorSmallType p01(res_size), p02(res_size), p03(res_size), p12(res_size), p13(res_size), p23(res_size);
	RSSVectorSmallType t012(res_size), t013(res_size), t023(res_size), t123(res_size);
	RSSVectorSmallType q0123(res_size), outputMask(res_size), temp(res_size);
	vector<smallType> Delta(res_size, 0);

	for (size_t i = 0; i < res_size; ++i)
	{
		l0[i] = c_1[4*i].second;
		l1[i] = c_1[4*i + 1].second;
		l2[i] = c_1[4*i + 2].second;
		l3[i] = c_1[4*i + 3].second;
	}

	PrecomputeObject.getProductShares(l0, l1, p01, res_size);
	PrecomputeObject.getProductShares(l0, l2, p02, res_size);
	PrecomputeObject.getProductShares(l0, l3, p03, res_size);
	PrecomputeObject.getProductShares(l1, l2, p12, res_size);
	PrecomputeObject.getProductShares(l1, l3, p13, res_size);
	PrecomputeObject.getProductShares(l2, l3, p23, res_size);
	PrecomputeObject.getProductShares(p01, l2, t012, res_size);
	PrecomputeObject.getProductShares(p01, l3, t013, res_size);
	PrecomputeObject.getProductShares(p02, l3, t023, res_size);
	PrecomputeObject.getProductShares(p12, l3, t123, res_size);
	PrecomputeObject.getProductShares(p01, p23, q0123, res_size);
	PrecomputeObject.getRandomShares(outputMask, res_size);

	for (size_t i = 0; i < res_size; ++i)
	{
		smallType m0 = c_1[4*i].first;
		smallType m1 = c_1[4*i + 1].first;
		smallType m2 = c_1[4*i + 2].first;
		smallType m3 = c_1[4*i + 3].first;

		auto mul2 = [](smallType x, smallType y) { return multiplicationModPrime[x][y]; };
		auto mul3 = [&](smallType x, smallType y, smallType z) { return mul2(mul2(x, y), z); };
		auto addScaled = [&](RSSSmallType &acc, const RSSSmallType &share, smallType scalar) {
			acc.first = additionModPrime[acc.first][multiplicationModPrime[share.first][scalar]];
			acc.second = additionModPrime[acc.second][multiplicationModPrime[share.second][scalar]];
		};

		temp[i].first = subtractModPrime[q0123[i].first][outputMask[i].first];
		temp[i].second = subtractModPrime[q0123[i].second][outputMask[i].second];

		addScaled(temp[i], l0[i], mul3(m1, m2, m3));
		addScaled(temp[i], l1[i], mul3(m0, m2, m3));
		addScaled(temp[i], l2[i], mul3(m0, m1, m3));
		addScaled(temp[i], l3[i], mul3(m0, m1, m2));
		addScaled(temp[i], p01[i], mul2(m2, m3));
		addScaled(temp[i], p02[i], mul2(m1, m3));
		addScaled(temp[i], p03[i], mul2(m1, m2));
		addScaled(temp[i], p12[i], mul2(m0, m3));
		addScaled(temp[i], p13[i], mul2(m0, m2));
		addScaled(temp[i], p23[i], mul2(m0, m1));
		addScaled(temp[i], t012[i], m3);
		addScaled(temp[i], t013[i], m2);
		addScaled(temp[i], t023[i], m1);
		addScaled(temp[i], t123[i], m0);

		smallType publicProduct = mul2(mul2(m0, m1), mul2(m2, m3));
		if (partyNum == PARTY_A)
			temp[i].first = additionModPrime[temp[i].first][publicProduct];
		if (partyNum == PARTY_C)
			temp[i].second = additionModPrime[temp[i].second][publicProduct];
	}

	funcReconstruct(temp, Delta, res_size, "4-input small 4-neighbormultiply Reconst", false);

	for(int i =0; i < res_size; i++)
	{
		c_2[i] = make_pair(Delta[i], outputMask[i]);

	}

}

void Meteor_funcCrunchMultiply(const MEVectorSmallType &c, vector<smallType> &betaPrime, size_t size)
{
	MEVectorSmallType gamma(size);
	PrecomputeObject.getNonZeroPrimeJShares(gamma, size);
	Meteor_funcCrunchMultiply(c, gamma, betaPrime, size);
}

void Meteor_funcCrunchMultiply(const MEVectorSmallType &c, const MEVectorSmallType &gamma, vector<smallType> &betaPrime, size_t size)
{
	size_t sizeLong = size * BIT_SIZE;
	size_t rounds = size_t(log2(BIT_SIZE)/2);
	assert(c.size() == sizeLong && "CrunchMultiply input size mismatch");
	assert(gamma.size() == size && "CrunchMultiply gamma size mismatch");
	assert(betaPrime.size() == size && "CrunchMultiply output size mismatch");

	vector<MEVectorSmallType> c_tau(rounds);
	for(int i =0; i < rounds; i++){
		c_tau[i] = MEVectorSmallType(sizeLong >> (2*(i+1)), make_pair(0, make_pair(0,0)));
	}
	Meteor_funcMultiplyNeighbors(c, c_tau[0], sizeLong);
	for(int i = 1; i < rounds; i++){
		Meteor_funcMultiplyNeighbors(c_tau[i-1], c_tau[i], sizeLong >> (2*i));
	}

	MEVectorSmallType maskedProduct(size);
	RSSVectorSmallType maskedProductMask(size);
	vector<smallType> openedMask(size, 0);

	Meteor_funcDotProduct(c_tau[rounds-1], gamma, maskedProduct, size);

	for(int i=0; i < size; i++){
		maskedProductMask[i] = maskedProduct[i].second;
	}

	funcReconstruct(maskedProductMask, openedMask, size, "masked k-mult product mask", false);

	for(int i = 0; i < size; i++){
		smallType openedProduct = additionModPrime[maskedProduct[i].first][openedMask[i]];
		betaPrime[i] = (openedProduct != 0);
	}
}

void Meteor_funcSecMSB(const MEVectorType &x, const MEVectorSmallType &lambdaBits, const MEVectorSmallType &lambdaMSB, const MEVectorSmallType &sBits,
					   const MEVectorSmallType &betaBit, const MEVectorSmallType &betaPrime,
					   const MEVectorSmallType &gamma, MEVectorSmallType &msb, size_t size)
{
	assert(x.size() == size && "SecMSB input size mismatch");
	assert(lambdaBits.size() == size*BIT_SIZE && "SecMSB lambdaBits size mismatch");
	assert(lambdaMSB.size() == size && "SecMSB lambdaMSB size mismatch");
	assert(sBits.size() == size*BIT_SIZE && "SecMSB sBits size mismatch");
	assert(betaBit.size() == size && "SecMSB betaBit size mismatch");
	assert(betaPrime.size() == size && "SecMSB betaPrime size mismatch");
	assert(gamma.size() == size && "SecMSB gamma size mismatch");
	assert(msb.size() == size && "SecMSB output size mismatch");

	MEVectorSmallType carry(size, make_pair(smallType(0), make_pair(smallType(0), smallType(0))));
	MEVectorSmallType propagate(size), generate(size), carryProduct(size);

	auto notJBit = [](const MESmallType &bit) {
		MESmallType ret = bit;
		ret.first ^= 1;
		return ret;
	};
	auto xorJBit = [](const MESmallType &left, const MESmallType &right) {
		MESmallType ret;
		ret.first = left.first ^ right.first;
		ret.second = left.second ^ right.second;
		return ret;
	};

	for (size_t bit = 0; bit + 1 < BIT_SIZE; ++bit)
	{
		for (size_t i = 0; i < size; ++i)
		{
			size_t index = i*BIT_SIZE + (BIT_SIZE - 1 - bit);
			smallType mxBit = smallType((x[i].first >> bit) & 1);
			if (mxBit == 0)
			{
				generate[i] = make_pair(smallType(0), make_pair(smallType(0), smallType(0)));
				propagate[i] = lambdaBits[index];
			}
			else
			{
				generate[i] = lambdaBits[index];
				propagate[i] = notJBit(lambdaBits[index]);
			}
		}

		Meteor_funcDotProductBits(propagate, carry, carryProduct, size);
		for (size_t i = 0; i < size; ++i)
			carry[i] = xorJBit(generate[i], carryProduct[i]);
	}

	for (size_t i = 0; i < size; ++i)
	{
		msb[i].first = getMSB(x[i].first) ^ lambdaMSB[i].first ^ carry[i].first;
		msb[i].second = lambdaMSB[i].second ^ carry[i].second;
	}
}


void Meteor_funcPrivateCompare(const MEVectorSmallType &share_m, const vector<myType> &r, const MEVectorSmallType &beta, vector<smallType> &betaPrime, size_t size)
{
	log_print("Meteor private compare functionality");
	assert(share_m.size() == size * BIT_SIZE && "input error share_m");
	assert(r.size() == size && "Input error r");
	assert(beta.size() == size && "Input error beta");

	size_t sizeLong = size * BIT_SIZE;
	size_t index3, index2;
	MEVectorSmallType c(sizeLong), diff(sizeLong), twoBetaMinusOne(sizeLong), xMinusR(sizeLong);
	MESmallType a, tempM, tempN;
	smallType bit_r;

	for(int index2 = 0; index2 < size; index2++)
	{
		twoBetaMinusOne[index2*BIT_SIZE].first =  subtractModPrime[additionModPrime[beta[index2].first][beta[index2].first]][1];
		// (beta[index2].first + beta[index2].first - 1) % PRIME_NUMBER;
		twoBetaMinusOne[index2*BIT_SIZE].second.first = additionModPrime[beta[index2].second.first][beta[index2].second.first];
		//(beta[index2].second.first + beta[index2].second.first) % PRIME_NUMBER;
		twoBetaMinusOne[index2*BIT_SIZE].second.second = additionModPrime[beta[index2].second.second][beta[index2].second.second];
		//(beta[index2].second.second + beta[index2].second.second) % PRIME_NUMBER;

		for(size_t k = 0; k < BIT_SIZE; k++)
		{
			index3 = index2 * BIT_SIZE + k;
			twoBetaMinusOne[index3] = twoBetaMinusOne[index2*BIT_SIZE];

			bit_r = (smallType)( (r[index2] >> (BIT_SIZE-1-k)) & 1 );
			diff[index3] = share_m[index3];

			if(bit_r == 1){
				diff[index3].first = subtractModPrime[diff[index3].first][1];
			}
		}
	}

		// (-1)^beta * (x[i]-r[i]) in vectorization.
	Meteor_funcDotProduct(diff, twoBetaMinusOne, xMinusR, sizeLong);


	for(int index2 = 0; index2 < size; index2++)
	{
		a = make_pair(0, make_pair(0,0));
		for(size_t k = 0; k < BIT_SIZE; k++)
		{
			index3 = index2 * BIT_SIZE + k;
			c[index3] = a;
			tempM = share_m[index3];

			bit_r = (smallType)((r[index2] >> (BIT_SIZE-1-k)) & 1);

			if (bit_r == 0)
			{
				tempN = tempM;
			}
			else
			{
				tempN.first = subtractModPrime[1][tempM.first];
				tempN.second.first = subtractModPrime[0][tempM.second.first];
				tempN.second.second = subtractModPrime[0][tempM.second.second];
			}

			a.first = additionModPrime[a.first][tempN.first];
			a.second = addModPrime(a.second, tempN.second);

			c[index3].first = additionModPrime[additionModPrime[c[index3].first][xMinusR[index3].first]][1];
			c[index3].second = addModPrime(c[index3].second, xMinusR[index3].second);
		}

	}
	Meteor_funcCrunchMultiply(c, betaPrime, size);

}


void Meteor_funcRELUPrime(const MEVectorType &a, MEVectorSmallType &b, size_t size)
{
	log_print("Meteor funcRELUPrime");
	assert(a.size() == size && "input a is size error!");
	assert(b.size() == size && "output b is size error!");

	size_t bitSize = size * BIT_SIZE;
	RSSVectorMyType lambda(size), maskedInput(size);
	MEVectorSmallType lambdaBits(bitSize), lambdaMSB(size), sBits(bitSize), betaBit(size), betaPrimeMask(size), gamma(size);
	MEVectorSmallType carry(size);
	vector<myType> mPrime(size, 0), r(size, 0), carryThreshold(size, 0);
	vector<smallType> maskedCarry(size, 0);

	PrecomputeObject.getPreMSBObjects(lambda, lambdaBits, lambdaMSB, sBits, betaBit, betaPrimeMask, gamma, size);

	for(size_t i = 0; i < size; i++)
	{
		maskedInput[i] = a[i].second - lambda[i];
		if (partyNum == PARTY_A)
			maskedInput[i].first += a[i].first;
		if (partyNum == PARTY_C)
			maskedInput[i].second += a[i].first;
	}
	funcReconstruct(maskedInput, mPrime, size, "relu remasked input", false);

	for(size_t i = 0; i < size; i++)
	{
		carryThreshold[i] = myType(0) - 2*mPrime[i];
		r[i] = (carryThreshold[i] == 0) ? 0 : carryThreshold[i] - 1;
	}

	Meteor_funcPrivateCompare(sBits, r, betaPrimeMask, maskedCarry, size);

	for(size_t i = 0; i < size; i++)
	{
		if (carryThreshold[i] == 0)
		{
			carry[i].first = 0;
			carry[i].second = make_pair(smallType(0), smallType(0));
		}
		else
		{
			carry[i].first = betaBit[i].first ^ maskedCarry[i] ^ 1;
			carry[i].second = betaBit[i].second;
		}
		b[i].first = getMSB(mPrime[i]) ^ lambdaMSB[i].first ^ carry[i].first ^ 1;
		b[i].second = lambdaMSB[i].second ^ carry[i].second;
	}
}


void Meteor_funcRELU(const MEVectorType &a, MEVectorSmallType &temp, MEVectorType &b, size_t size)
{
	log_print("Meteor funcRELU");
	assert(a.size() == size && "input a is size error");
	assert(temp.size() == size && "input temp is size error");
	assert(b.size() == size && "output b is size error");

	Meteor_funcRELUPrime(a, temp, size);

	MEVectorType reluPrimeArithmetic(size);
	Meteor_funcBit2A(temp, reluPrimeArithmetic, size);
	Meteor_funcDotProduct(a, reluPrimeArithmetic, b, size, false, 0);
}


//Chunk wise maximum of a vector of size rows*columns and maximum is caclulated of every
//column number of elements. max is a vector of size rows, maxPrime, of rows*columns*columns;
void Meteor_funcMaxpool(const MEVectorType &a, MEVectorType &max,
						 size_t rows, size_t columns)
{
	log_print("funcMaxpool");
	assert(columns < 256 && "Pooling size has to be smaller than 8-bits");

	size_t size = rows*columns;
	MEVectorType diff(rows);
	MEVectorSmallType rp(rows); // dmpIndexShares(columns*size), temp(size);

	for (size_t i = 0; i < rows; ++i)
	{
		max[i] = a[i*columns];
		//max[i].second = a[i*columns].second;
	}
	for (size_t i = 1; i < columns; ++i)
	{
		for (size_t	j = 0; j < rows; ++j)
		{
			diff[j].first = max[j].first - a[j*columns + i].first;
			diff[j].second = max[j].second - a[j*columns + i].second;
		}

		Meteor_funcRELU(diff, rp, max, rows);
		//funcSelectBitShares(maxPrime, dmpIndexShares, rp, temp, rows, columns, i);

		//for (size_t i = 0; i < size; ++i)
		//	maxPrime[i] = temp[i];

		for (size_t	j = 0; j < rows; ++j)
		{
			max[j].first = max[j].first + a[j*columns + i].first;
			max[j].second = max[j].second + a[j*columns + i].second;
		}
	}
}




/******************************** Test ********************************/


void testMeteorPC(size_t size, size_t iter)
{
	for(int runs = 0; runs < iter; runs++){
		size_t bitSize = size * BIT_SIZE;
		RSSVectorMyType lambda(size);
		MEVectorSmallType lambdaBits(bitSize), lambdaMSB(size), sBits(bitSize), betaBit(size), betaPrimeMask(size), gamma(size);
		vector<smallType> betaPrime(size, 0);
		vector<myType> lambdaPlain(size), r(size);
		RSSVectorSmallType betaBitMask(size);
		vector<smallType> betaBitPlain(size);

		PrecomputeObject.getPreMSBObjects(lambda, lambdaBits, lambdaMSB, sBits, betaBit, betaPrimeMask, gamma, size);
		funcReconstruct(lambda, lambdaPlain, size, "pc lambda", false);
		for (size_t i = 0; i < size; ++i)
			betaBitMask[i] = betaBit[i].second;
		funcReconstructBit(betaBitMask, betaBitPlain, size, "pc beta bit mask", false);

		for (size_t i = 0; i < size; ++i)
			r[i] = (myType(1) << ((i + runs) % (BIT_SIZE - 1))) | myType(1);

		Meteor_funcPrivateCompare(sBits, r, betaPrimeMask, betaPrime, size);

		for (size_t i = 0; i < size; ++i)
		{
			myType s = lambdaPlain[i] << 1;
			smallType beta = betaBit[i].first ^ betaBitPlain[i];
			smallType expected = smallType(s >= r[i]) ^ beta ^ 1;
			assert(betaPrime[i] == expected && "Meteor PrivateCompare mismatch");
		}
	}
}


void testMeteorfuncCrunchMultiply(size_t size, size_t iter)
{
	for(int runs = 0; runs < iter; runs++){
		size_t sizeLong = size * BIT_SIZE;
		MEVectorSmallType a(sizeLong);
		RSSVectorSmallType aMask(sizeLong);
		vector<smallType> aMaskPlain(sizeLong), betaPrime(size, 0), expected(size, 1);

		PrecomputeObject.getRandomShares(aMask, sizeLong);
		funcReconstruct(aMask, aMaskPlain, sizeLong, "crunch multiply mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			for (size_t bit = 0; bit < BIT_SIZE; ++bit)
			{
				size_t index = i*BIT_SIZE + bit;
				smallType actual = ((i + runs) % 2 == 1 && bit == 5) ? 0 : smallType(1 + ((i + bit) % 2));
				if (actual == 0)
					expected[i] = 0;
				a[index] = make_pair(subtractModPrime[actual][aMaskPlain[index]], aMask[index]);
			}
		}

		Meteor_funcCrunchMultiply(a, betaPrime, size);
		for (size_t i = 0; i < size; ++i)
			assert(betaPrime[i] == expected[i] && "Meteor CrunchMultiply non-zero mismatch");
	}
}

void testThunderNMult(size_t size, size_t N,  size_t iter)
{
	assert(N == 4 && "ThunderNMult unit test currently exercises the 4-input online gate");
	size_t totalSize = size*N;
	for(int runs = 0; runs < iter; runs++){
		MEVectorType a(totalSize), c(size);
		RSSVectorMyType aMask(totalSize), cMask(size);
		vector<myType> aPlain(totalSize), cPlain(size);
		PrecomputeObject.getRandomShares(aMask, totalSize);
		for (size_t i = 0; i < totalSize; ++i)
			a[i] = make_pair(myType(i + 1), aMask[i]);

		Thunder_funcNMultGateOnline(a, c, totalSize);
		for (size_t i = 0; i < size; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(aMask, aPlain, totalSize, "thunder nmult a mask", false);
		funcReconstruct(cMask, cPlain, size, "thunder nmult c mask", false);
		for (size_t i = 0; i < size; ++i)
		{
			myType expected = 1;
			for (size_t j = 0; j < N; ++j)
			{
				size_t index = N*i + j;
				myType actual = a[index].first + aPlain[index];
				expected *= actual;
			}
			myType actualC = c[i].first + cPlain[i];
			assert(actualC == expected && "Thunder NMult mismatch");
		}
	}
}



void testMeteorNeighborMultly(size_t size, size_t iter)
{

	assert((size % 4) == 0 && "neighbor multiply test size must be a multiple of 4");
	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorSmallType a(size), c(size/4);
		RSSVectorSmallType aMask(size), cMask(size/4);
		vector<smallType> aPlain(size), cPlain(size/4);
		PrecomputeObject.getRandomShares(aMask, size);
		for (int i = 0; i < size; ++i)
			a[i] = make_pair(smallType((i + 1) % PRIME_NUMBER), aMask[i]);

		Meteor_funcMultiplyNeighbors(a, c, size);
		for (int i = 0; i < size/4; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(aMask, aPlain, size, "neighbor a mask", false);
		funcReconstruct(cMask, cPlain, size/4, "neighbor c mask", false);
		for (int i = 0; i < size/4; ++i)
		{
			smallType expected = 1;
			for (int j = 0; j < 4; ++j)
			{
				int index = 4*i + j;
				smallType actual = additionModPrime[a[index].first][aPlain[index]];
				expected = multiplicationModPrime[expected][actual];
			}
			smallType actualC = additionModPrime[c[i].first][cPlain[i]];
			assert(actualC == expected && "Meteor neighbor multiply mismatch");
		}
	}
}

void testMeteorBitProduct(size_t size, size_t iter)
{

	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorSmallType a(size), b(size), c(size);
		RSSVectorSmallType aMask(size), bMask(size), cMask(size);
		vector<smallType> aPlain(size), bPlain(size), cPlain(size);
		PrecomputeObject.getRandomBitShares(aMask, size);
		PrecomputeObject.getRandomBitShares(bMask, size);
		for (int i = 0; i < size; ++i)
		{
			a[i] = make_pair(smallType(i & 1), aMask[i]);
			b[i] = make_pair(smallType((i + 1) & 1), bMask[i]);
		}

		Meteor_funcDotProductBits(a, b, c, size);
		for (int i = 0; i < size; ++i)
			cMask[i] = c[i].second;
		funcReconstructBit(aMask, aPlain, size, "a mask", false);
		funcReconstructBit(bMask, bPlain, size, "b mask", false);
		funcReconstructBit(cMask, cPlain, size, "c mask", false);
		for (int i = 0; i < size; ++i)
		{
			smallType actualA = a[i].first ^ aPlain[i];
			smallType actualB = b[i].first ^ bPlain[i];
			smallType actualC = c[i].first ^ cPlain[i];
			assert(actualC == (actualA & actualB) && "Meteor bit product mismatch");
		}
	}
}

void testMeteorBit2A(size_t size, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorSmallType bits(size);
		MEVectorType arith(size);
		RSSVectorSmallType bitMask(size);
		RSSVectorMyType arithMask(size);
		vector<smallType> bitMaskPlain(size);
		vector<myType> arithMaskPlain(size);

		PrecomputeObject.getRandomBitShares(bitMask, size);
		for (size_t i = 0; i < size; ++i)
			bits[i] = make_pair(smallType(i & 1), bitMask[i]);

		Meteor_funcBit2A(bits, arith, size);
		for (size_t i = 0; i < size; ++i)
			arithMask[i] = arith[i].second;

		funcReconstructBit(bitMask, bitMaskPlain, size, "bit2a bit mask", false);
		funcReconstruct(arithMask, arithMaskPlain, size, "bit2a arithmetic mask", false);
		for (size_t i = 0; i < size; ++i)
		{
			smallType expected = bits[i].first ^ bitMaskPlain[i];
			myType actual = arith[i].first + arithMaskPlain[i];
			assert(actual == myType(expected) && "Meteor Bit2A mismatch");
		}
	}
}

void testPreMSBObjects(size_t size, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		size_t bitSize = size * BIT_SIZE;
		RSSVectorMyType lambda(size);
		MEVectorSmallType lambdaBits(bitSize), lambdaMSB(size), sBits(bitSize), betaBit(size), betaPrime(size), gamma(size);
		PrecomputeObject.getPreMSBObjects(lambda, lambdaBits, lambdaMSB, sBits, betaBit, betaPrime, gamma, size);

		vector<myType> lambdaPlain(size);
		RSSVectorSmallType lambdaBitsMask(bitSize), lambdaMSBMask(size), sBitsMask(bitSize), betaBitMask(size), betaPrimeMask(size), gammaMask(size);
		vector<smallType> lambdaBitsPlain(bitSize), lambdaMSBPlain(size), sBitsPlain(bitSize), betaBitPlain(size), betaPrimePlain(size), gammaPlain(size);

		for (size_t i = 0; i < size; ++i)
		{
			lambdaMSBMask[i] = lambdaMSB[i].second;
			betaBitMask[i] = betaBit[i].second;
			betaPrimeMask[i] = betaPrime[i].second;
			gammaMask[i] = gamma[i].second;
		}
		for (size_t i = 0; i < bitSize; ++i)
		{
			lambdaBitsMask[i] = lambdaBits[i].second;
			sBitsMask[i] = sBits[i].second;
		}

		funcReconstruct(lambda, lambdaPlain, size, "premsb lambda", false);
		funcReconstructBit(lambdaBitsMask, lambdaBitsPlain, bitSize, "premsb lambda bits mask", false);
		funcReconstructBit(lambdaMSBMask, lambdaMSBPlain, size, "premsb lambda msb mask", false);
		funcReconstruct(sBitsMask, sBitsPlain, bitSize, "premsb s bits mask", false);
		funcReconstructBit(betaBitMask, betaBitPlain, size, "premsb beta bit mask", false);
		funcReconstruct(betaPrimeMask, betaPrimePlain, size, "premsb beta prime mask", false);
		funcReconstruct(gammaMask, gammaPlain, size, "premsb gamma mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			smallType actualMSB = lambdaMSB[i].first ^ lambdaMSBPlain[i];
			assert(actualMSB == getMSB(lambdaPlain[i]) && "PreMSB lambda MSB mismatch");

			myType shifted = lambdaPlain[i] << 1;
			for (size_t bit = 0; bit < BIT_SIZE; ++bit)
			{
				size_t index = i*BIT_SIZE + bit;
				smallType actualLambdaBit = lambdaBits[index].first ^ lambdaBitsPlain[index];
				smallType expectedLambdaBit = smallType((lambdaPlain[i] >> (BIT_SIZE - 1 - bit)) & 1);
				assert(actualLambdaBit == expectedLambdaBit && "PreMSB lambda bit mismatch");

				smallType actualSBit = additionModPrime[sBits[index].first][sBitsPlain[index]];
				smallType expectedSBit = smallType((shifted >> (BIT_SIZE - 1 - bit)) & 1);
				assert(actualSBit == expectedSBit && "PreMSB shifted lambda bit mismatch");
			}

			smallType actualBetaBit = betaBit[i].first ^ betaBitPlain[i];
			smallType actualBetaPrime = additionModPrime[betaPrime[i].first][betaPrimePlain[i]];
			assert(actualBetaPrime == actualBetaBit && "PreMSB beta representation mismatch");

			smallType actualGamma = additionModPrime[gamma[i].first][gammaPlain[i]];
			assert(actualGamma != 0 && "PreMSB gamma must be non-zero");
		}
	}
}

void testMeteorSecMSB(size_t size, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		size_t bitSize = size * BIT_SIZE;
		RSSVectorMyType lambda(size);
		MEVectorSmallType lambdaBits(bitSize), lambdaMSB(size), sBits(bitSize), betaBit(size), betaPrime(size), gamma(size), msb(size);
		MEVectorType x(size);
		vector<myType> lambdaPlain(size);
		RSSVectorSmallType msbMask(size);
		vector<smallType> msbMaskPlain(size);

		PrecomputeObject.getPreMSBObjects(lambda, lambdaBits, lambdaMSB, sBits, betaBit, betaPrime, gamma, size);
		funcReconstruct(lambda, lambdaPlain, size, "secmsb lambda", false);

		for (size_t i = 0; i < size; ++i)
		{
			myType actual;
			switch (i % 6)
			{
				case 0: actual = 0; break;
				case 1: actual = 1; break;
				case 2: actual = LARGEST_NEG - 1; break;
				case 3: actual = LARGEST_NEG; break;
				case 4: actual = MINUS_ONE; break;
				default: actual = myType(0x123456789abcdef0ULL + i + runs); break;
			}
			x[i] = make_pair(actual - lambdaPlain[i], lambda[i]);
		}

		Meteor_funcSecMSB(x, lambdaBits, lambdaMSB, sBits, betaBit, betaPrime, gamma, msb, size);

		for (size_t i = 0; i < size; ++i)
			msbMask[i] = msb[i].second;
		funcReconstructBit(msbMask, msbMaskPlain, size, "secmsb msb mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			myType actual = x[i].first + lambdaPlain[i];
			smallType actualMSB = msb[i].first ^ msbMaskPlain[i];
			assert(actualMSB == getMSB(actual) && "Meteor SecMSB mismatch");
		}
	}
}

void testMeteorSmallDotProduct(size_t size, size_t iter)
{

	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorSmallType a(size), b(size), c(size);
		RSSVectorSmallType aMask(size), bMask(size), cMask(size);
		vector<smallType> aPlain(size), bPlain(size), cPlain(size);
		PrecomputeObject.getRandomShares(aMask, size);
		PrecomputeObject.getRandomShares(bMask, size);
		for (int i = 0; i < size; ++i)
		{
			a[i] = make_pair(smallType((i + 1) % PRIME_NUMBER), aMask[i]);
			b[i] = make_pair(smallType((i + 2) % PRIME_NUMBER), bMask[i]);
		}

		Meteor_funcDotProduct(a, b, c, size);
		for (int i = 0; i < size; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(aMask, aPlain, size, "a mask", false);
		funcReconstruct(bMask, bPlain, size, "b mask", false);
		funcReconstruct(cMask, cPlain, size, "c mask", false);
		for (int i = 0; i < size; ++i)
		{
			smallType actualA = additionModPrime[a[i].first][aPlain[i]];
			smallType actualB = additionModPrime[b[i].first][bPlain[i]];
			smallType actualC = additionModPrime[c[i].first][cPlain[i]];
			assert(actualC == multiplicationModPrime[actualA][actualB] && "Meteor field product mismatch");
		}
	}
}


void testMeteorDotProduct(size_t size, size_t iter)
{

	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorType a(size), b(size), c(size);
		RSSVectorMyType aMask(size), bMask(size), cMask(size);
		vector<myType> aPlain(size), bPlain(size), cPlain(size);
		PrecomputeObject.getRandomShares(aMask, size);
		PrecomputeObject.getRandomShares(bMask, size);
		for (int i = 0; i < size; ++i)
		{
			a[i] = make_pair(myType(i + 1), aMask[i]);
			b[i] = make_pair(myType(i + 2), bMask[i]);
		}

		Meteor_funcDotProduct(a, b, c, size, false, FLOAT_PRECISION);
		for (int i = 0; i < size; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(aMask, aPlain, size, "a mask", false);
		funcReconstruct(bMask, bPlain, size, "b mask", false);
		funcReconstruct(cMask, cPlain, size, "c mask", false);
		for (int i = 0; i < size; ++i)
		{
			myType actualA = a[i].first + aPlain[i];
			myType actualB = b[i].first + bPlain[i];
			myType actualC = c[i].first + cPlain[i];
			assert(actualC == actualA * actualB && "Meteor integer product mismatch");
		}
	}
}

void testMeteorMatMul(size_t rows, size_t common_dim, size_t columns, size_t iter)
{

	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorType a(rows*common_dim), b(common_dim*columns), c(rows*columns);
		RSSVectorMyType aMask(rows*common_dim), bMask(common_dim*columns), cMask(rows*columns);
		vector<myType> aPlain(rows*common_dim), bPlain(common_dim*columns), cPlain(rows*columns);
		PrecomputeObject.getRandomShares(aMask, rows*common_dim);
		PrecomputeObject.getRandomShares(bMask, common_dim*columns);
		for (size_t i = 0; i < rows*common_dim; ++i)
			a[i] = make_pair(myType(i + 1), aMask[i]);
		for (size_t i = 0; i < common_dim*columns; ++i)
			b[i] = make_pair(myType(i + 2), bMask[i]);

		Meteor_funcMatMul(a, b, c, rows, common_dim, columns, 0, 0, 0);
		for (size_t i = 0; i < rows*columns; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(aMask, aPlain, rows*common_dim, "matmul a mask", false);
		funcReconstruct(bMask, bPlain, common_dim*columns, "matmul b mask", false);
		funcReconstruct(cMask, cPlain, rows*columns, "matmul c mask", false);

		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < columns; ++j)
			{
				myType expected = 0;
				for (size_t k = 0; k < common_dim; ++k)
				{
					myType actualA = a[i*common_dim + k].first + aPlain[i*common_dim + k];
					myType actualB = b[k*columns + j].first + bPlain[k*columns + j];
					expected += actualA * actualB;
				}
				myType actualC = c[i*columns + j].first + cPlain[i*columns + j];
				assert(actualC == expected && "Meteor MatMul mismatch");
			}
		}
	}
}

void testMeteorTruncation(size_t size, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorType a(size), b(size);
		RSSVectorMyType aMask(size), bMask(size);
		vector<myType> aMaskPlain(size), bMaskPlain(size), actual(size);
		vector<int64_t> signedCases = {
			0, 1, -1,
			(1 << FLOAT_PRECISION) - 1,
			(1 << FLOAT_PRECISION),
			(1 << FLOAT_PRECISION) + 1,
			-((1 << FLOAT_PRECISION) - 1),
			-(1 << FLOAT_PRECISION),
			-((1 << FLOAT_PRECISION) + 1),
			1234567, -1234567,
			(1LL << 30), -(1LL << 30)
		};
		PrecomputeObject.getRandomShares(aMask, size);
		funcReconstruct(aMask, aMaskPlain, size, "truncation input mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			int64_t signedValue = signedCases[(i + runs) % signedCases.size()];
			actual[i] = static_cast<myType>(signedValue);
			a[i] = make_pair(actual[i] - aMaskPlain[i], aMask[i]);
		}

		Meteor_funcTruncate(a, b, size, FLOAT_PRECISION);
		for (size_t i = 0; i < size; ++i)
			bMask[i] = b[i].second;
		funcReconstruct(bMask, bMaskPlain, size, "truncation output mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			myType opened = b[i].first + bMaskPlain[i];
			myType expected = signedRightShift(actual[i], FLOAT_PRECISION);
			assert((opened == expected || opened == expected + 1) && "Meteor faithful truncation mismatch");
		}
	}
}

void testMeteorTruncatingDotProduct(size_t size, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorType a(size), b(size), c(size);
		RSSVectorMyType aMask(size), bMask(size), cMask(size);
		vector<myType> aMaskPlain(size), bMaskPlain(size), cMaskPlain(size), actualA(size), actualB(size);
		PrecomputeObject.getRandomShares(aMask, size);
		PrecomputeObject.getRandomShares(bMask, size);
		funcReconstruct(aMask, aMaskPlain, size, "trunc dot a mask", false);
		funcReconstruct(bMask, bMaskPlain, size, "trunc dot b mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			int64_t signedA = int64_t((i*17 + runs*3) % 401) - 200;
			int64_t signedB = int64_t((i*29 + runs*5) % 307) - 153;
			actualA[i] = static_cast<myType>(signedA);
			actualB[i] = static_cast<myType>(signedB);
			a[i] = make_pair(actualA[i] - aMaskPlain[i], aMask[i]);
			b[i] = make_pair(actualB[i] - bMaskPlain[i], bMask[i]);
		}

		Meteor_funcDotProduct(a, b, c, size, true, FLOAT_PRECISION);
		for (size_t i = 0; i < size; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(cMask, cMaskPlain, size, "trunc dot c mask", false);

		for (size_t i = 0; i < size; ++i)
		{
			myType product = actualA[i] * actualB[i];
			myType expected = signedRightShift(product, FLOAT_PRECISION);
			myType opened = c[i].first + cMaskPlain[i];
			assert((opened == expected || opened == expected + 1) && "Meteor truncating dot product mismatch");
		}
	}
}

void testMeteorTruncatingMatMul(size_t rows, size_t common_dim, size_t columns, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorType a(rows*common_dim), b(common_dim*columns), c(rows*columns);
		RSSVectorMyType aMask(rows*common_dim), bMask(common_dim*columns), cMask(rows*columns);
		vector<myType> aMaskPlain(rows*common_dim), bMaskPlain(common_dim*columns), cMaskPlain(rows*columns);
		vector<myType> actualA(rows*common_dim), actualB(common_dim*columns);
		PrecomputeObject.getRandomShares(aMask, rows*common_dim);
		PrecomputeObject.getRandomShares(bMask, common_dim*columns);
		funcReconstruct(aMask, aMaskPlain, rows*common_dim, "trunc matmul a mask", false);
		funcReconstruct(bMask, bMaskPlain, common_dim*columns, "trunc matmul b mask", false);

		for (size_t i = 0; i < rows*common_dim; ++i)
		{
			int64_t signedValue = int64_t((i*11 + runs*7) % 257) - 128;
			actualA[i] = static_cast<myType>(signedValue);
			a[i] = make_pair(actualA[i] - aMaskPlain[i], aMask[i]);
		}
		for (size_t i = 0; i < common_dim*columns; ++i)
		{
			int64_t signedValue = int64_t((i*19 + runs*13) % 193) - 96;
			actualB[i] = static_cast<myType>(signedValue);
			b[i] = make_pair(actualB[i] - bMaskPlain[i], bMask[i]);
		}

		Meteor_funcMatMul(a, b, c, rows, common_dim, columns, 0, 0, FLOAT_PRECISION);
		for (size_t i = 0; i < rows*columns; ++i)
			cMask[i] = c[i].second;
		funcReconstruct(cMask, cMaskPlain, rows*columns, "trunc matmul c mask", false);

		for (size_t i = 0; i < rows; ++i)
		{
			for (size_t j = 0; j < columns; ++j)
			{
				myType sum = 0;
				for (size_t k = 0; k < common_dim; ++k)
					sum += actualA[i*common_dim + k] * actualB[k*columns + j];

				myType expected = signedRightShift(sum, FLOAT_PRECISION);
				myType opened = c[i*columns + j].first + cMaskPlain[i*columns + j];
				assert((opened == expected || opened == expected + 1) && "Meteor truncating MatMul mismatch");
			}
		}
	}
}


void testMeteorRelu(size_t size, size_t iter)
{
	for (int runs = 0; runs < iter; ++runs)
	{
		MEVectorType a(size), b(size);
		MEVectorSmallType reluPrime(size);
		RSSVectorMyType aMask(size), bMask(size);
		vector<myType> aMaskPlain(size), bMaskPlain(size), actual(size);

		PrecomputeObject.getRandomShares(aMask, size);
		funcReconstruct(aMask, aMaskPlain, size, "relu input mask", false);
		for (size_t i = 0; i < size; ++i)
		{
			switch (i % 6)
			{
				case 0: actual[i] = 0; break;
				case 1: actual[i] = 1; break;
				case 2: actual[i] = LARGEST_NEG - 1; break;
				case 3: actual[i] = LARGEST_NEG; break;
				case 4: actual[i] = MINUS_ONE; break;
				default: actual[i] = myType(12345 + runs + i); break;
			}
			a[i] = make_pair(actual[i] - aMaskPlain[i], aMask[i]);
		}

		Meteor_funcRELU(a, reluPrime, b, size);
		for (size_t i = 0; i < size; ++i)
			bMask[i] = b[i].second;
		funcReconstruct(bMask, bMaskPlain, size, "relu output mask", false);
		for (size_t i = 0; i < size; ++i)
		{
			myType expected = getMSB(actual[i]) ? 0 : actual[i];
			myType opened = b[i].first + bMaskPlain[i];
			assert(opened == expected && "Meteor ReLU mismatch");
		}
	}
}



void testMeteorRELUPrime(size_t size, size_t iter)
{
	for(int runs = 0; runs < iter; runs++)
	{
		MEVectorType a(size);
		MEVectorSmallType b(size);
		RSSVectorMyType aMask(size);
		RSSVectorSmallType bMask(size);
		vector<myType> aMaskPlain(size), actual(size);
		vector<smallType> bMaskPlain(size);

		PrecomputeObject.getRandomShares(aMask, size);
		funcReconstruct(aMask, aMaskPlain, size, "relu prime input mask", false);
		for (size_t i = 0; i < size; ++i)
		{
			switch (i % 6)
			{
				case 0: actual[i] = 0; break;
				case 1: actual[i] = 1; break;
				case 2: actual[i] = LARGEST_NEG - 1; break;
				case 3: actual[i] = LARGEST_NEG; break;
				case 4: actual[i] = MINUS_ONE; break;
				default: actual[i] = myType(67890 + runs + i); break;
			}
			a[i] = make_pair(actual[i] - aMaskPlain[i], aMask[i]);
		}

		Meteor_funcRELUPrime(a, b, size);
		for (size_t i = 0; i < size; ++i)
			bMask[i] = b[i].second;
		funcReconstructBit(bMask, bMaskPlain, size, "relu prime output mask", false);
		for (size_t i = 0; i < size; ++i)
		{
			smallType opened = b[i].first ^ bMaskPlain[i];
			smallType expected = getMSB(actual[i]) ^ 1;
			assert(opened == expected && "Meteor ReLUPrime mismatch");
		}
	}
}


void testMeteorMaxpool(size_t ih, size_t iw, size_t Din, size_t f, size_t S, size_t B, size_t iter)
{
	size_t ow 		= (((iw-f)/S)+1);
	size_t oh		= (((ih-f)/S)+1);

	size_t inputSize = iw*ih*Din*B;
	size_t outputSize = ow*oh*Din*B;
	size_t windowSize = f*f;
	MEVectorType a(inputSize), b(outputSize), temp1(outputSize*windowSize);
	RSSVectorMyType inputMask(inputSize), outputMask(outputSize);
	vector<myType> inputMaskPlain(inputSize), outputMaskPlain(outputSize);
	vector<myType> actual(inputSize);
	vector<int64_t> actualSigned(inputSize);
	size_t sizeBeta = iw;
	size_t sizeD 	= sizeBeta*ih;
	size_t sizeB 	= sizeD*Din;
	size_t counter 	= 0;

	for (int runs = 0; runs < iter; ++runs)
	{
		PrecomputeObject.getRandomShares(inputMask, inputSize);
		funcReconstruct(inputMask, inputMaskPlain, inputSize, "maxpool input mask", false);

		for (size_t i = 0; i < inputSize; ++i)
		{
			actualSigned[i] = int64_t((i*13 + runs*7) % 31) - 15;
			actual[i] = static_cast<myType>(actualSigned[i]);
			a[i] = make_pair(actual[i] - inputMaskPlain[i], inputMask[i]);
		}

		counter = 0;
		for (int b = 0; b < B; ++b)
			for (size_t r = 0; r < Din; ++r)
				for (size_t beta = 0; beta < ih-f+1; beta+=S)
					for (size_t alpha = 0; alpha < iw-f+1; alpha+=S)
						for (int q = 0; q < f; ++q)
							for (int p = 0; p < f; ++p)
							{
								temp1[counter].first =
								a[b*sizeB + r*sizeD + (beta + q)*sizeBeta + (alpha + p)].first;
								temp1[counter].second.first =
								a[b*sizeB + r*sizeD + (beta + q)*sizeBeta + (alpha + p)].second.first;
								temp1[counter].second.second =
								a[b*sizeB + r*sizeD + (beta + q)*sizeBeta + (alpha + p)].second.second;
								counter++;
							}

		Meteor_funcMaxpool(temp1, b, outputSize, windowSize);
		for (size_t i = 0; i < outputSize; ++i)
			outputMask[i] = b[i].second;
		funcReconstruct(outputMask, outputMaskPlain, outputSize, "maxpool output mask", false);

		counter = 0;
		for (int batch = 0; batch < B; ++batch)
			for (size_t r = 0; r < Din; ++r)
				for (size_t beta = 0; beta < ih-f+1; beta+=S)
					for (size_t alpha = 0; alpha < iw-f+1; alpha+=S)
					{
						int64_t expectedSigned = -1000;
						for (int q = 0; q < f; ++q)
							for (int p = 0; p < f; ++p)
							{
								size_t inputIndex = batch*sizeB + r*sizeD + (beta + q)*sizeBeta + (alpha + p);
								int64_t candidate = actualSigned[inputIndex];
								if (candidate > expectedSigned)
									expectedSigned = candidate;
							}

						myType opened = b[counter].first + outputMaskPlain[counter];
						assert(opened == static_cast<myType>(expectedSigned) && "Meteor Maxpool mismatch");
						counter++;
					}
	}
}

void testMeteorConvolution(size_t iw, size_t ih, size_t Din, size_t Dout, size_t f, size_t S, size_t P, size_t B, size_t iter)
{
	size_t ow 		= (((iw-f+2*P)/S)+1);
	size_t oh		= (((ih-f+2*P)/S)+1);
	size_t tempSize = ow*oh;

	MEVectorType a(iw*ih*Din*B, make_pair(0, make_pair(0,0)));
	MEVectorType b(f*f*Din*Dout, make_pair(0, make_pair(0,0)));
	MEVectorType ans(ow*oh*Dout*B, make_pair(0, make_pair(0,0)));
	MEVectorType c(Dout, make_pair(0, make_pair(0,0)));

	for (int runs = 0; runs < iter; ++runs)
	{
	// 	//Reshape activations
	 	MEVectorType temp1((iw+2*P)*(ih+2*P)*Din*B, make_pair(0, make_pair(0,0)));
	 	zeroPad(a, temp1, iw, ih, P, Din, B);

	 	//Reshape for convolution
	 	MEVectorType temp2((f*f*Din) * (ow * oh * B));
	 	//convToMult(temp1, temp2, (iw+2*P), (ih+2*P), f, Din, S, B);
		{
		size_t loc_input, loc_output;
		for (size_t i = 0; i < B; ++i)
			for (size_t j = 0; j < oh; j++)
				for (size_t k = 0; k < ow; k++)
				{
					loc_output = (i*ow*oh + j*ow + k);
					for (size_t l = 0; l < Din; ++l)
					{
						loc_input = i*(iw+2*P)*(ih+2*P)*Din + l*(iw+2*P)*(ih+2*P) + j*S*(iw+2*P) + k*S;
						for (size_t a = 0; a < f; ++a)			//filter height
							for (size_t b = 0; b < f; ++b){		//filter width
								temp2[(l*f*f + a*f + b)*ow*oh*B + loc_output].first = temp1[loc_input + a*(iw+2*P) + b].first;
								temp2[(l*f*f + a*f + b)*ow*oh*B + loc_output].second = temp1[loc_input + a*(iw+2*P) + b].second;
							}
					}
				}
		}
		MEVectorType temp3(Dout * (ow*oh*B));
		Meteor_funcMatMul(b, temp2, temp3, Dout, (f*f*Din), (ow*oh*B), 0, 1, FLOAT_PRECISION);


	// 	//Add biases and meta-transpose
		for (size_t i = 0; i < B; ++i)
			for (size_t j = 0; j < Dout; ++j)
				for (size_t k = 0; k < tempSize; ++k){
					ans[i*Dout*tempSize + j*tempSize + k].first = temp3[j*B*tempSize + i*tempSize + k].first + c[j].first;
					ans[i*Dout*tempSize + j*tempSize + k].second = temp3[j*B*tempSize + i*tempSize + k].second + c[j].second;
				 }
	 }
}

void test_MeteorBatchNorm(size_t numBatches, size_t inputSize, size_t iter)
{
	size_t B = numBatches;
	size_t m = inputSize;

	MEVectorType gamma(numBatches, make_pair(0, make_pair(0,0))), g_repeat(B*m), input(B*m, make_pair(0, make_pair(0,0))), activations(B*m, make_pair(0, make_pair(0,0))), beta(B, make_pair(0, make_pair(0,0)));

	for(int runs = 0; runs < iter; runs++){
		for (int i = 0; i < B; ++i){
			for (int j = 0; j < m; ++j){
				g_repeat[i*m+j] = gamma[i];
			}
		}

		Meteor_funcDotProduct(g_repeat, input, activations, B*m, true, FLOAT_PRECISION);

		for (int i = 0; i < B; ++i){
			for (int j = 0; j < m; ++j){
				activations[i*m+j].first = activations[i*m+j].first + beta[i].first;
				activations[i*m+j].second = activations[i*m+j].second + beta[i].second;
			}
		}
	}

}
