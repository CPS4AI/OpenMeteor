
#pragma once
#ifdef METEOR_AESNI_BACKEND
#include "TedKrovetzAesNiWrapperC.h"
#endif
#include <cstring>
#include <iostream>
#include <fstream>
#include "AESObject.h"

using namespace std;


AESObject::AESObject(char* filename)
{
	ifstream f(filename);
	string str { istreambuf_iterator<char>(f), istreambuf_iterator<char>() };
	f.close();
	unsigned char common_aes_key[32];
	memset(common_aes_key, 0, sizeof(common_aes_key));
	memcpy(common_aes_key, str.data(), std::min<size_t>(str.size(), sizeof(common_aes_key)));
#ifdef METEOR_AESNI_BACKEND
	AES_set_encrypt_key((unsigned char*)common_aes_key, 256, &aes_key);
#else
	private_AES_set_encrypt_key(common_aes_key, 256, &aes_key);
	random8BitNumber.fill(0);
	random64BitNumber.fill(0);
#endif
}

#ifdef METEOR_AESNI_BACKEND
__m128i AESObject::newRandomNumber()
{
	rCounter++;
	if (rCounter % RANDOM_COMPUTE == 0)//generate more random seeds
	{
		for (int i = 0; i < RANDOM_COMPUTE; i++)
			tempSecComp[i] = _mm_set1_epi32(rCounter+i);//not exactly counter mode - (rcounter+i,rcouter+i,rcounter+i,rcounter+i)
		AES_ecb_encrypt_chunk_in_out(tempSecComp, pseudoRandomString, RANDOM_COMPUTE, &aes_key);
	}
	return pseudoRandomString[rCounter%RANDOM_COMPUTE];
}
#else
std::array<uint8_t, 16> AESObject::newRandomNumber()
{
	rCounter++;
	std::array<uint8_t, 16> input;
	std::array<uint8_t, 16> output;
	uint32_t counter32 = static_cast<uint32_t>(rCounter);
	for (int lane = 0; lane < 4; ++lane)
	{
		input[lane*4 + 0] = static_cast<uint8_t>(counter32);
		input[lane*4 + 1] = static_cast<uint8_t>(counter32 >> 8);
		input[lane*4 + 2] = static_cast<uint8_t>(counter32 >> 16);
		input[lane*4 + 3] = static_cast<uint8_t>(counter32 >> 24);
	}
	AES_encrypt(input.data(), output.data(), &aes_key);
	return output;
}
#endif


myType AESObject::get64Bits()
{
	myType ret;

	if (random64BitCounter == 0)
		random64BitNumber = newRandomNumber();
	
	int x = random64BitCounter % 2;
#ifdef METEOR_AESNI_BACKEND
	uint64_t *temp = (uint64_t*)&random64BitNumber;

	switch(x)
	{
    case 0 : ret = (myType)temp[1];
             break;
    case 1 : ret = (myType)temp[0];
             break;
	}
#else
	uint64_t low = 0, high = 0;
	memcpy(&low, random64BitNumber.data(), sizeof(uint64_t));
	memcpy(&high, random64BitNumber.data() + sizeof(uint64_t), sizeof(uint64_t));
	ret = (x == 0) ? static_cast<myType>(high) : static_cast<myType>(low);
#endif

	random64BitCounter++;	
	if (random64BitCounter == 2)
		random64BitCounter = 0;

	return ret;
}

smallType AESObject::get8Bits()
{
	smallType ret;

	if (random8BitCounter == 0)
		random8BitNumber = newRandomNumber();
	
#ifdef METEOR_AESNI_BACKEND
	uint8_t *temp = (uint8_t*)&random8BitNumber;
	ret = (smallType)temp[random8BitCounter];
#else
	ret = static_cast<smallType>(random8BitNumber[random8BitCounter]);
#endif

	random8BitCounter++;	
	if (random8BitCounter == 16)
		random8BitCounter = 0;

	return ret;
}

smallType AESObject::randModPrime()
{
	smallType ret;
	
	do
	{
		ret = get8Bits();
	} while (ret >= BOUNDARY);

	return (ret % PRIME_NUMBER);
}

smallType AESObject::randNonZeroModPrime()
{
	smallType ret;
	do
	{
		ret = randModPrime();
	} while (ret == 0);

	return ret;
}


smallType AESObject::AES_random(int i)
{
	smallType ret;
	do
	{
		ret = get8Bits();
	} while (ret >= ((256/i) * i));

	return (ret % i);
}
