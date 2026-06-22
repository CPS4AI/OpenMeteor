
#ifndef AESOBJECT_H
#define AESOBJECT_H

#pragma once
#include <algorithm>
#include <array>
#include <stdint.h>
#ifdef METEOR_AESNI_BACKEND
#include "TedKrovetzAesNiWrapperC.h"
#else
#include "../util/aes.h"
#endif
#include "globals.h"


class AESObject
{
private:
	//AES variables
#ifdef METEOR_AESNI_BACKEND
	__m128i pseudoRandomString[RANDOM_COMPUTE];
	__m128i tempSecComp[RANDOM_COMPUTE];
	unsigned long rCounter = -1;
	AES_KEY_TED aes_key;

	//Extraction variables
	__m128i random8BitNumber {0};
	uint8_t random8BitCounter = 0;
	__m128i random64BitNumber {0};
	uint8_t random64BitCounter = 0;

	//Private extraction functions
	__m128i newRandomNumber();
#else
	std::array<uint8_t, 16> random8BitNumber;
	std::array<uint8_t, 16> random64BitNumber;
	uint8_t random8BitCounter = 0;
	uint8_t random64BitCounter = 0;
	uint64_t rCounter = UINT64_MAX;
	AES_KEY aes_key;

	//Private extraction functions
	std::array<uint8_t, 16> newRandomNumber();
#endif

	//Private helper functions
	smallType AES_random(int i);

public:
	//Constructor
	AESObject(char* filename);
	
	//Randomness functions
	myType get64Bits();
	smallType get8Bits();

	//Other randomness functions
	smallType randModPrime();
	smallType randNonZeroModPrime();
};



#endif
