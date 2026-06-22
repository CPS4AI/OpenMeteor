#include <iostream>
#include <string>
#include "AESObject.h"
#include "Precompute.h"
#include "secondary.h"
#include "connect.h"
#include "NeuralNetConfig.h"
#include "NeuralNetwork.h"
#include "unitTests.h"


int partyNum;
AESObject* aes_indep;
AESObject* aes_next;
AESObject* aes_prev;
Precompute PrecomputeObject;


int main(int argc, char** argv)
{
/****************************** PREPROCESSING ******************************/ 
	parseInputs(argc, argv);
	NeuralNetConfig* config = new NeuralNetConfig(NUM_ITERATIONS);
	string network, dataset, runMode, unitTest;
	bool PRELOADING = true;

/****************************** SELECT NETWORK ******************************/
	//Run mode {unit, inference, preloaded, train}
	//Network {SecureML, Sarda, MiniONN, LeNet, AlexNet, and VGG16}
	//Dataset {MNIST, CIFAR10, and ImageNet}
	runMode = "unit";
	unitTest = "MeteorRELU";
	if (argc == 8)
	{network = argv[6]; dataset = argv[7];}
	else if (argc == 9)
	{runMode = argv[6]; network = argv[7]; dataset = argv[8];}
	else if (argc == 10)
	{runMode = argv[6]; network = argv[7]; dataset = argv[8]; unitTest = argv[9];}
	else
	{
		network = "SecureML";
		dataset = "MNIST";
	}
	if (runMode.compare("unit") != 0 and
		runMode.compare("inference") != 0 and
		runMode.compare("preloaded") != 0 and
		runMode.compare("train") != 0)
	{
		cout << "Run mode must be unit, inference, preloaded, or train" << endl;
		return -1;
	}
	selectNetwork(network, dataset, config);
	config->checkNetwork();
	NeuralNetwork* net = new NeuralNetwork(config);

/****************************** AES SETUP and SYNC ******************************/ 
	aes_indep = new AESObject(argv[3]);
	aes_next = new AESObject(argv[4]);
	aes_prev = new AESObject(argv[5]);

	initializeCommunication(argv[2], partyNum);
	synchronize(2000000);
	cout << "connected success" << endl;

/****************************** RUN NETWORK/UNIT TESTS ******************************/ 
	start_m();
	if (runMode.compare("unit") == 0)
	{
		runTest("Test", unitTest, network);
	}
	else if (runMode.compare("inference") == 0)
	{
		PRELOADING = false;
		network += " test";
		test(PRELOADING, network, net);
	}
	else if (runMode.compare("preloaded") == 0)
	{
		// Supported for preloaded MNIST networks generated with MINI_BATCH_SIZE == 128.
		PRELOADING = true;
		network += " preloaded";
		preload_network(PRELOADING, network, net);
		test(PRELOADING, network, net);
	}
	else if (runMode.compare("train") == 0)
	{
		network += " train";
		train(net);
	}

	end_m(network);
	cout << "----------------------------------------------" << endl;
	cout << "Run details: " << NUM_OF_PARTIES << "PC (P" << partyNum
		 << "), " << NUM_ITERATIONS << " iterations, batch size " << MINI_BATCH_SIZE << endl
		 << "Running semi-honest " << runMode << " mode for " << network << " on " << dataset << " dataset" << endl;
	cout << "----------------------------------------------" << endl << endl;  

	printNetwork(net);

/****************************** CLEAN-UP ******************************/ 
	delete aes_indep;
	delete aes_next;
	delete aes_prev;
	delete config;
	delete net;
	deleteObjects();

	return 0;
}




