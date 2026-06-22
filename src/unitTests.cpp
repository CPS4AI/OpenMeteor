
#include "Functionalities.h"


void runTest(string str, string whichTest, string &network)
{
	if (str.compare("Debug") == 0)
	{
		cout << "Debug mode is not available in unit test runner" << endl;
	}
	else if (str.compare("Test") == 0)
	{

		if(whichTest.compare("MeteorRELU") == 0){
			network = "Test Meteor RELU";
			testMeteorRelu(8, 3);
		}
		else if(whichTest.compare("MeteorRELUPrime") == 0){
			network = "Test Meteor RELUPrime";
			testMeteorRELUPrime(8, 3);
		}
		else if(whichTest.compare("MeteorPC") == 0){
			network = "Test Meteor PC";
			testMeteorPC(3, 3);
		}
		else if(whichTest.compare("BitProduct")==0){
			network = "Test BitProduct";
			testMeteorBitProduct(10, NUM_ITERATIONS);
		}
		else if(whichTest.compare("MeteorBit2A")==0){
			network = "Test Meteor Bit2A";
			testMeteorBit2A(10, 3);
		}
		else if(whichTest.compare("PreMSBObjects")==0){
			network = "Test PreMSB Objects";
			testPreMSBObjects(2, 3);
		}
		else if(whichTest.compare("MeteorSecMSB")==0){
			network = "Test Meteor SecMSB";
			testMeteorSecMSB(8, 3);
		}
		else if(whichTest.compare("MeteorCrunchMultiply")==0){
			network = "Test Meteor CrunchMultiply";
			testMeteorfuncCrunchMultiply(3, 3);
		}
		else if (whichTest.compare("MeteorDotProduct") == 0){
			network = "Test Meteor DotProduct";
			testMeteorDotProduct(10, 3);
		}
		else if (whichTest.compare("MeteorTruncation") == 0){
			network = "Test Meteor faithful truncation";
			testMeteorTruncation(13, 5);
		}
		else if (whichTest.compare("MeteorTruncatingDotProduct") == 0){
			network = "Test Meteor truncating DotProduct";
			testMeteorTruncatingDotProduct(10, 3);
		}
		else if (whichTest.compare("MeteorTruncatingMatMul") == 0){
			network = "Test Meteor truncating MatMul";
			testMeteorTruncatingMatMul(2, 3, 2, 3);
		}
		else if (whichTest.compare("MeteorSmallDotProduct") == 0){
			network = "Test Meteor Small DotProduct";
			testMeteorSmallDotProduct(10, 3);
		}
		else if (whichTest.compare("MeteorMatMul") == 0){
			network = "Test Meteor MatMul";
			testMeteorMatMul(2, 3, 2, 3);
		}
		else if (whichTest.compare("MeteorNeighborMultiply") == 0){
			network = "Test Meteor Neighbor Multiply";
			testMeteorNeighborMultly(8, 3);
		}
		else if (whichTest.compare("ThunderNMult") == 0){
			network = "Test Thunder N Mult";
			testThunderNMult(2, 4, 3);
		}
		else if(whichTest.compare("MeteorMaxpool") == 0){
			network = "Test Meteor Maxpool";
			testMeteorMaxpool(4, 4, 2, 2, 2, 2, 2);
		}
		else if (whichTest.compare("Conv") == 0)
		{
			network = "Test Meteor Conv";
			testMeteorConvolution(28, 28, 1, 20, 3, 1, 0, MINI_BATCH_SIZE, NUM_ITERATIONS);
		}
		else if (whichTest.compare("BN") == 0)
		{
			network = "Test Meteor BN";
			test_MeteorBatchNorm(MINI_BATCH_SIZE, 784, NUM_ITERATIONS);
		}
		else
			assert(false && "Unknown test mode selected");
	}
	else
		assert(false && "Only Debug or Test mode supported");
}
