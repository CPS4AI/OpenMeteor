/*
 * basicSockets.h
 *
 *  Created on: Aug 3, 2015
 *      Author: froike (Roi Inbar) 
 * 	Modified: Aner Ben-Efraim
 * 
 */
#include <stdio.h>
//#include <stropts.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#ifdef __linux__
#include <linux/netdevice.h>
#else
#include <net/if.h>
#endif
#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string>
#include <array>
#include "globals.h"
using namespace std;

#ifndef BMRNET_H_
#define BMRNET_H_

#ifdef _WIN32
 #include<winsock2.h>
#else
	#include <sys/socket.h>
	#include <netinet/in.h>
	#include <arpa/inet.h>
	#include <stdbool.h>
#endif
/*GLOBAL VARIABLES - LIST OF IP ADDRESSES*/
extern char** localIPaddrs;
extern int numberOfAddresses;
#define NUMCONNECTIONS 3


//gets the list of IP addresses
char** getIPAddresses();
int getPartyNum(char* filename);


class CommunicationObject
{
private:
	uint64_t bytesSent = 0;
	uint64_t bytesReceived = 0;
	uint64_t numberOfSends = 0;
	uint64_t numberOfRecvs = 0;	
	std::array<uint64_t, NUM_MEASUREMENT_PHASES> phaseBytesSent;
	std::array<uint64_t, NUM_MEASUREMENT_PHASES> phaseBytesReceived;
	std::array<uint64_t, NUM_MEASUREMENT_PHASES> phaseSends;
	std::array<uint64_t, NUM_MEASUREMENT_PHASES> phaseRecvs;
	MeasurementPhase phase = ONLINE_PHASE;
	bool measurement = false;	

public: 
	CommunicationObject()
	{
		reset();
	}

	void reset()
	{
		bytesSent = 0;
		bytesReceived = 0;
		numberOfSends = 0;
		numberOfRecvs = 0;
		phaseBytesSent.fill(0);
		phaseBytesReceived.fill(0);
		phaseSends.fill(0);
		phaseRecvs.fill(0);
		phase = ONLINE_PHASE;
		measurement = false;
	}

	void setMeasurement(bool a)
	{
		measurement = a;
	}

	void setPhase(MeasurementPhase p)
	{
		phase = p;
	}

	void incrementSent(int size)
	{
		if (measurement)
		{
			bytesSent += size;
			numberOfSends++;
			phaseBytesSent[phase] += size;
			phaseSends[phase]++;
		}
	}

	void incrementRecv(int size)
	{
		if (measurement)
		{
			bytesReceived += size;
			numberOfRecvs++;
			phaseBytesReceived[phase] += size;
			phaseRecvs[phase]++;
		}
	}

	uint64_t getSent() {return bytesSent;}
	uint64_t getRecv() {return bytesReceived;}
	uint64_t getRoundsSent() {return numberOfSends;}
	uint64_t getRoundsRecv() {return numberOfRecvs;}
	uint64_t getSent(MeasurementPhase p) {return phaseBytesSent[p];}
	uint64_t getRecv(MeasurementPhase p) {return phaseBytesReceived[p];}
	uint64_t getRoundsSent(MeasurementPhase p) {return phaseSends[p];}
	uint64_t getRoundsRecv(MeasurementPhase p) {return phaseRecvs[p];}
	bool getMeasurement() {return measurement;}
};


class BmrNet {
private:
	char * host;
	unsigned int port;
	bool is_JustServer;
	int socketFd[NUMCONNECTIONS];
	#ifdef _WIN32
	    PCSTR Cport;
		WSADATA wsa;
		DWORD dwRetval;
	#endif


public:
	/**
	 * Constructor for servers and clients, got the host and the port for connect or listen.
	 * After creation call listenNow() or connectNow() function.
	 */
	BmrNet(char * host, int port);

	/**
	 * Constructor for servers only. got the port it will listen to.
	 * After creation call listenNow() function.
	 */
	BmrNet(int portno);

	/**
	 * got data and send it to the other side, wait for response and return it.
	 * return pointer for the data that recived.
	 */
	void* sendAndRecive(const void* data, int get_size, int send_size);

	
	virtual ~BmrNet();

	/**
	 * Start listen on the given port.
	 */
	bool listenNow();

	/**
	 * Try to connect to server by given host and port.
	 * return true for success or false for failure.
	 */
	bool connectNow();

	/**
	 * Send Data to the other side.
	 * return true for success or false for failure.
	 */
	bool sendMsg(const void* data, int size, int conn);

	/**
	 * Recive data from other side.
	 * return true for success or false for failure.
	 */
	bool receiveMsg(void* buff, int buffSize, int conn);


};



#endif /* BMRNET_H_ */
