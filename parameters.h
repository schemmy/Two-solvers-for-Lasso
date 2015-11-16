#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>

using namespace std;

int parseDistributedOptions(char* &fileName, int &method, double &tau, double &lambda,
		int argc, char *argv[]) {

	char c;

	while ((c = getopt(argc, argv, "A:M:t:l:")) != -1) {
		switch (c) {
		case 'A':
			fileName = optarg;
			break;

		case 'M':
			method = atoi(optarg);
			break;

		case 't':
			tau = atof(optarg);
			break;

		case 'l':
			lambda = atof(optarg);
			break;

		default:
			break;
		}

	}

	return 0;
}

#endif