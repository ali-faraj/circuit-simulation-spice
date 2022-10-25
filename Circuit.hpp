#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <iostream>
#include <fstream> 
#include <vector>
#include "Device.hpp"
#include <Eigen/Dense>

using namespace std;
using std::vector;
using Eigen::MatrixXd;

class Circuit{
private :
	unsigned int numRes, numInd, numCap, numConstVsource, numConstIsource, numVsource, numIsource, n, m;
	double tstop, timeConst;
	bool trans, NonNegPole;
	vector<ConstElem> resistor, inductor, capacitor; //Vectors of passive elements
	vector<ConstElem> constVsource, constIsource; //Vectors of constant voltage and intensity sources
	vector<TimeDepSource> Vsource, Isource; //Vectors of time-dependent voltage and intensity sources
	MatrixXd A, D; //Matrices of the algebro-differential equation
	MatrixXd I, Iconst, E; //Vectors corresponding to intensity and voltage sources (Iconst is the constant part of the intensity source)
	MatrixXd Z; //Right-hand-side of the discretized system
	std::ofstream logFile;
	
public :
	
	void display();
	
	void buildCirc(const char *inputFileName);
	
	void openLogFile(const char *logFileName);
	
	void closeLogFile();
	
	unsigned int getNumNode();
	
	void buildMat();
	
	void dispMat();

	void simulate(const char *waveformFileName);

	void computeParam(unsigned int &nt, double &h);

	void computePole();
	
	Circuit();
	
};

#endif
