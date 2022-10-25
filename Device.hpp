#ifndef DEVICE_H
#define DEVICE_H

#include <iostream>
#include <vector>

using namespace std;
using std::vector;

//Super-class from which will inherites the circuit devices
class Device{
private :
	unsigned int n1, n2;
public :
        unsigned int getN1();
	
        unsigned int getN2();
	
        void display();
	
        Device();
	
        Device(unsigned int node1, unsigned int node2);
	
};

//Constant device (resitor, capacitor, inductor, constant source)
class ConstElem : public Device{
private :
        double val;
public :
        double getVal();
	
        void display();
	
        ConstElem();
	
        ConstElem(unsigned int node1, unsigned int node2, double value);

	ConstElem(unsigned int node1, unsigned int node2, string word);
	
};

//Time dependent source (voltage or intensity)
class TimeDepSource : public Device{
private :
	unsigned int s;
        vector<double> t, val;
public :
        const vector<double> &getVal();
	
        const vector<double> &getT();
	
        void display();

	void computeNumParam(double &h, double &tmax);

	double eval(double t0);
	
        TimeDepSource();
	
        TimeDepSource(unsigned int node1, unsigned int node2, const vector<double> &tm, const vector<double> &value);

	TimeDepSource(unsigned int node1, unsigned int node2, const vector<string> &lineWord, unsigned int numWord);
	
};

#endif
