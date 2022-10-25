#include <iostream>
#include <vector>
#include "Device.hpp"

using namespace std;
using std::vector;

unsigned int Device::getN1(){
	return n1;
}

unsigned int Device::getN2(){
	return n2;
}

void Device::display(){
	cout << n1 << " " << n2 << " ";
}

Device::Device(){ //Default constructor
	n1 = 0; //node 1
	n2 = 0; //node 2
}

Device::Device(unsigned int node1, unsigned int node2){ //Constructor
	n1 = node1;
	n2 = node2;
}

double ConstElem::getVal(){
	return val;
}

void ConstElem::display(){
	Device::display();
	cout << " " << val << endl;
}

ConstElem::ConstElem():Device(){ //Default constructor
	val = 0.; //Value
}

ConstElem::ConstElem(unsigned int node1, unsigned int node2, double value):Device(node1, node2){ //Constructor
	val = value;
}

ConstElem::ConstElem(unsigned int node1, unsigned int node2, string word):Device(node1, node2){ //Constructor using a word read from input file
	//Convert word to a double and set the value val
	val = atof( word.c_str() );
}

const vector<double> &TimeDepSource::getVal(){
	return val;
}

const vector<double> &TimeDepSource::getT(){
	return t;
}

void TimeDepSource::display(){
	Device::display();
	cout << "time(";
	for (unsigned int i = 0; i < t.size(); i++){
		cout << " " << t[i];
	}
	cout << " ) value(";
	for (unsigned int i = 0; i < val.size(); i++){
		cout << " " << val[i];
	}
	cout << " )" << endl;
}

//Compute the time step and final time corresponding to a time dependent source
void TimeDepSource::computeNumParam(double &h, double &tmax){
	if (s == 1){
		tmax = min(tmax,t[0]);
	}
	else{
		for (unsigned int k = 0; k < s-1; k++){
			if (t[k+1]-t[k] <= 0.){
				cout << "Input time must but in increasing order" << endl;
				exit(1);
			}
			else{
				h = min(h, t[k+1]-t[k]);
			}
		}
		tmax=min(tmax, t[s-1]);
	}
}

//Gives the value of a source corresponding to a time t0 (constant interpolation is used)
double TimeDepSource::eval(double t0){
	double v0=0.;
	if(s == 1 && t0 <= t[0]){
		v0 = val[0];
	}
	else if (t0 == t[s-1]){
		v0 = val[s-1];
	}
	else{
		for (unsigned int k=0; k < s-1; k++){
			if (t0 >= t[k] && t0 < t[k+1]){
				v0 = val[k];
				break;
			}
		}
	}
	return v0;
}

TimeDepSource::TimeDepSource():Device(){ //Default constructor
	s = 0; //Number of values of time dependent source
}

TimeDepSource::TimeDepSource(unsigned int node1, unsigned int node2, const vector<double> &tm, const vector<double> &value):Device(node1, node2){ //Constructor
	t = tm;
	val = value;
	if (tm.size() != value.size()){
		cout << "input TimeDepSource not correct form" << endl;
		exit(1);
		//return 1;
	}
	else{
		s = t.size();
	}
}

TimeDepSource::TimeDepSource(unsigned int node1, unsigned int node2, const vector<string> &lineWord,unsigned int numWord):Device(node1, node2){ //Constructor using the string read from the inputfile
	int lastWordSize;
	string word; //Temporary string
	t.clear();
	val.clear();
	//removes the useless "PWL(" 
	word = lineWord[3].substr(4, lineWord[3].size()-4);
	//convert words read from input file to double and add to the source times and values
	t.push_back( atof( word.c_str() ) );
	val.push_back( atof( lineWord[4].c_str() )  );
	for (unsigned int i = 5; i < numWord-3; i = i+2){
		t.push_back( atof( lineWord[i].c_str() ) );
		val.push_back( atof( lineWord[i+1].c_str() ) );
	}
	t.push_back( atof( lineWord[numWord-2].c_str() ) );
	lastWordSize = lineWord[numWord-1].size();
	word = lineWord[numWord-1].substr(0, lastWordSize-1);//remove useless ")"
	val.push_back( atof( word.c_str() ) );
	//set the size of the time and value vectors 
	s = t.size();
}
