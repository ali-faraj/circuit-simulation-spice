//Commande de compilation "g++ -I D:\Documents\eigen-3.4.0 Device.hpp Circuit.hpp Circuit.cpp main.cpp -o main.exe"
#include <iostream>
#include "Circuit.hpp"

using namespace std;
using std::vector;

int main(int argc, char * argv[])
{
	Circuit circ;
	if (argc < 2){
		cout << "Input file name is required when calling executable main" << endl;
	}
	else{
		circ.openLogFile("logfile.log");
		circ.buildCirc(argv[1]);
		//circ.display();
		circ.buildMat();
		//circ.dispMat();
		circ.simulate("waveformfile.cvs");
		circ.closeLogFile();
	}
	return 0;
}
