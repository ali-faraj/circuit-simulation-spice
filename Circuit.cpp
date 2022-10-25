#include <iostream>
#include <algorithm>
#include "Circuit.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <math.h> 
#include <Eigen/Eigenvalues>

using namespace std;

void Circuit::display(){
	cout << "** Circuit:" << endl;
	for (unsigned int i = 0; i < numRes; i++){
                cout << "R" << i << " ";
                resistor[i].display();
	}
	for (unsigned int i = 0; i < numInd; i++){
                cout << "L" << i << " ";
                inductor[i].display();
	}
	for (unsigned int i = 0; i < numCap; i++){
                cout << "C" << i << " ";
                capacitor[i].display();
	}
	for (unsigned int i = 0; i < numConstVsource; i++){
                cout << "V" << i << " ";
                constVsource[i].display();
	}
	for (unsigned int i = 0; i < numConstIsource; i++){
                cout << "I" << i << " ";
                constIsource[i].display();
	}
	for (unsigned int i = 0; i < numVsource; i++){
		cout << "V" << i << " ";
                Vsource[i].display();
	}
	for (unsigned int i = 0; i < numIsource; i++){
                cout << "I" << i << " ";
                Isource[i].display();
	}
	logFile << "* Display circuit proceded" << endl << endl;
}

//Read the devices form input file and build the circuit
void Circuit::buildCirc(const char *inputFileName){
	ifstream inputFile(inputFileName); //input stream corresponding to input file (with name inputFileName)
	string line, word;
	vector<string> lineWord; //vector of the words in one line 
	unsigned int n1, n2;
	unsigned int numWord; //number of words in one line
	char firstChar; //first character of the line
	istringstream issLine; //input string stream corresponding to one line
	if (inputFile){ //if the input file is open
		while (getline(inputFile, line)){ //read the lines and and attribute to variable line
			firstChar = line.at(0);
			issLine.clear();
			lineWord.clear();
			issLine.str(line);
			//Get the words and add them to the string vector lineWord
			while( getline(issLine, word, ' ') ){
				lineWord.push_back(word);
			}
			//number of words in the line
			numWord = lineWord.size();
			if (numWord < 3){
				cout << "wrong format in input file 1" << endl;
				exit(1);
			}
			//Indicate if the simulation is transient and reads tstop in the transient case
			if (lineWord[0] == ".tran"){
				tstop = atof( lineWord[2].c_str() );
				trans = true;
			}
			else{
				//get the nodes of the device
				n1 = stoi(lineWord[1]);
				n2 = stoi(lineWord[2]);
				//update the number of nodes
				n = max(n, max(n1, n2));
				switch(firstChar) {
					//Add a device depending on the first character and increase the corresponding number of devices
				case 'R':
					resistor.push_back( ConstElem( n1, n2, lineWord[3] ) );
					numRes = numRes+1;
					break;
				case 'L':
					inductor.push_back( ConstElem( n1, n2, lineWord[3] ) );
					numInd = numInd+1;
					break;
				case 'C':
					capacitor.push_back( ConstElem( n1, n2, lineWord[3] ) );
					numCap = numCap+1;
					break;
					//Add a voltage source by removing useless strings from input file and increase the corresponding number of sources
				case 'V':
					if (lineWord[3].compare(0, 3, "dc=") == 0){
						word = lineWord[3].substr(3, lineWord[3].size()-3);
						constVsource.push_back( ConstElem( n1, n2, word ) );
						numConstVsource = numConstVsource+1;
					}
					//if the source is time dependent, the source times and values are obtained by reading the line
					else if (lineWord[3].compare(0, 4, "PWL(") == 0){
						Vsource.push_back(TimeDepSource(n1, n2, lineWord, numWord));
						numVsource = numVsource+1;
					}
					else{
						cout << "wrong format in input file 2" << endl;
						exit(1);
					}
					break;
					//Add an intensity source by removing useless strings from input file and increase the corresponding number of sources
				case 'I':
					if (lineWord[3].compare(0, 3, "dc=") == 0){
						word = lineWord[3].substr(3, lineWord[3].size()-3);
						constIsource.push_back( ConstElem( n1, n2, word ) );
						numConstIsource = numConstIsource+1;
					}
					//if the source is time dependent, the source times and values are obtained by reading the line
					else if (lineWord[3].compare(0, 4, "PWL(") == 0){
						Isource.push_back(TimeDepSource(n1, n2, lineWord, numWord));
						numIsource = numIsource+1;
					}
					else{
						cout << "wrong format in input file 3" << endl;
						exit(1);
					}
					break;
				default:
					cout << "wrong format in input file 4" << endl;
					cout << line << endl;
					exit(1);
				}
			}
		}
		logFile << "* Input data reading complete. Circuit is built." << endl << endl;
		inputFile.close();
	}
	else{
		cout << "Could not open input file" << endl;
		exit(1);
	}
}

void Circuit::openLogFile(const char *logFileName){
	logFile.open(logFileName);
	if (!logFile.is_open()){
		cout << "Could not open log file" << endl;
		exit(1);
	}
}

void Circuit::closeLogFile(){
	logFile.close();
}

unsigned int Circuit::getNumNode(){
	return n;
}

//Build the matrices of the algebro-differential equation D*dX/dt+AX=Z
void Circuit::buildMat(){
	int n1, n2;
	double val;
	MatrixXd G, B; //Block matrices for the matrix A
	MatrixXd C, L; //Block matrices for the matrix D
	//Size of matrix blocks corresponding to voltage sources and inductors (equal to the number of intensities to be computed)
	m = numConstVsource + numVsource + numInd;
	//Preallocate full arrays
	A.resize(n+m,n+m);
	A = MatrixXd::Zero(n+m,n+m);
	D.resize(n+m,n+m);
	D = MatrixXd::Zero(n+m,n+m);
	Z.resize(n+m,1);
	Z = MatrixXd::Zero(n+m,1);
	//Preallocate the first blocks
	G.resize(n,n);
	G = MatrixXd::Zero(n,n);
	C.resize(n,n);
	C = MatrixXd::Zero(n,n);
	Iconst.resize(n,1);
	Iconst = MatrixXd::Zero(n,1);
	I.resize(n,1);
	I = MatrixXd::Zero(n,1);
	//Fill the matrix G using resistors
	for (unsigned int i = 0; i < numRes; i++){
		n1 = resistor[i].getN1() - 1;
		n2 = resistor[i].getN2() - 1;
		val = resistor[i].getVal();
		if (n1 >= 0){
			G(n1,n1) = G(n1,n1) + 1./val;
		}
		if (n2 >= 0){
			G(n2,n2) = G(n2,n2) + 1./val;
		}
		if (n1 >= 0 && n2 >= 0){
			G(n1,n2) = G(n1,n2) - 1./val;
			G(n2,n1) = G(n2,n1) - 1./val;
		}
	}
	//Fill the matrix C using capacitors
	for (unsigned int i = 0; i < numCap; i++){
		n1 = capacitor[i].getN1() - 1;
		n2 = capacitor[i].getN2() - 1;
		val = capacitor[i].getVal();
		if (n1 >= 0){
			C(n1,n1) = C(n1,n1) + val;
		}
		if (n2 >= 0){
			C(n2,n2) = C(n2,n2) + val;
		}
		if (n1 >= 0 && n2 >= 0){
			C(n1,n2) = C(n1,n2) - val;
			C(n2,n1) = C(n2,n1) - val;
		}
	}
	//Fill the matrix I using intensity sources
	for (unsigned i = 0; i < numConstIsource; i++){ //Constant intensity sources
		val = constIsource[i].getVal();
		n1 = constIsource[i].getN1() - 1;
		n2 = constIsource[i].getN2() - 1;
		if (n1 >= 0){
			Iconst(n1) = Iconst(n1) - val;
		}
		if (n2 >= 0){
			Iconst(n2) = Iconst(n2) + val;
		}
	}
	I = Iconst;
	for (unsigned i = 0; i < numIsource; i++){ //Time-dependent intensity sources
		val = Isource[i].eval(0.);
		n1 = Isource[i].getN1() - 1;
		n2 = Isource[i].getN2() - 1;
		if (n1 >= 0){
			I(n1) = I(n1) - val;
		}
		if (n2 >= 0){
			I(n2) = I(n2) + val;
		}
	}
	//Replace first blocks in the full array
	A.block(0,0,n,n) = G;
	D.block(0,0,n,n) = C;
	Z.block(0,0,n,1) = I;
	//////////////////////
	if (m != 0){
		//Preallocate the other blocks
		B.resize(n,m);
		B = MatrixXd::Zero(n,m);
		E.resize(m,1);
		E = MatrixXd::Zero(m,1);
		L.resize(m,m);
		L = MatrixXd::Zero(m,m);
		//Fill the matrix B and E corresponding to voltage sources
		for (unsigned i = 0; i < numConstVsource; i++){//Constant voltage sources
			n1 = constVsource[i].getN1() - 1;
			n2 = constVsource[i].getN2() - 1;
			if (n1 >= 0){
				B(n1,i) = 1.;
			}
			if (n2 >= 0){
				B(n2,i) = -1.;
			}
			E(i) = constVsource[i].getVal();
		}
		for (unsigned i = 0; i < numVsource; i++){//Time-dependent voltage sources
			n1 = Vsource[i].getN1() - 1;
			n2 = Vsource[i].getN2() - 1;
			if (n1 >= 0){
				B(n1,numConstVsource+i) = 1.;
			}
			if (n2 >= 0){
				B(n2,numConstVsource+i) = -1.;
			}
			E(numConstVsource+i) = Vsource[i].eval(0.);
		}
		//Fill the matrix B and L corresponding to inductors
		for (unsigned i = 0; i < numInd; i++){
			n1 = inductor[i].getN1() - 1;
			n2 = inductor[i].getN2() - 1;
			if (n1 >= 0){
				B(n1,numConstVsource+numVsource+i) = 1.;
			}
			if (n2 >= 0){
				B(n2,numConstVsource+numVsource+i) = -1.;
			}
			L(numConstVsource+numVsource+i,numConstVsource+numVsource+i) = -inductor[i].getVal();
		}
		//Replace other blocks in the full array
		A.block(0,n,n,m) = B;
		A.block(n,0,m,n) = B.transpose();
		D.block(n,n,m,m) = L;
		Z.block(n,0,m,1) = E;	
	}
	if (abs(A.determinant()) < 1.e-14){
		cout << "Matrix A is singular" << endl;
		exit(1);
	}
	if (numCap > 0 || numInd > 0){ //If there is a time derivative in the algebro-differential equation
		//Computation of the time constant of the system
		computePole();
		logFile << "* Time constant of the system computed" << endl << endl;
		if (!NonNegPole){
			logFile << "** All the poles have a negative real part **" << endl << endl;	
		}
	}
	logFile << "* Matrix A, D and Z completed" << endl <<endl;
}

//Computes the voltages and intensities (steady or time evolution) solution to D*dX/dt+AX=Z and writes the output in a waveform file
void Circuit::simulate(const char *waveformFileName){
	int n1, n2;
	unsigned int nt; //number of time steps
	double h; //time step
	double val;
	vector<double> t;
	MatrixXd M; //Matrix corresponding to the discrete time evolution
	MatrixXd X0; //Steady solution or solution at previous time step
	MatrixXd X1; //Solution at local time step
	MatrixXd X; //Matrix containing the solution for all time steps
	std::ofstream waveformFile;
	//Open waveform file
	waveformFile.open(waveformFileName);
	if ( !waveformFile.is_open() ){
		cout << "Could not open waveform file" << endl;
		exit(1);
	}
	//Compute solution X0 (initial for transient or steady-state for time-independent)
	X0.resize(n+m,1);
	X0 = MatrixXd::Zero(n+m,1);
	X0 = A.colPivHouseholderQr().solve(Z);
	if (trans){//In the case of the time-dependent simulation
		logFile << "** Simulation is time-dependent **" << endl << endl;
		logFile << "* Initial solution computed" << endl << endl;
		//Compute the numerical parameters
		computeParam(nt,h);
		logFile << "* Numerical parameters computed" << endl << endl;
		//Preallocate full arrays
		M.resize(n+m,n+m);
		M = MatrixXd::Zero(n+m,n+m);
		X1.resize(n+m,1);
		X1 = MatrixXd::Zero(n+m,1);
		X.resize(n+m,nt);
		X = MatrixXd::Zero(n+m,nt);
		//Compute the matrix M of the discretized evolution
		M = 1./h*D + A;
		if (abs(M.determinant()) < 1.e-14 && nt > 1){
			cout << "Matrix M is singular" << endl;
			exit(1);
		}
		//Add initial solution
		X.block(0,0,n+m,1) = X0;
		t.push_back(0.);
		//Time loop
		for (unsigned int k = 1; k < nt; k++){
			// Add time value corresponding to time step
			t.push_back(k*h);
			// Update intensity sources (only the time dependent ones)
			I = Iconst;
			for (unsigned int i = 0; i < numIsource; i++){
				val = Isource[i].eval(t[k]);
				n1 = Isource[i].getN1() - 1;
				n2 = Isource[i].getN2() - 1;
				if (n1 >= 0){
					I(n1) = I(n1) - val;
				}
				if (n2 >= 0){
					I(n2) = I(n2) + val;
				}
			}
			Z.block(0,0,n,1) = I;
			if (m!= 0){//If there are voltage sources or inductors
				for (unsigned int i = 0; i < numVsource; i++){
					n1 = Vsource[i].getN1() - 1;
					n2 = Vsource[i].getN2() - 1;
					E(numConstVsource+i) = Vsource[i].eval(t[k]);
				}
				Z.block(n,0,m,1) = E;
			}
			//Udate the right hand side and solve the linear system corresponding to the discrete evolution
			Z = Z + 1./h*D*X0;
			X1 = M.colPivHouseholderQr().solve(Z);
			//Update solution at previous time step
			X0 = X1;
			//Add local solution
			X.block(0,k,n+m,1) = X1;
		}
		logFile << "* Transient simulation is finished" << endl << endl;
		if (NonNegPole){
			logFile << "**** Warning: Time evolution of steady-state solution is not necessarly meaningfull (frequency analysis might be required) ****" << endl << endl;
		}
		//Write solution in waveform file (computed voltages and intensities at all time steps)
		waveformFile << ".HEADER" << endl;
		waveformFile << "Transient simulation" << endl;
		//Write voltages
		for (unsigned int i = 0; i < n; i++){
			waveformFile << "..NAMES" << endl;
			waveformFile << "Time,v(" << i + 1 << ")" << endl;
			waveformFile << "..UNITS" << endl;
			waveformFile << "s,Voltage(V)" << endl;
			waveformFile << ".DATA" << endl;
			//For all time values
			for (unsigned int k = 0; k < nt; k++){
				waveformFile << t[k] << "," << X(i,k) << endl;
			}
		}
		//Write intensities
		for (unsigned int i = n; i < n+m; i++){
			waveformFile << "..NAMES" << endl;
			waveformFile << "Time,i(" << i + 1 - n << ")" << endl;
			waveformFile << "..UNITS" << endl;
			waveformFile << "s,Current(A)" << endl;
			waveformFile << ".DATA" << endl;
			//For all time values
			for (unsigned int k = 0; k < nt; k++){
				waveformFile << t[k] << "," << X(i,k) << endl;
			}
		}
	}
	else{//correspond to time-independent simulation
		logFile << "** Simulation is time-independent **" << endl << endl;
		logFile << "* Steady-state solution computed" << endl << endl;
		if (NonNegPole){
			logFile << "**** Warning: Steady-state solution is not necessarly meaningfull (frequency analysis might be required) ****" << endl << endl;
		}
		//Write steady state solution in waveform file (voltages and intensities)
		waveformFile << ".HEADER" << endl;
		waveformFile << "Time independent" << endl;
		//Write voltages
		for (unsigned int i = 0; i < n; i++){
			waveformFile << "..NAMES" << endl;
			waveformFile << "Time,v(" << i + 1 << ")" << endl;
			waveformFile << "..UNITS" << endl;
			waveformFile << "s,Voltage(V)" << endl;
			waveformFile << ".DATA" << endl;
			waveformFile << 0 << "," << X0(i) << endl;
		}
		//Write intensities
		for (unsigned int i = n; i < n+m; i++){
			waveformFile << "..NAMES" << endl;
			waveformFile << "Time,i(" << i + 1 - n << ")" << endl;
			waveformFile << "..UNITS" << endl;
			waveformFile << "s,Current(A)" << endl;
			waveformFile << ".DATA" << endl;
			waveformFile << 0 << "," << X0(i) << endl;
		}
	}
	waveformFile.close();
	logFile << "* Output written in waveform file" << endl << endl;
}

//Computation of the time step and number of time steps for the simulator
void Circuit::computeParam(unsigned int &nt, double &h){
	double tmax; //maximal time corresponding to a time dependent source
	
	if (tstop < 0.){
		cout << "Negative values are not accepted for tstop" <<endl;
		exit(1);
	}
	else if (tstop == 0.){ //Only the initial solution is computed
		nt = 1;
		h = 1.;
	}
	else{
		//Initialize the time step and maximal time corresponding to time dependent sources (if any)
		tmax = 1.e12;
		h = tstop/100.;
		//Update the time step and maximal time corresponding to time dependent sources (if any)
		for (unsigned int i = 0; i < numVsource; i++){
			Vsource[i].computeNumParam(h,tmax);
		}
		for (unsigned int i = 0; i < numIsource; i++){
			Isource[i].computeNumParam(h,tmax);
		}
		if (tstop > tmax){
			cout << "tstop must be smaller then last input time" << endl;
			exit(1);
		}
		if (numInd > 0 || numCap > 0){//if there is a time derivative the time step is refined using the time constant of the system
			h = min(h,timeConst/20.);
		}
		//compute the number of time steps using the time step h and the final simulation time
		nt = floor(tstop/h)+1;
	}
}

//If there are time derivatives in the algebro-differential equation, computes the time constant of the system and examines the presence of poles with non negative real part
void Circuit::computePole(){
	double absval;
	Eigen::EigenSolver<MatrixXd> es(A.inverse()*D); //Eigenvalues of A^{-1}*D
	for (unsigned i = 0; i < n+m; i++){
		absval = abs(es.eigenvalues()(i)); //Compute the inverse of the absolute value of the poles
		if (absval > 1.e-14){
			timeConst = min( timeConst , 1./absval ); // update the time constant using the absolute value of the poles
			if ( real( -1./es.eigenvalues()(i) ) > -1.e-14 ){//if the real part of a pole is non negative
				NonNegPole = true;
			}
		}
	}
}

//Display the matrices
void Circuit::dispMat(){
	cout << "** Matrix A" << endl;
	cout  <<  A  <<  endl;
	cout << "** Matrix D" << endl;
	cout << D << endl;
	cout << "** Matrix Z" << endl;
	cout << Z << endl;
	logFile << "* Display matrices proceded" << endl << endl;
}

//Default constructor
Circuit::Circuit(){
	numRes = 0;//Number of resistors
	numInd = 0;//Number of inductors
	numCap = 0;//Number of capacitors
	numConstVsource = 0;//Number of constant voltage sources
	numConstIsource = 0;//Number of constant intensity sources
	numVsource = 0;//Number of time-dependent voltage sources
	numIsource = 0;//Number of time dependent intensity sources
	n = 0;//Number of nodes
	m = 0;//Number of intensities corresponding to voltage sources and inductors
	trans = false;//Indicate if the simulation is transient
	NonNegPole = false;//Indicates if there are poles with nonnegative Real part (this variable is used only if there is a time derivative in the algebro-differential equation)
	tstop = 0.;//Final simulation time
	timeConst = 1.e12;//Time constant of the system (this variable is used only if there is a time derivative in the algebro-differential equation)
}
