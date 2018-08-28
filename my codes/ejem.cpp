#include <iostream>
#include "matrix.h"

using namespace std;



int main()
{


// testing the matrix class and basic operation + * and inverse
	int r = 4;
	int c = 4;
	matrix matrixA(r,c);
	matrix matrixB(r,c);
	matrix matrixC(r,c);
	matrix matrixD(r,c);

	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			matrixA.setElement(i,j,1);
			matrixB.setElement(i,j,2);
		}
	}
	
/*	cout << "matrixA" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixA.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;

	cout << "matrixB" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixB.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;
	
	cout << "matrixC before addition" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixC.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;	
	
	matrixC = matrixA + matrixB;
	
	cout << "matrixC after addition" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixC.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;	


	matrixC = matrixA * matrixB;
	
	cout << "matrixC after multiplication" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixC.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;	*/


matrixD.setElement(0,0,5);
matrixD.setElement(0,1,6);
matrixD.setElement(0,2,6);
matrixD.setElement(0,3,8);


matrixD.setElement(1,0,2);
matrixD.setElement(1,1,2);
matrixD.setElement(1,2,2);
matrixD.setElement(1,3,8);

matrixD.setElement(2,0,6);
matrixD.setElement(2,1,6);
matrixD.setElement(2,2,2);
matrixD.setElement(2,3,8);

matrixD.setElement(3,0,2);
matrixD.setElement(3,1,3);
matrixD.setElement(3,2,6);
matrixD.setElement(3,3,7);


	cout << "matrixD" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixD.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;

	matrixC = matrixD.inv();

	cout << "matrixD after inv" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixC.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;	


	matrixA = matrixC * matrixD;
	
	cout << "matrixC after multiplication" << endl;
	for (int i=0; i<r; i++){
		for (int j=0; j<c; j++){
			cout << matrixA.getElement(i,j) << "\t";
		}
		cout << endl;
	}
	cout << endl;	








}
