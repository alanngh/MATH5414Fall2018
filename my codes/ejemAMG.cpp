#include "aggragation.h"


using namespace std;



int main()
{
	int n=9;
	matrix A(n,n);

	for (int i=0;i<n;i++)	A.array[i][i] = 4;
	for (int i=1;i<n;i++){
		A.array[i][i-1] = -1;
		A.array[i-1][i] = -1;
	}
	for (int i=3;i<n;i++){
		A.array[i][i-3] = -1;
		A.array[i-3][i] = -1;
	}

	A.array[2][3]=0;  	
	A.array[3][2]=0;
	A.array[5][6]=0;
	A.array[6][5]=0;

	A.print();


	Aggragation(A,0.2);

	return 0;
}
