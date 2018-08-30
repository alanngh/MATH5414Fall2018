#include <stdio.h> 
#include <stdlib.h>
#include <math.h> 
#include "matrix.h"




int Find(int* N_i, int ind , int key ){
	for (int i=0;i<ind;i++){
		if (N_i[i]==key)	return i;
	}
	return -1;
}





//void Aggagation(matrix A, double eps, matrix &AF, matrix &PF, matrix &P){
void Aggragation(matrix A, double eps){	

	int **N = new int*[A.row];        
	int N_index[A.row];
	int V[A.row];
	int k;

	// compute strong conections
	for (int i=0;i<A.row;i++){
		k=0;
		for (int j=0;j<A.row;j++){
			if (fabs(A.array[i][j]) >= eps*sqrt(A.array[i][i]*A.array[j][j])){
				V[k] = j;
				k++;				
			}
		}
		N[i] = new int[k];
		for (int l=0;l<k;l++)      N[i][l]=V[l];
		N_index[i] = k;
	}


/*      for (int i=0;i<A.row;i++){
		printf("\n strongly coupled neighborhood of node %d\n",i);
		for (int j=0;j<N_index[i];j++)  	printf("%d   ",N[i][j]);
		printf("\n");

	}  */


	// compute filtered matrix
	printf("\nFiltered Matrix\n");
	matrix AF(A.row,A.col);

	for (int i=0;i<A.row;i++){
		for (int j=0;j<A.col;j++){
			if (j!=i){
				if (Find(N[i],N_index[i],j)!=-1)	AF.array[i][j] = A.array[i][j];
			}
		}	
	
	
		double S = 0;
		for (int j=0;j<A.col;j++){
			if (j!=i)	S = A.array[i][j] - AF.array[i][j];
		}

		AF.array[i][i] = A.array[i][i] - S;
	}

	AF.print();


	// compute R (nonaisled nodes?)


	int R[A.row];
	int R_index = 0;

	for (int i=0 ; i<A.row ; i++){
		if (N_index[i] > 1){  //index does not start at zero
			R[R_index] = i;
			R_index++;
		}
	}

	k = R_index;
	int mA = 0;
	
	int **C = new int*[A.row];        
	int C_index[A.row];


	for (int i=0;i<R_index;i++){   // select the neigborhood

		bool ok = true;

		for (int j=0;j<N_index[i];j++){ // select elements in the neigborhood i 
			bool Enc = false;
			for (int l=0;l<k;l++){ // compare with each element in R
				if (N[i][j] == R[l] ) 	Enc = true;
				if (Enc == false && l==k-1 ){
					ok = false;
					break;
				}
			}			
			if (ok == false) break;
		}


		//printf(" prueba ok = %d \n",ok);

		if (ok){
			C[mA] = new int[N_index[i]];  // create each component
			for (int j=0;j<N_index[i];j++){
				C[mA][j] = N[i][j];
				C_index[mA] = N_index[i]; 
			}
			
		/*	// imprimir C
			for (int l=0;l<C_index[mA];l++){
				printf("%d\t",C[mA][l]);
			}
			printf("\n"); */


			// remove C from R
			for (int j = 0;j<C_index[mA];j++){
				for (int l =0;l<R_index;l++){
					if (C[mA][j] == R[l]){
						R[l] = -1;
						break;
					}
				}
			} 

			/*// imprimir R
			for (int l=0;l<R_index;l++){
				printf("%d\t",R[l]);
			}
			printf("\n"); */

			R_index = R_index - C_index[mA];  
// 			printf("R_index = %d ,  C_index[%d] = %d  ,  New R_index = %d  \n"  ,k,mA,C_index[mA],R_index);

			int j=0;
			while (j < R_index){  //delet C from R
				if (R[j]==-1){
					for (int l=j+1;l<k;l++){
						if (R[l]!= -1){
							R[j] = 	R[l];
							R[l] = -1;
							break;
						}
					}		
				}
				j++;					
			}

			/* // imprimir R
			for (int l=0;l<R_index;l++){
				printf("%d\t",R[l]);
			}
			printf("\n");*/

			k = R_index;
			mA++; 
		}
		//else printf("\n neigborhood %d is no in R\n",i);
	}
	
	// print C initials
	printf("\nInitial aggregates : \n");
	for (int i =0;i<mA;i++){
		for (int j=0;j<C_index[i];j++){
			printf("%d\t",C[i][j]);
		}
		printf("\n");
	}

	printf("\nResidual nodes  R_index = %d \n",k);
	for (int i=0;i<R_index;i++){
		printf("%d\t",R[i]);
	}
	printf("\n");

	for (int i=0;i<k;i++){
		int Nod;
		for (int j=0;j < N_index[R[i]]; j++){
			if (R[i] != N[R[i]][j]){
				Nod = N[R[i]][j];
				break;
			}
		}


		for (int j=0;j<mA;j++){
			int P = Find(C[j],C_index[j],Nod);
			if (P > -1){
				int temp[C_index[j]+1];  // create new component
				for (int l=0;l<C_index[j];l++){
					temp[l] = C[j][l];
				}
				temp[C_index[j]] = R[i];

				delete C[j];   // this is not the best but works
				C_index[j] = C_index[j]+1;
				C[j] = new int[C_index[j]];
				for (int l=0;l<C_index[j];l++){
					C[j][l] = temp[l];
				}

			}
		}
	}

	// print C final
	printf("\nFinal aggregates : \n");
	for (int i =0;i<mA;i++){
		for (int j=0;j<C_index[i];j++){
			printf("%d\t",C[i][j]);
		}
		printf("\n");
	}




	// Time to construct P
	matrix P(A.row,mA);
	for (int i=0;i<mA;i++){
		for (int j=0;j<C_index[i];j++){
			P.array[C[i][j]][i] = 1;
		}		
	}

	P.print();

	matrix D(A.row,A.col);
	for (int i=0;i<A.row;i++)	D.array[i][i] = 1/A.array[i][i];  // already compute D inv


	//yeah~! power iteration to compute largest eigenvalue
	matrix b(A.col,1);
	b.array[0][0]=1;
	matrix bdif(b);
	k = 0;
//	for (int i=0;i<100;i++){  
	while (bdif.norm() > 0.000001){
		bdif = b;
		b = AF*b;
		b = b*(1/b.norm());
//		b.transpose().print();
		k++;
		bdif = b - bdif;
	} 

	printf("\nPower iteration converged after %d iterations \n",k);
	matrix prod1(b.transpose()*A*b);
	matrix prod2(b.transpose()*b);

	double lambda =prod1.array[0][0]/prod2.array[0][0];
	double omega = 3/(4*lambda);

//	printf("\nlambda = %g \t omega = %g \n",lambda,omega);
			

	P = (D.identity() - D*AF*omega)*P;
	P.print();









} 


	
