#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class matrix {
	
	public:
	int row;
	int col;
	float array[100][100];
	
		matrix(int,int); //constuctor
		matrix(const matrix&); //constructor to copy a matrix
		matrix(matrix*); //constructor to copy a matrix
		void setElement(int r, int c, float e);
		float getElement (int r, int c);		
		matrix operator+(const matrix&);
		matrix operator-(const matrix&);
		matrix operator*(const matrix&);
		matrix operator*(const double);
		matrix inv(void);
		matrix identity(void);
		matrix transpose(void);
		void print(void);
		double norm(void);

};

matrix::matrix(int r, int c){ //construct a zero matrix
	row = r;
	col = c;
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			array[i][j] = 0;
		}
	}
};

matrix::matrix(const matrix& m){ //defines the copying constructor
	row = m.row;
	col = m.col;
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			array[i][j] = m.array[i][j];
		}
	}
}

matrix::matrix(matrix* m){ //defines the copying constructor
	row = m->row;
	col = m->col;
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			array[i][j] = m->array[i][j];
		}
	}
}



void matrix::setElement(int r, int c, float e){
	array[r][c] = e;
};

matrix matrix::operator+(const matrix& m){
	// first, make sure matrices can be added. if not, return original matrix
	if (row != m.row || col != m.col){
		printf("Matrix sizes do not match. Mission impossible.\n");
		return (*this);
	}
	matrix new_mat(row,col);
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			new_mat.array[i][j] = array[i][j] + m.array[i][j];
		}
	}
	return new_mat;
}

matrix matrix::transpose(void){
	matrix At(col,row);
	for (int i=0;i<row;i++){
		for (int j=0;j<col;j++){
			At.array[j][i] = array[i][j];
		}
	}
	return At;
}



matrix matrix::operator*(const double k){
	matrix new_mat(row,col);
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			new_mat.array[i][j] = array[i][j]*k;
		}
	}
	return new_mat;
}


matrix matrix::operator-(const matrix& m){
	// first, make sure matrices can be substracted. if not, return original matrix
	if (row != m.row || col != m.col){
		printf("Matrix sizes do not match. Mission impossible.\n");
		return (*this);
	}
	matrix new_mat(row,col);
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			new_mat.array[i][j] = array[i][j] - m.array[i][j];
		}
	}
	return new_mat;
}




float matrix::getElement(int r, int c){
	return array[r][c];
}


matrix matrix::operator*(const matrix& m){
	// first, make sure matrices can be multyplied. if not, return original matrix
	if (col != m.row){
		printf("Matrix sizes do not match. Mission impossible.\n");
		return (*this);
	}
	matrix new_mat(row,m.col);
	for (int i=0; i<row; i++){
		for (int j=0; j<m.col; j++){
			new_mat.array[i][j] = 0;
			for (int k=0;k<col;k++){
				new_mat.array[i][j] = new_mat.array[i][j] + array[i][k]*m.array[k][j];
			}
		}
	}
	return new_mat;
}



double matrix::norm(void){ // only for vector i.e. column or row matrices
	double S=1;
	if(col == 1){
		S = 0;
		for (int i=0 ; i<row ; i++ ){
			S = S + array[i][0];
		}
		S = sqrt(S);
	}
	if (row == 1){
		S = 0;
		for (int i=0 ; i<row ; i++ ){
			S = S + array[0][i];
		}
		S = sqrt(S);
	}
	return S;
}




void matrix::print(void){
	printf("\n");
	for (int i=0;i<row;i++){
		for (int j=0;j<col;j++){
			printf("%.10g \t ",array[i][j]);
		}
	printf("\n");
	}				
	printf("\n");
}




matrix matrix::identity(void){
	// first, make sure matrix is square. if not, return original matrix
	if (col != row){
		printf("Matrix sizes do not match. Mission impossible.\n");
		return (*this);
	}
	matrix aid(row,col);
	for (int i=0;i<row;i++)
		aid.array[i][i]=1;

	return aid;
}

matrix matrix::inv(void){
	// first, make sure matrix is square. if not, return original matrix
	if (col != row){
		printf("Matrix sizes do not match. Mission impossible.\n");
		return (*this);
	}
	matrix ainv(row,col);
	matrix a(this);
	ainv = ainv.identity();


	for (int k = 0; k < row - 1; k++) {		
		for (int i = k+1; i < row;  i++) {
			double factor = a.array[i][k]/a.array[k][k]; 
			for (int j = 0; j < row ; j++) {
				a.array[i][j] = a.array[i][j] - factor * a.array[k][j];
				//if (j < k+1)    a.array[i][j] = 0;
				ainv.array[i][j] = ainv.array[i][j] - factor * ainv.array[k][j];
			}
		}
	}

	for (int k = row-1; k > 0; k--) {		
		for (int i = k-1; i >= 0;  i--) {
			double factor = a.array[i][k]/a.array[k][k]; 
			for (int j = row-1; j >= 0 ; j--) {
				a.array[i][j] = a.array[i][j] - factor * a.array[k][j];
				//if (j > k-1)    a.array[i][j] = 0;
				ainv.array[i][j] = ainv.array[i][j] - factor * ainv.array[k][j];
			}
		}
	} 

	for (int i=0;i<row;i++){
		for (int j=0;j<col;j++){
			ainv.array[i][j] = ainv.array[i][j]/a.array[i][i];
		}
	}

	return ainv;
}



