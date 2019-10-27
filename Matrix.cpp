#include "Matrix.h"
#include "Vector.h"

//using namespace std;

mat::mat() {
	rows = 0;
	cols = 0;
}

mat::mat(uint8_t m, uint8_t n) {
	rows = m;
	cols = n;
}

// mat Addition
mat mat::operator+(const mat& m) {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j] + m.data[i][j];
		}
	}
	return b;
}

void mat::operator+=(const mat& m) {
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] += m.data[i][j];
		}
	}
}

// mat Subtraction
mat mat::operator-(const mat& m) {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j] - m.data[i][j];
		}
	}
	return b;
}

void mat::operator-=(const mat& m) {
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] -= m.data[i][j];
		}
	}
}

// Element-by-element multiplication
mat mat::operator^(const mat& m) {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j] * m.data[i][j];
		}
	}
	return b;
}

void mat::operator^=(const mat& m) {
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] *= m.data[i][j];
		}
	}
}

// Element-by-element division
mat mat::operator/(const mat& m) {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j] / m.data[i][j];
		}
	}
	return b;
}

void mat::operator/=(const mat& m) {
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] /= m.data[i][j];
		}
	}
}

// mat multiplication
mat mat::operator*(const mat& m) {
	mat b(rows, m.cols);
	for(uint8_t i=0; i<b.rows; ++i) {
		for(uint8_t j=0; j<b.cols; ++j) {
			double val = 0;
			for(uint8_t k=0; k<cols; k++) {
				val += data[i][k]*m.data[k][j];
			}
			b.data[i][j] = val;
		}
	}
	return b;
}

void mat::operator*=(const mat& m) {
	mat b(rows, m.cols);
	for(uint8_t i=0; i<b.rows; ++i) {
		for(uint8_t j=0; j<b.cols; ++j) {
			double val = 0.0f;
			for(uint8_t k=0; k<cols; ++k) {
				val += data[i][k]*m.data[k][j];
			}
			b.data[i][j] = val;
		}
	}
	copy(b);
}

// scalar multiplication
mat mat::operator*(const double k) {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j] * k;
		}
	}
	return b;
}

void mat::operator*=(const double k) {
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] *= k;
		}
	}
}

// scalar division
mat mat::operator/(const double k) {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j] / k;
		}
	}
	return b;
}

void mat::operator/=(const double k) {
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] /= k;
		}
	}
}

// Matrix Initializers
mat mat::zeros(uint8_t m, uint8_t n) {
	mat b(m,n);
	for(uint8_t i=0; i<b.rows; ++i) {
		for(uint8_t j=0; j<b.cols; ++j) {
			b.data[i][j] = 0;
		}
	}
	return b;
}

mat mat::ones(uint8_t m, uint8_t n) {
	mat b(m,n);
	for(uint8_t i=0; i<b.rows; ++i) {
		for(uint8_t j=0; j<b.cols; ++j) {
			b.data[i][j] = 1;
		}
	}
	return b;
}

mat mat::identity(uint8_t n) {
	mat b(n,n);
	for(uint8_t i=0; i<b.rows; ++i) {
		for(uint8_t j=0; j<b.cols; ++j) {
			if (i == j) {
				b.data[i][j] = 1;
			} else {
				b.data[i][j] = 0;
			}
		}
	}
	return b;
}

// Copy matrix
void mat::copy(const mat& m) {
	rows = m.rows;
	cols = m.cols;
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			data[i][j] = m.data[i][j];
		}
	}
}

void mat::resize(uint8_t m, uint8_t n) {
	assert(m<=MAX_LEN_MAT && n<=MAX_LEN_MAT);
    rows = m;
    cols = n;
}

// Compute determinant - FIXME: Needs work!
double mat::det() {
    if(rows == 1) {
        return data[0][0];
    } else {
        double d = 0;
        for(uint8_t j=0; j<cols; ++j) {
            if (j%2 == 0) {
                d += cofactor(0,j).det()*data[0][j];
            } else {
                d -= cofactor(0,j).det()*data[0][j];
            }
        }
        return d;
    }
}

// Compute transpose
mat mat::t() {
	mat b(cols, rows);
	for(uint8_t i=0; i<cols; ++i) {
		for(uint8_t j=0; j<rows; ++j) {
			b.data[i][j] = data[j][i];
		}
	}
	return b;
}

// Cofactor matrix
mat mat::cofactor(uint8_t m, uint8_t n) {
    mat b(rows-1, cols-1);

    uint8_t ii=0, jj=0;
	uint8_t i, j;
    // Top-left square
	for(i=0; i<m; ++i) {
		for(j=0; j<n; ++j) {
            b.data[ii][jj] = data[i][j];
            ++jj;
		}
		++ii;
		jj=0;
	}

    ii=m; jj=n;
    // Bottom-right square
    for(i=m+1; i<rows; ++i) {
		for(j=n+1; j<cols; ++j) {
            b.data[ii][jj] = data[i][j];
            ++jj;
		}
		++ii;
		jj=n;
	}

    ii=0; jj=n;
    // Top-right rectangle
    for(i=0; i<m; ++i) {
        for(j=n+1; j<cols; ++j) {
            b.data[ii][jj] = data[i][j];
            ++jj;
		}
		++ii;
		jj=n;
    }

    ii=m; jj=0;
    // Bottom-left rectangle
    for(i=m+1; i<rows; ++i) {
		for(j=0; j<n; ++j) {
            b.data[ii][jj] = data[i][j];
            ++jj;
		}
		++ii;
		jj=0;
	}

    return b;
}


// Compute inverse
mat mat::inv() {
	// Gauss-Jordan Elimination Method - based on my MATLAB implementation
	// NOTE - A and I and made to be separate matrices (not concatenated) to save
	// memory space (so, a lower MAX_LEN_MAT can be used)
	mat A = get_mat();
	mat I = identity(rows);
	vec tmp(A.get_cols());

	int i,j;
	double k;

	// STEP 1 - Makes lower triangular zeros
	for(i=1; i<(int)rows; ++i) {
		for(j=0; j<i; ++j) {
			k = A.data[i][j]/A.data[j][j];
			
			A.get_row(j, tmp);
			tmp*=k;
			A.row_substract(i, tmp);
			
			I.get_row(j, tmp);	
			tmp*=k;
			I.row_substract(i, tmp);
		}
	}
	// STEP 2 - Makes upper triangular zeros
	for(i=(int)rows-2; i>=0; --i) {
		for(j=(int)rows-1; j>=i+1; --j) {
			k = A.data[i][j]/A.data[j][j];
			
			A.get_row(j, tmp);
			tmp*=k;
			A.row_substract(i, tmp);

			I.get_row(j, tmp);
			tmp*=k;
			I.row_substract(i, tmp);
		}
	}

	// STEP 3 - Make the left side an identity matrix
	for(i=0; i<(int)rows; ++i) {
		I.row_divide(i, A.data[i][i]);
	}

	return I; // The identity matrix transforms into the inverse
}


// Matrix concatenation
mat mat::cols_cat(const mat& m) {
	mat b(rows, cols + m.cols);
	uint8_t j;
	for(uint8_t i=0; i<rows; ++i) {
		for(j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j];
		}
		for(j=cols; j<cols+m.cols; ++j) {
			b.data[i][j] = m.data[i][j-cols];
		}
	}
	return b;
}

mat mat::rows_cat(const mat& m) {
	mat b(rows + m.rows, cols);
	uint8_t i;
	for(uint8_t j=0; i<cols; ++j) {
		for(i=0; i<rows; ++i) {
			b.data[i][j] = data[i][j];
		}
		for(i=rows; i<rows+m.rows; ++i) {
			b.data[i][j] = m.data[i-rows][j];
		}
	}
	return b;
}

// Matrix access functions


mat mat::get_subm(uint8_t m1, uint8_t m2, uint8_t n1, uint8_t n2) {
	mat b(m2-m1+1, n2-n1+1);
	for (uint8_t i=0; i<=m2-m1; ++i) {
		for (uint8_t j=0; j<=n2-n1; ++j) {
			b.data[i][j] = data[m1+i][n1+j];
		}
	}
	return b;
}

void mat::get_subm(uint8_t m1, uint8_t m2, uint8_t n1, uint8_t n2, mat& b) {
	b.resize(m2-m1+1, n2-n1+1);
	for (uint8_t i=0; i<=m2-m1; ++i) {
		for (uint8_t j=0; j<=n2-n1; ++j) {
			b.data[i][j] = data[m1+i][n1+j];
		}
	}
}

vec mat::get_row(uint8_t m) {
	vec b(cols);
	for (uint8_t j=0; j<cols; ++j) {
		b.data[j] = data[m][j];
	}
	return b;	
}

void mat::get_row(uint8_t m, vec& b) {
	assert(b.dim == cols);
	
	for (uint8_t j=0; j<cols; ++j) {
		b.data[j] = data[m][j];
	}
}

vec mat::get_col(uint8_t n) {
	vec b(rows);
	for (uint8_t i=0; i<rows; ++i) {
		b.data[i] = data[i][n];
	}
	return b;
	
}

void mat::set_subm(uint8_t m, uint8_t n, const mat& subm) {
	for (uint8_t i=0; i<subm.rows; ++i) {
		for (uint8_t j=0; j<subm.cols; ++j) {
			data[m+i][n+j] = subm.data[i][j];
		}
	}
}

void mat::set_row(uint8_t m, const vec& row) {
	assert(row.dim == cols);
	
	for (uint8_t j=0; j<cols; ++j) {
		data[m][j] = row.data[j];
	}
}

void mat::row_substract(uint8_t m, const vec& row) {
	assert(row.dim == cols);
	
	for (uint8_t j=0; j<cols; ++j) {
		data[m][j] -= row.data[j];
	}
}

void mat::row_divide(uint8_t m, const double v ) {
	for (uint8_t j=0; j<cols; ++j) {
		data[m][j] /= v;
	}
}

void mat::set_col(uint8_t n, const vec& col) {
	assert(col.dim == rows);
	
	for (uint8_t i=0; i<rows; ++i) {
		data[i][n] = col.data[i];
	}
}

mat mat::get_mat() {
	mat b(rows, cols);
	for(uint8_t i=0; i<rows; ++i) {
		for(uint8_t j=0; j<cols; ++j) {
			b.data[i][j] = data[i][j];
		}
	}
	return b;
}

double& mat::operator()(uint8_t m, uint8_t n) {
	assert(m<rows && n<cols);
    return data[m][n];
}

void mat::print(String & str) {
  str="Rows:";
  str+=rows;
  str+="\nCols:";
  str+= cols;
  str += "\n";
  for (uint8_t i=0;i<cols;++i) {
    if (i==0) {
      str+="[";
    } else {
      str+= "|";
    }
    for (uint8_t j=0;j<rows;++j) {
        if (data[i][j] >=0.0f){
          str+=" ";
        }
        str+=String(data[i][j],6);
        str+=" ";
    }
    if (i==rows-1) {
      str+="]";
    } else {
      str+= "|";
    }
    str+="\n";
  }
}
