#include "Benchmarks.h"

Benchmarks::Benchmarks(RunParameter runParam){
	cout<<"Benchmarks Class initialization"<<endl;
	dimension = runParam.dimension;		
	nonSeparableGroupSize = runParam.nonSeparableGroupSize;
	MASK = ((L(1)) << (L(48))) - (L(1));
	m_havenextGaussian = 0;

	// allocate the memory
	anotherz = new double[dimension];
	anotherz1= new double[nonSeparableGroupSize];
	anotherz2= new double[dimension - nonSeparableGroupSize];

	// Runtime Parameters setting
	setOvectorToZero = false;

	functionInitRandomSeed = L(runParam.initRandomSeed);
	m_seed= functionInitRandomSeed;
	M  = 0x5DEECE66D;
	A  = 0xB;
}

Benchmarks::~Benchmarks(){
	delete[] anotherz;
	delete[] anotherz1;
	delete[] anotherz2;
	cout<<"Benchmarks Class Destroyed"<<endl;
}

int Benchmarks::next(int bits) {
  int64_t s;
  int64_t result;
  m_seed = s = (((m_seed * M) + A) & MASK);
  result = (s >> (L(48 - bits)));
  return((int)result);
}

int Benchmarks::nextInt(int n) {
	int bits, val;

	if ((n & (-n)) == n) {
		return((int) ((n * L(next(31))) >> L(31)));
	}

	do {
		bits = next(31);
		val  = bits % n;
	} while (bits - val + (n - 1) < 0);

	return(val);
}

double Benchmarks::nextDouble(){
  return ((((L(next(26))) <<
            (L(27))) + (L(next(27)))) / (double) ((L(1)) << (L(53))));
}

double Benchmarks::nextGaussian(){
  double multiplier, v1, v2, s;

  if (m_havenextGaussian) {
    m_havenextGaussian = 0;
    return(m_nextGaussian) ;
  }

  do {
    v1 = ((D(2.0) * nextDouble()) - D(1.0));
    v2 = ((D(2.0) * nextDouble()) - D(1.0));
    s  = ((v1 * v1) + (v2 * v2));
  } while ((s >= D(1.0)) || (s <= D(0.0)));

  multiplier = sqrt(D(-2.0) * log(s) / s);
  m_nextGaussian    = (v2 * multiplier);
  m_havenextGaussian = 1;
  return (v1 * multiplier);
}

double* Benchmarks::createShiftVector(int dim, double min,double max) {
  double* d;
  double  hw, middle;
  double  s;
  int     i;
  hw     = (0.5 * (max - min));
  middle = (min + hw);
  d      = (double*)malloc(sizeof(double) * dim);



  for (i = (dim - 1); i >= 0; i--) {
	  if (setOvectorToZero == true){
		  d[i] = s;
	  }else{
		  do {
			  s = (middle + (nextGaussian() * hw));
		  } while ((s < min) || (s > max));
		  d[i] = s;
	  }
  }

  return(d);
}

int* Benchmarks::createPermVector(int dim){
  int* d;
  int  i, j, k, t;
  d = (int*)malloc(sizeof(int) * dim);

  for (i = (dim - 1); i >= 0; i--) {
    d[i] = i;
  }

  for (i = (dim << 3); i >= 0; i--) {
    j = nextInt(dim);

    do {
      k = nextInt(dim);
    } while (k == j);

    t    = d[j];
    d[j] = d[k];
    d[k] = t;
  }

  return(d);
}

//Create a random rotation matrix
double** Benchmarks::createRotMatrix(int dim){
  double** m;
  int      i, j, k;
  double   dp, t;
  m = (double**)malloc(sizeof(double*) * dim);

  for(i = 0; i < dim; i++) {
    m[i] = (double*)malloc(sizeof(double) * dim);
  }

loop:

  for (;;) {
    for (i = (dim - 1); i >= 0; i--) {
      for (j = (dim - 1); j >= 0; j--) {
        m[i][j] = nextGaussian();
      }
    }

    // main loop of gram/schmidt
    for (i = (dim - 1); i >= 0; i--) {
      for (j = (dim - 1); j > i; j--) {
        // dot product
        dp = 0;

        for (k = (dim - 1); k >= 0; k--) {
          dp += (m[i][k] * m[j][k]);
        }

        // subtract
        for (k = (dim - 1); k >= 0; k--) {
          m[i][k] -= (dp * m[j][k]);
        }
      }

      // normalize
      dp = 0;

      for (k = (dim - 1); k >= 0; k--) {
        t   = m[i][k];
        dp += (t * t);
      }

      // linear dependency -> restart
      if (dp <= 0) {
        goto loop;
      }

      dp = (1 / sqrt(dp));

      for (k = (dim - 1); k >= 0; k--) {
        m[i][k] *= dp;
      }
    }

    return(m) ;
  }
}

/**
 * Create a random rotation matrix
 */
double* Benchmarks::createRotMatrix1D(int dim){
  double** a;
  double*  b;
  int      i, j, k;
  a = createRotMatrix(dim);
  b = (double*)malloc(sizeof(double) * (dim * dim));
  k = 0;

  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      b[k++] = a[i][j];
    }
  }

  return(b);
}


void Benchmarks::lookupprepare() {
  double pownum;
  int    i;
  i         = (dimension - 1);
  pownum    = (1.0 / i);
  lookup    = (double*)malloc(dimension * sizeof(double));
  lookup[i] = 1.0e6;
  lookup[0] = 1.0;

  for (--i; i > 0; i--) {
    lookup[i] = pow(1.0e6, i * pownum);
  }
}

/* 
 * Basic Mathematical Functions' Implementation
 */
double Benchmarks::elliptic(double*x,int dim) {
  double result = 0.0;
  int    i;

  for(i = dim - 1; i >= 0; i--) {
    result += lookup[i] * x[i] * x[i];
    //printf("%d %lf %lf\n",i,lookup[i],x[i]);
  }

  return(result);
}


double Benchmarks::rastrigin(double*x,int dim){
  double sum = 0;
  int    i;

  for(i = dim - 1; i >= 0; i--) {
    sum += x[i] * x[i] - 10.0 * cos(2 * PI * x[i]) + 10.0;
  }

  return(sum);
}

double Benchmarks::ackley(double*x,int dim){
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum;
  int    i;

  for(i = dim - 1; i >= 0; i--) {
    sum1 += (x[i] * x[i]);
    sum2 += cos(2.0 * PI * x[i]);
  }

  sum = -20.0 * exp(-0.2 * sqrt(sum1 / dim)) - exp(sum2 / dim) + 20.0 + E;
  return(sum);
}

double* Benchmarks::multiply(double*vector, double*matrix,int dim){
  int    i,j;
  double*result = (double*)malloc(sizeof(double) * dim);

  for(i = dim - 1; i >= 0; i--) {
    result[i] = 0;

    for(j = dim - 1; j >= 0; j--) {
      result[i] += vector[j] * matrix[dim * j + i];
    }
  }

  return(result);
}

//Create a random permutation vector of the numbers 0 to dim-1
double Benchmarks::rot_elliptic(double*x,int dim){
  double result = 0.0;
  int    j;
  double*z = multiply(x,RotMatrix,dim);

  for(j = dim - 1; j >= 0; j--) {
    result = elliptic(z,dim);
  }

  free(z);
  return(result);
}
