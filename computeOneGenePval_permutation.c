#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
//#include <cblas.h>
//#include <clapack.h>
#include <zlib.h>
//#include "lapack_wrapper.h"
#include <sys/stat.h>
#include <getopt.h>
#include <gsl/gsl_cdf.h>

#include "mkl.h"

// Constants for I/O routines
#define DEFAULT_TFAM_NCOLS 6
#define DEFAULT_TPED_NCOLS 4
#define DEFAULT_NDIGITS 10
#define DEFAULT_DELIMS " \t\r\n"
#define SZ_LONG_BUF 1000000
#define DEFAULT_SIZE_MATRIX 1000000
#define DEFAULT_SIZE_HEADER 100000
#define MAX_NUM_MARKERS 20000000 // 20M max # snps
#define SZBUF 1024
#define DBL_MISSING -1e99 
#define DEFAULT_TPED_NUM_HEADER_COLS 4
#define DEFAULT_TFAM_NUM_HEADER_COLS 6
#define DEFAULT_TPED_SNPID_INDEX 1
#define DEFAULT_PHENO_NUM_HEADER_COLS 2
#define KINSHIP_IBS_MEAN 1
#define KINSHIP_IBS_RAND 2
#define KINSHIP_BALDING_NICHOLS 3
#define PROBE_NCOLS 3

// Constants for numerical routines
#define FPMIN 1e-30
#define MAXIT 1000
#define EIGEN_EPS 1e-10
#define EIGEN_ADD 1.
#define TOL 1e-6
#define EPS 1e-10
#define DEFAULT_NGRIDS 100
#define DEFAULT_LLIM -10
#define DEFAULT_ULIM 10
#define XACCU 1e-4
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define MAXRSID 100
#define EIGEN_FAIL -100000
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define TOO_FEW_SAMPLE 20

#define NUM_COV 2
#define INTERCEPT_INDEX 0
#define SNP_SLOPE_INDEX 1

/* 2x2 matrix - diagonals are 0, 3
0 2
1 3
*/

#define COV_X_SNP_POS 3

static int g_verbose = 1;
static int very_verbose = 0;
static char ct = 'T';
static char cn = 'N';
static char cv = 'V';
static char cl = 'L';
static double onef = 1.0;
static double zerof = 0.0;
static double minusonef = -1.0;
static int onen = 1;

void swap (double *a, double *b){
double temp = *a;
*a = *b;
*b = temp;
}

void randomize ( double arr[], int n ){
int i = 0;
for (i = n-1; i > 0; i--)
	{
		int j = rand() % (i+1);
		swap(&arr[i], &arr[j]);
	}
}

void vector_copy(int n, double* X, double* Y) {
/*
	Copies a vector X of lenght n to vector Y
*/
int i;
for (i=0; i<n; i++) Y[i]=X[i];
}

int matrix_invert(int n, double *X, double *Y) {
/*  
		Calculates the inverse of the n*n matrix X: Y = X^-1
		Does not change the value of X, unless Y=X
*/

int info=0;

/*  When X!=Y we want to keep X unchanged. Copy to Y and use this as working variable  */
if (X!=Y) vector_copy(n*n, X, Y);

/*  We need to store the pivot matrix obtained by the LU factorisation  */
int *ipiv;
ipiv=malloc(n*sizeof(int));
if (ipiv==NULL) {
	printf("malloc failed in matrix_invert\n"); 
	return 2;
}

/*  Turn Y into its LU form, store pivot matrix  */
//info = clapack_dgetrf (CblasColMajor, n, n, Y, n, ipiv);
dgetrf(&n,&n,Y,&n,ipiv,&info);

/*  Don't bother continuing when illegal argument (info<0) or singularity (info>0) occurs  */
if (info!=0) return info;

int lwork = n*16;
double* work = malloc(sizeof(double)*lwork);
/*  Feed this to the lapack inversion routine.  */
// info = clapack_dgetri (CblasColMajor, n, Y, n, ipiv);
dgetri(&n,Y,&n,ipiv,work,&lwork,&info);

/*  Cleanup and exit  */
free(ipiv);
free(work);
return info;
}

void transpose(double *src, double *dst, int N, int M) {
	int n = 0; 
	for(n = 0; n<N*M; n++) {
			int i = n/N;
			int j = n%N;
			dst[n] = src[M*j + i];
	}
}

void matrixMultiply ( double *A, double *B, int N, int K , int M, double *C) {
// A is N by M ( N = row_A ; M = col_A ) 
// B is M by K ( K = col_B) 
int i; 
int j; 
int l; 
for(i=0; i<N; i++) {
	for( j=0; j<K; j++) {
		double tmp = 0;
		for( l=0; l<M; l++) {
				tmp += A[M*i+l]*B[K*l+j];
		}
		C[K*i + j] = tmp;
	}
}
}

double betacf(double a, double b, double x) {
int m,m2; 
double aa,c,d,del,h,qab,qam,qap;

qab=a+b; 
qap=a+1.0; 
qam=a-1.0; 
c=1.0; 
d=1.0-qab*x/qap; if (fabs(d) < FPMIN) d=FPMIN; d=1.0/d;
h=d; 
for (m=1;m<=MAXIT;m++) {
	m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2)); 
	d=1.0+aa*d; 
	if (fabs(d) < FPMIN) d=FPMIN; 
	c=1.0+aa/c; 
	if (fabs(c) < FPMIN) c=FPMIN; 
	d=1.0/d; 
	h *= d*c; 
	aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
	d=1.0+aa*d; 
	if (fabs(d) < FPMIN) d=FPMIN; 
	c=1.0+aa/c; 
	if (fabs(c) < FPMIN) c=FPMIN; 
	d=1.0/d; 
	del=d*c; 
	h *= del; 
	if (fabs(del-1.0) < EPS) break;
} 
if (m > MAXIT) {
	fprintf(stderr,"a or b too big, or MAXIT too small in betacf %lf %lf %lf",a,b,x); 
	abort();
}
return h;
}

double gammln(double xx) { // estimates the "gamma" in the t-distribution CDF ... somehow...
double x,y,tmp,ser; 
static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5}; 
int j;
y=x=xx; 
tmp=x+5.5; 
tmp -= (x+0.5)*log(tmp); 
ser=1.000000000190015; 
for (j=0;j<=5;j++) 
	ser += cof[j]/++y; 
return -tmp+log(2.5066282746310005*ser/x);
}

double betai(double a, double b, double x) {
// not sure what's happening here, but it looks like it computes the part in t-CDF 
double gammln(double xx); 
void nrerror(char error_text[]); 
double bt;
if (x < 0.0 || x > 1.0) {
	fprintf(stderr,"Bad x in routine betai"); 
	abort();
}
if (x == 0.0 || x == 1.0) bt=0.0; 
else bt=exp((gammln(a+b)-gammln(a)-gammln(b))+(a*log(x))+(b*log(1.0-x)));
if (x < (a+1.0)/(a+b+2.0))
	return bt*betacf(a,b,x)/a; 
else
	return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double tcdf(double t, double nu) { // the CDF of the t-distribution 
if ( isnan(t) ) return 1.;
else return betai(nu/2.,0.5,nu/(nu+t*t));
}

void doLinearModel ( double *myarr, int numInd, int numCov, double *y , double *beta ,double *snp_stderr, double *pvalue) {

double tX [numCov][numInd]; // X is [col][row], so we now transpose it 
transpose(myarr, *tX, numInd, numCov); // transpose is done 

double tXX [numCov][numCov];  // t(x) %*% x 
matrixMultiply ( *tX, myarr, numCov, numCov , numInd, *tXX) ; 

double itXX [numCov][numCov]; // inverse 
matrix_invert(numCov, *tXX, *itXX); 

double itXXtX [numCov] [numInd] ; // (X^T X)^{-1} X^T 
matrixMultiply ( *itXX, *tX, numCov, numInd , numCov, *itXXtX) ; 

double itXXtXy [numCov] [1] ; // (X^T X)^{-1} X^T y = beta 
matrixMultiply ( *itXXtX, y, numCov, 1 , numInd, *itXXtXy) ; // since y is 1-D, not need the *y 

double yHat [numInd]; // this is estimated yHat
matrixMultiply ( myarr, *itXXtXy, numInd, 1 , numCov, yHat) ; 

double sigma_hat = 0;
int people = 0; 
for (people =0; people < numInd; people++ ){
	sigma_hat += pow( yHat[people] - y[people] , 2.0); 
}
sigma_hat = sigma_hat/(1.0*(numInd-numCov)); // this is sigma^2_e 

double stderr; 
stderr = sqrt ( itXX [(numCov-1)][(numCov-1)] * sigma_hat ); 
*snp_stderr = stderr; 

double betaAsIs = itXXtXy[numCov-1][0]; 
double stdBeta = betaAsIs/stderr; // standardized beta 
*beta = betaAsIs; 
double snp_ps = tcdf( -1.0*fabs(stdBeta) , (numInd-numCov));
*pvalue = snp_ps;	

}

void print_help(void) {
fprintf(stderr,"See instruction: github/datduong/FUNC-eGene.");
}

double normPdf_ALT (double x, double mean);
double normPdf_NULL (double x);
double computeGradient (double testStatistics, double mean, double wts ) ; 
float * read_wts ( char* file, float* arr, int stop ); 

double MEAN = 3.5 ; // to be replaced by input 

int main(int argc, char** argv) {

clock_t tic = clock(); 

float num_holder = 0; // temporary holder 
char str_holder [100]; // temporary holder

char *covf, *phenof, *genof, *lr_outFile, *obsValues_outFile, *wtsf, *pval, *bothpval; 
int c, numInd, numSNP, numCov, optimize, numPermu, seed; 

optimize = 1; // if 1, then do optimization (by default)
covf = phenof = genof = lr_outFile = obsValues_outFile = wtsf = pval = bothpval = NULL; 
numInd = numSNP = numCov = numPermu = seed = -1; 

static struct option long_options[] = {
	{"var", required_argument, 0, 'a' },
	{"expr", required_argument, 0, 'b' },
	{"geno", required_argument, 0, 'c' },
	{"lr", required_argument, 0, 'd' },
	{"detail", required_argument, 0, 'e' },
	{"numSNP", required_argument, 0, 'f' },
	{"numCov", required_argument, 0, 'g' },
	{"help", no_argument, 0, 'h' },
	{"w", required_argument, 0, 'i' }, 
	{"alt", required_argument, 0, 'j' },
	{"pval", required_argument, 0, 'k' },
	{"final", required_argument,0,'m'},
	{"reuse", required_argument,0,'n'},
	{"permutation ", required_argument,0,'o'},
	{"seed",required_argument,0,'p'},
	{0, 0, 0, 0 }
	};

int long_index = 0;
	while ((c = getopt_long_only(argc, argv, "a:b:c:d:e:f:g:h:i:j:k:m:n:o:p",long_options, &long_index)) != -1 ) {
	switch(c) {
	case 'a':
		covf = optarg;
	break;
	case 'b':
		phenof = optarg;
	break;
	case 'c':
		genof = optarg;
	break;
	case 'd':
		lr_outFile = optarg;
	break;
	case 'e':
		obsValues_outFile = optarg;
	break;
	case 'f':
		numSNP = atoi(optarg);
	break;
	case 'g':
		numCov = atoi(optarg);
	break;
	case 'i': 
		wtsf = optarg;
	break;
	case 'j':
		MEAN = atof(optarg);
	break;
	case 'k': 
		pval = optarg;
	break; 
	case 'm': 
		bothpval = optarg;
		break;
	case 'n':
		optimize = atoi(optarg);
		break;
	case 'o':
		numPermu = atoi(optarg);
		break;
	case 'p':
		seed = atoi(optarg);
		break;
	default:
			fprintf(stderr,"Error : %c",c);
			print_help();
			abort();
		}
}

int k = 0; 
int i = 0; 
double lr = 0; 
double maxLR = -INFINITY; 
double minPval = INFINITY; 
double maxBeta = -INFINITY;

//double signif_pos = 0; // position 
int snpps = 0;
int geneps = 0 ; 
int distance = 0; 

printf ("Read genotypes.\n"); 
FILE *genoFile = fopen(genof,"r"); // genotype 
fscanf(genoFile, "%s", str_holder); // header chr
fscanf(genoFile, "%d", &numInd); // header number people 

printf ("Read phenotypes.\n"); 
FILE *phenoFile = fopen(phenof,"r");  // gene expression 
fscanf(phenoFile, "%s", str_holder); // gene name 
fscanf(phenoFile, "%f", &num_holder); // chromosome 
fscanf(phenoFile, "%d", &geneps); // gene position 
double y [numInd] ; // gene expression 
for (k = 0; k <numInd; k++) {	
	fscanf(phenoFile, "%f", &num_holder);
	y[k] = num_holder; 
}
fclose(phenoFile);

printf ("Read covariates.\n"); 
FILE *covariateFile = fopen(covf,"r"); // other variables 
char line[SZ_LONG_BUF];
	fgets(line, SZ_LONG_BUF, covariateFile ); // header	
double myarr [numInd] [numCov+2]; // add 2 because: assign genoFile into last column, assign intercept to col#1 

for (k = 0; k<numInd; k++ ){
	myarr[k][0] =1; // intercept 
}
for (i = 1; i < (numCov+1); i++) { // variable 
	fscanf(covariateFile, "%s", str_holder); // read string name 
	for (k = 0; k < numInd; k++) { // person 
		fscanf(covariateFile, "%f", &num_holder);
		myarr[k][i] = num_holder; // filled in the person, byrow=TRUE 
	}
}
fclose(covariateFile);

// read in the priors here....
printf ("Read priors.\n"); 
float *snp_weights ; 
float array [numSNP]; // use &array to get the pointer address 
snp_weights = read_wts( wtsf, (float*)&array , numSNP); // casting a pointer . 
int w = 0; 
float sum_w = 0; 
for (w=0; w<numSNP; w++) { // compute sum of priors 
	sum_w += *(snp_weights+w); 
}
if (sum_w != 1) {
	for (w=0; w<numSNP; w++) { // normalized priors to sum up to 1 
		float temp = *(snp_weights+w)/sum_w;
		*(snp_weights+w) = temp; 
	}
}

// begin compute lr
printf ("Compute weighted_likelihood.\n"); 
FILE *ob = fopen(obsValues_outFile,"w+"); 
fprintf (ob, "pos\tdistance\tbeta\tsd_beta\tweighted_sd_beta\tweighted_lr\n");

// store SNP, reduce number read/write to files 
double snp_genotype [numInd][numSNP]; 
int min_index = 0 ; 
int min_snp_pos = 0; 
double best_beta = 0; 
double beststderr = 0; 
 
for ( i = 0 ; i < numSNP; i++){ // over each SNP 

	fscanf(genoFile, "%s", str_holder); // name 
	fscanf(genoFile, "%f", &num_holder); // chromosome, not too important 
	fscanf(genoFile, "%d", &snpps); // position 
	
	for (k=0; k<numInd; k++){
		fscanf(genoFile, "%f", &num_holder); // actual snps 
		myarr[k][numCov+1] = num_holder; 
		snp_genotype[k][i] = num_holder; 
		//printf ("%f ", num_holder); 
	}

	double beta = 0 ; // NOT standardized beta.  
	double snp_stderr = 0;
	double pvalue = 0; 	
	int numvar = numCov + 2; 
	doLinearModel ( *myarr, numInd, numvar, y, &beta, &snp_stderr, &pvalue ); // pass by ref. 
	
	double sdzBeta = beta/snp_stderr; // this beta ~ t-density 

	double tempps = 1.0-pvalue/2.0;
			double this_beta = 0;
			vdCdfNormInv(1, &tempps, &this_beta); // okay, the inv(norm-density) is 1-to-1. 

			if (optimize==1){
				// the underlying "source" of the llh is the normal-density 
				// sup power. so we need to set the mu as beta_hat
				lr = computeGradient ( this_beta, this_beta, *(snp_weights+i)) ; // STANDARDIZED BETA, SD=1
			}
			else if (optimize==0){
					lr = computeGradient ( this_beta, MEAN, *(snp_weights+i)) ; // STANDARDIZED BETA, SD=1
			}
			else {
					printf("\noptimize indication not recognized.\n");
			}	
			if (lr > maxLR){ // most significant
				maxLR = lr;  
	}
	if (pvalue < minPval ) {
		minPval = pvalue; 
		min_index = i; 
		min_snp_pos = snpps; 
	}
	if (fabs(sdzBeta) > maxBeta ) {
		maxBeta = fabs(sdzBeta); 
	}
	distance = geneps - snpps ; // tss - snp 
	fprintf (ob, "%d\t%d\t%f\t%f\t%f\t%f\n", snpps, distance, beta,snp_stderr, this_beta, lr);
}

printf ("snp with min p-value: pos %d, index: %d, min p-value: %f\n", min_snp_pos, min_index, minPval); 
fclose(genoFile); 
fclose(ob);

FILE *fout = fopen(lr_outFile,"w+"); // most significant SNP 
fprintf (fout, "%f\n", maxLR);
fclose(fout); 

FILE *fout2 = fopen(pval,"w+"); 
fprintf(fout2, "%f\n", minPval);
fclose(fout2); 

/* set seed so we can do many parallel jobs */
/*VSLStreamStatePtr stream;
vslNewStream(&stream, VSL_BRNG_MT19937, time(NULL)+seed); 
*/
srand(time(NULL)+seed);

if (numPermu > 0){

	printf("begin permutations.\n");

	int compareCounter = 0; 
	int compareCounterGrd = 0; 
	int compareCounter_beta = 0; 

	int numP = 0;
	for ( numP = 0; numP < numPermu; numP++) {
		
		if (numP >= 100) {
				if (numP % (numPermu/10) == 0) {
				printf("%d # permutations done\n",numP);
					}
			}

		randomize(y,numInd); 
		double bestBetaSim = -INFINITY;
		double bestGradSim = -INFINITY; 
		double bestPvalSim = 1; 

		for ( i = 0 ; i < numSNP; i++){ // over each SNP 

			for (k=0; k<numInd; k++){
				double tempv = snp_genotype[k][i];
				myarr[k][numCov+1] = tempv; 
			}

			double beta = 0 ; // NOT standardized beta.  
			double snp_stderr = 0;
			double pvalue = 0; 	
			int numvar = numCov + 2; 
			doLinearModel ( *myarr, numInd, numvar, y, &beta, &snp_stderr, &pvalue ); // pass by ref. 
		
			double sdzBeta = beta/snp_stderr; 
	
			double tempps = 1.0-pvalue/2.0;
					double this_beta = 0;
					vdCdfNormInv(1, &tempps, &this_beta); // okay, the inv(norm-density) is 1-to-1. 

			if (optimize==1){
				// sup power. so we need to set the mu as beta_hat 
				lr = computeGradient ( this_beta, this_beta, *(snp_weights+i)) ; // STANDARDIZED BETA, SD=1 
			}
			else if (optimize==0){
				lr = computeGradient ( this_beta, MEAN, *(snp_weights+i)) ; // STANDARDIZED BETA, SD=1 
			}
			else {
				printf("\noptimize indication not recognized.\n");
			}
			if (lr > bestGradSim){ // most significant  
				bestGradSim = lr;  
			}
			if (pvalue < bestPvalSim) {
				bestPvalSim = pvalue; 
			}
		}

		if (bestPvalSim < minPval) {
			compareCounter = compareCounter + 1; 
		}
		if (bestGradSim > maxLR) {
			compareCounterGrd = compareCounterGrd + 1; 
		}

	}	

	double p1 = (double)compareCounter/numPermu; 
	double p2 = (double)compareCounterGrd/numPermu;
	printf ("pvalue based on min-p %f\n", p1 ); 
	printf("pvalue based on lr %f\n", p2 );

	FILE *bothp = fopen(bothpval,"w+");
	fprintf(bothp, "uniform_pvalue %f nonuniform_pvalue %f\n", p1, p2);
	fclose(bothp); 
} else {
	printf ("eGene p-value not computed because numPermu is less than 0.\n");
}

clock_t toc = clock();

printf("\nElapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	
return 0;

}


float * read_wts ( char* file, float* arr, int stop ){
FILE *fp; 
if((fp=fopen(file,"r"))==NULL){
		printf("priors not found...\n");
		exit(1);
}
int i = 0; 
for (i=0; i<stop; i++){
	fscanf(fp,"%f", &arr[i]); // the value will be store to array. 
}
return arr; // return the pointer to the array 
}

double computeGradient (double testStatistics, double mean, double wts ) {
	double part1 = 0; 
	part1 = log( 0.5 * ( exp (testStatistics*mean) + exp(-1.0*testStatistics*mean) ) ) - 0.5*mean*mean;
	return ( log(wts) + part1 ); 
}
