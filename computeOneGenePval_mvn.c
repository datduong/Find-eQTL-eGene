#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <zlib.h>
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


#define NUM_COV 4
#define INTERCEPT_INDEX 0
#define ENV_SLOPE_INDEX 1
#define SNP_SLOPE_INDEX 2
#define GXE_SLOPE_INDEX 3

/* 4x4 matrix - diagonals are 0, 5, 10, 15
0  4  8  12
1  5  9  13
2  6  10 14
3  7  11 15
*/

#define COV_X_SNP_POS 10
#define COV_X_GXE_POS 15

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


double NormalCDFInverse(double p);
double RationalApproximation(double t);

// Input routines
void print_help(void);
void save_matrix(double* kins, int nr, int nc, int ndigits, char* filename);

int matrix_invert(int n, double *X, double *Y);

// stat
double tcdf(double t, double nu);

double getMinPvalue(int n, double* pvalue);

int readMAFFile(char* file, int* mafInt, int* snpPos, char* delims);

int readNumThreshold(char* fprdir);

void readFPRRatio(char* fprdir, int numMAF, int numThreshold, double** mafFPR, double** mafRatio, int* endThreshold, char* delims);

int minDiff(double value, double* A, int* candidates, int numCandidate);

int binarySearch(double value, double* A, int imax);

double correctPvalue(double pvalue, double* fpr, double* ratio, int endThreshold);

void sampleMVN(double maxBetaObserved, double minPvalueObserved, double maxGrdnObserved, int seed, double** corrWindow, int numSample, int* mafInt, double** mafFPR, double** mafRatio, int numThreshold, int* endThreshold, char* outf, int noPvalCor, int numWindow, int* numSNPInWindow, int doExPvalCor, double corFactor, float *snp_weights, int optimize);

int readCorrelation(char* corrf, int numSNP, double* correlation, char* delims);

double readPvalFile(char* file);

void splitCorrelationMatrix(double* correlation, double** corrWindow, int* snpPos, int numWindow, int numSNP, int* numSNPInWindowArray);

float * read_wts ( char* file, float* arr, int stop);

double normPdf_ALT (double x, double mean);
double normPdf_NULL (double x);
double computeGradient (double testStatistics, double mean, double wts ) ; 

double extractPermutatePval(char* file) ; 

int numSNP = 1;
int numInd = 0;
int windowSize = 1000000000;

double MEAN = 3.5 ; // this will be replaced by the input 

int main(int argc, char** argv) {
	
  clock_t tic = clock(); 
	
  char *corrf, *outf, *delims, *genof, *fprdir, *pvalf, *prior_file , *grdf;
  int c, numSample, seed, noPvalCor, doExPvalCor, optimize;
  double corFactor;

  corrf = outf = genof = fprdir = pvalf = prior_file = grdf = NULL; 
  numSNP = numSample = seed = optimize = -1; //numInd = 
  noPvalCor = doExPvalCor = 0;
  corFactor = -1;

  delims = DEFAULT_DELIMS;

  static struct option long_options[] = {
    {"corf", required_argument, 0, 'a' },
    {"output", required_argument, 0, 'b' },
    {"numSample", required_argument, 0, 'c' },
    {"numSNP", required_argument, 0, 'd' },
    {"seed", required_argument, 0, 'e' },
    {"geno", required_argument, 0, 'f' },
    {"fprdir", required_argument, 0, 'g' },
    {"help", no_argument, 0, 'h' },
    {"pval", required_argument, 0, 'i' },
    {"nocor", no_argument, 0, 'j' },
    {"windowSize", required_argument, 0, 'k' },
    {"extremePCor", no_argument, 0, 'l' },
    {"corFactor", required_argument, 0, 'm' },
		{"w", required_argument, 0, 'n' }, 
    {"alt", required_argument, 0, 'o'},
    {"lr", required_argument, 0, 'p'},
    {"reuse",required_argument, 0,'r'},
    //{"detail",required_argument, 0,'s'},
    {0, 0, 0, 0 }
  };

  int long_index = 0;
  while ((c = getopt_long_only(argc, argv, "a:b:c:d:e:f:g:hi:jk:lm:n:o:q:r:s:t",long_options, &long_index)) != -1 ) {
    switch(c) {
    case 'a':
      corrf = optarg;
      break;
    case 'b':
      outf = optarg;
      break;
    case 'c':
      numSample = atoi(optarg);
      break;
    case 'd':
      numSNP = atoi(optarg);
      break;
    case 'e':
      seed = atoi(optarg);
      break;
    case 'f':
      genof = optarg;
      break;
    case 'g':
      fprdir = optarg;
      break;
    case 'h':
      print_help();
      abort();
    case 'i':
      pvalf = optarg;
      break;
    case 'j':
      noPvalCor = 1;
      break;
    case 'k':
      windowSize = atoi(optarg);
      break;
    case 'l':
      doExPvalCor = 1; // change to 1 for "yes" do correct, 0 for "no"
      break;
    case 'm':
      corFactor = atof(optarg);
      break;
		case 'n':
      prior_file = optarg;
      break;
    case 'o':
      MEAN = atof (optarg);
      break;
    case 'p':
      grdf = optarg;
      break;
		case 'r':
      optimize = atoi(optarg);
      break;
    /*case 's':
      detail_out = optarg;
      break;*/
    default:
      fprintf(stderr,"Error : Unknown option character %c",c);
      print_help();
      abort();
    }
  }
  if (corrf == NULL) {
    printf("genotype correlation file must be specified\n");
    print_help();
    exit(-1);
  }
  if (outf == NULL) {
    printf("output prefix must be specified\n");
    print_help();
    exit(-1);
  }
  if (genof == NULL) {
    printf("genotype file must be specified\n");
    print_help();
    exit(-1);
  }
  if (fprdir == NULL) {
    printf("fpr directory must be specified\n");
    print_help();
    exit(-1);
  }
  if (pvalf == NULL) {
    printf("pval must be specified\n");
    print_help();
    exit(-1);
  }
  if (corFactor < 0) {
    printf("corFactor must be specified and > 0\n");
    print_help();
    exit(-1);
  }
  /*
  if (doExPvalCor == 1 && corFactor < 0) {
    printf("corFactor must be specified and > 0 with extremePCor option\n");
    print_help();
    exit(-1);
  }
  */
  if (doExPvalCor == 1 && numSample < 10000000) {
    printf("numSample >= 10M with extremePCor option\n");
    print_help();
    exit(-1);
  }
  if (numSNP < 0) {
    printf("numSNP must be specified and > 0\n");
    print_help();
    exit(-1);
  }
  if (numSample < 0) {
    printf("numSample must be specified and > 0\n");
    print_help();
    exit(-1);
  }
  if (seed < 0) {
    printf("seed must be specified and > 0\n");
    print_help();
    exit(-1);
  }
  noPvalCor = 0; 
  if (noPvalCor == 1) {
    printf("p-value correction is OFF\n");
  }
  else {
    printf("p-value correction is ON\n");
  }
  if (doExPvalCor == 1) {
    printf("Extreme p-value estimation is ON\n");
  }
  else {
    printf("Extreme p-value estimation is OFF\n");
  }

 
  MKL_Set_Num_Threads(1);

  // read priors here....
  float *snp_weights ;
  float array [numSNP]; // use &array to get the pointer address 
  snp_weights = read_wts( prior_file, (float*)&array , numSNP); // casting a pointer . 
 
  float sum_w = 0; 
  int w = 0;
	for (w=0; w<numSNP; w++) {
		sum_w += *(snp_weights+w); 
	}
	if (sum_w != 1) {
		for (w=0; w<numSNP; w++) {
			float temp = *(snp_weights+w)/sum_w;
			*(snp_weights+w) =  temp; 
		}
	}

  int* mafInt = (int*)malloc(sizeof(int)*numSNP);
  int* snpPos = (int*)malloc(sizeof(int)*numSNP);
  int numWindow = readMAFFile(genof, mafInt, snpPos, delims);

  double* correlation = (double*)malloc(sizeof(double)*numSNP*numSNP);
  int success = readCorrelation(corrf, numSNP, correlation, delims);

  double** corrWindow = (double**)malloc(sizeof(double*)*numWindow);
  int* numSNPInWindow = (int*)malloc(sizeof(int)*numWindow);
  printf("# of windows = %d, window size = %d\n",numWindow, windowSize);
  splitCorrelationMatrix(correlation, corrWindow, snpPos, numWindow, numSNP, numSNPInWindow);

  int numMAF = 50;
  int numThreshold = readNumThreshold(fprdir);

  double** mafFPR = NULL;
  double ** mafRatio = NULL;
  int* endThreshold = NULL;

  mafFPR = (double**)malloc(numMAF*sizeof(double*));
  mafRatio = (double**)malloc(numMAF*sizeof(double*));
  endThreshold = (int*)malloc(numMAF*sizeof(int));

  int pp = 0;
  for (pp = 0; pp < numMAF; pp++) {
    mafFPR[pp] = (double*)malloc(numThreshold*sizeof(double));
    mafRatio[pp] = (double*)malloc(numThreshold*sizeof(double));
  }

  readFPRRatio(fprdir, numMAF, numThreshold, mafFPR, mafRatio, endThreshold, delims);

  double minPvalueObserved = readPvalFile(pvalf); 
  printf("observed min pvalue = %e\n",minPvalueObserved);

  double tempps = 1.0-minPvalueObserved/2.0; // zscore gives identical answer as minp does. 
  double maxBetaObserved = 0;
  vdCdfNormInv(1, &tempps, &maxBetaObserved); // use inv(norm), and not inv(t-density), okay because the inv(norm) is 1-to-1 and monotonic
  printf ("observed max abs zscore:%f\n",maxBetaObserved);
  
  double maxGrdnObserved = readPvalFile(grdf); // read permutate max llh 
  printf ("observed max likelihood:%f\n",maxGrdnObserved);
  
  sampleMVN(maxBetaObserved, minPvalueObserved, maxGrdnObserved, seed, corrWindow, numSample, mafInt, mafFPR, mafRatio, numThreshold, endThreshold, outf, noPvalCor, numWindow, numSNPInWindow, doExPvalCor, corFactor, snp_weights,optimize);

  free(mafInt);
  free(endThreshold);
  free(correlation);
  for (pp = 0; pp < numMAF; pp++) {
    free(mafFPR[pp]);
    free(mafRatio[pp]);
  }
  free(mafFPR);
  free(mafRatio);
  for (pp = 0; pp < numWindow; pp++) {
    free(corrWindow[pp]);
  }
  free(corrWindow);
  free(numSNPInWindow);
	
	clock_t toc = clock();

  printf("\nElapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	
  return 0;
}


void splitCorrelationMatrix(double* correlation, double** corrWindow, int* snpPos, int numWindow, int numSNP, int* numSNPInWindowArray) {
  int i,j,k, m;
  int startPos = -1;
  int startIndex = -1;
  int endIndex = -1;
  int numSNPInWindow = 0;
  int windowIndex = -1;
  for (i = 0; i < numSNP; i++) {
    int pos = snpPos[i];
    if (i == 0) {
      startPos = pos;
      startIndex = i;
      numSNPInWindow++;
      windowIndex = 0;
    }
    else {
      if (pos >= startPos+windowSize) {
	endIndex = i-1;
	corrWindow[windowIndex] = (double*)malloc(sizeof(double)*numSNPInWindow*numSNPInWindow);
	numSNPInWindowArray[windowIndex] = numSNPInWindow;
	m = 0;
	for (j = startIndex; j <= endIndex; j++) {
	  for (k = startIndex; k <= endIndex; k++) {
	    int index = numSNP*j + k;
	    corrWindow[windowIndex][m] = correlation[index];
	    m++;
	  }
	}
	// printf("window %d, startPos = %d, startIndex = %d, endIndex = %d, # of SNPs in window = %d\n",windowIndex,startPos,startIndex,endIndex,numSNPInWindow);
	/*
	  char idbuf[SZ_LONG_BUF];
	  sprintf(idbuf,"temp/matrix_%d.txt",windowIndex);
	  save_matrix(corrWindow[windowIndex], numSNPInWindow, numSNPInWindow, 7, idbuf);
	*/
	startPos = pos;
	startIndex = i;
	numSNPInWindow = 1;
	windowIndex++;
      }
      else {
	numSNPInWindow++;
      }
    }
  }
  endIndex = numSNP-1;
  corrWindow[windowIndex] = (double*)malloc(sizeof(double)*numSNPInWindow*numSNPInWindow);
  numSNPInWindowArray[windowIndex] = numSNPInWindow;
  m = 0;
  for (j = startIndex; j <= endIndex; j++) {
    for (k = startIndex; k <= endIndex; k++) {
      int index = numSNP*j + k;
      corrWindow[windowIndex][m] = correlation[index];
      m++;
    }
  }
  printf("window %d, startPos = %d, startIndex = %d, endIndex = %d, # of SNPs in window = %d\n",windowIndex,startPos,startIndex,endIndex,numSNPInWindow);
	
  /*
    char idbuf[SZ_LONG_BUF];
    sprintf(idbuf,"temp/matrix_%d.txt",windowIndex);
    printf("window %d, startPos = %d, startIndex = %d, endIndex = %d, # of SNPs in window = %d\n",windowIndex,startPos,startIndex,endIndex,numSNPInWindow);
    save_matrix(corrWindow[windowIndex], numSNPInWindow, numSNPInWindow, 7, idbuf);
    exit(-1);
  */
}

double extractPermutatePval(char* file) {
// only use to read in the file created by the ./permutation 
  FILE* fp = fopen(file, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open file %s\n",file);
    exit(-1);
  }
  char line[SZ_LONG_BUF];
  fgets(line, SZ_LONG_BUF, fp );
  char* token = NULL;
  char* next;
  token = strtok(line,DEFAULT_DELIMS); 
  token = strtok(NULL,DEFAULT_DELIMS);// second 
  token = strtok(NULL,DEFAULT_DELIMS);// third
  token = strtok(NULL,DEFAULT_DELIMS);// four one is the number we want. 
  double pval = strtod(token, &next);
  return pval;
}

double readPvalFile(char* file) {
// only use to read in the file created by the ./permutation 
  FILE* fp = fopen(file, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open file %s\n",file);
    exit(-1);
  }
  char line[SZ_LONG_BUF];
  fgets(line, SZ_LONG_BUF, fp );
  char* token = NULL;
  char* next;
  token = strtok(line,DEFAULT_DELIMS);
  double pval = strtod(token, &next);
  return pval;
}

double getMinPvalue(int n, double* pvalue) {
  int i = 0;
  double minp = 1;
  for (i = 0; i < n; i++) {
    if (pvalue[i] < minp) { // takes a vector or something, and then iterator over. 
      minp = pvalue[i];
    }
  }
  return minp;
}

int readMAFFile(char* file, int* mafInt, int* snpPos, char* delims) {
  FILE* fp = fopen(file, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open file %s\n",file);
    exit(-1);
  }
  char line[SZ_LONG_BUF];
  fgets(line, SZ_LONG_BUF, fp ); // header
  char* token = NULL;
  token = strtok(line,delims); // # of tissues
  token = strtok(NULL,delims); // # of individuals
  numInd = atoi(token);
	 
  int snpIndex = 0;
  int startPos = -1;
  int numWindow = 1;
	
	char* next; // to convert to decimal 

  while (fgets(line, SZ_LONG_BUF, fp ) != NULL) {
    /*if (snpIndex > numSNP) {
      printf("expected %d SNPs, but observed more SNPs\n", numSNP);
      exit(-1);
    }*/
    token = strtok(line,delims); // rsid 
		token = strtok(NULL,delims); // chr
    token = strtok(NULL,delims); // pos
			
    int pos = atoi(token); // convert this token into int. 
    snpPos[snpIndex] = pos;
    if (snpIndex == 0) {
      startPos = pos;
    }
    else {
      if (pos > startPos+windowSize) {
				startPos = pos;
				numWindow++;
      }
    }

		double genoSum = 0; 
		int i = 0; 
    for (i = 0; i < numInd; i++) { // for each person. we go over their SNPs info 
      token = strtok(NULL,delims);
			double geno = strtod(token, &next);
			genoSum += geno;
			
    }
	
    double tempmaf = (double)genoSum/(2.0*(double)numInd); // compute Minor-allele frequency 
    	if (tempmaf > 0.5) {
      tempmaf = 1.0 - tempmaf;
    }
    mafInt[snpIndex] = (int)round(tempmaf*100.0); // store SNP info here. 
    if (mafInt[snpIndex] < 1) {
	   printf("%dth SNP has MAF < 1\n",snpIndex);
      mafInt[snpIndex] = 1;
    }
    snpIndex++; // which snp are we looking at 
  }
  fclose(fp);
  return numWindow;
}

int readCorrelation(char* corrf, int numSNP, double* correlation, char* delims) {
  //gzFile gzfp = gzopen(corrf,"r");
  FILE* corrfp = fopen(corrf, "r");
  //if (gzfp == NULL) {
  if (corrfp == NULL) {
    printf("ERROR: cannot open snps correlation file %s\n",corrf);
    exit(-1);
  }
  char cholline[SZ_LONG_BUF];
  int snpIndex = 0;
  char* next;
  //while (gzgets(gzfp, cholline, SZ_LONG_BUF) != NULL) {
  while (fgets(cholline, SZ_LONG_BUF, corrfp ) != NULL) {
    /*if (snpIndex >= numSNP) {
      printf("# of SNPss > %d\n",numSNP);
      exit(-1);
    }*/
    int i = 0;
    for (i = 0; i < numSNP-snpIndex; i++) {
      //for (i = 0; i < numSNP; i++) {
      char* token = NULL;
      if (i == 0) {
        token = strtok(cholline,delims);
      }
      else {
        token = strtok(NULL,delims);
      }
      //int index = i*numSNP+snpIndex;
      //correlation[index] = strtod(token, &next);
      int index1 = snpIndex*numSNP + i*numSNP+ snpIndex;
      double value = strtod(token, &next);
      correlation[index1] = value;
      if (i != 0) {
        int index2 = snpIndex*numSNP + snpIndex + i;
	correlation[index2] = value;
      }
    }
    snpIndex++;
  }
  int i = 0;
  for (i = 0; i < numSNP; i++) {
    correlation[i*numSNP+i] = correlation[i*numSNP+i]+EPS;
  }
  //save_matrix(correlation, numSNP, numSNP, 7, "temp/cor.matrix.txt");
  fclose(corrfp);
  //gzclose(gzfp);
  return 0;
}

int readNumThreshold(char* fprdir) {
  char file[SZBUF];
  sprintf(file,"%s/maf_1.txt", fprdir);
  FILE* fp = fopen(file, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open file %s\n",file);
    exit(-1);
  }
  char line[SZ_LONG_BUF];
  fgets(line, SZ_LONG_BUF, fp );
  int numLine = 0;
  while (fgets(line, SZ_LONG_BUF, fp ) != NULL) {
    numLine++;
  }
  fclose(fp);
  printf("# of thresholds = %d\n",numLine);
  return numLine;
}

void readFPRRatio(char* fprdir, int numMAF, int numThreshold, double** mafFPR, double** mafRatio, int* endThreshold, char* delims) {
  int i = 0;
  for (i = 0; i < numMAF; i++) {
    char file[SZBUF];
    sprintf(file,"%s/maf_%d.txt", fprdir,(i+1));
    FILE* fp = fopen(file, "r");
    if (fp == NULL) {
      printf("ERROR: cannot open file %s\n",file);
      exit(-1);
    }
    char line[SZ_LONG_BUF];
    // !!! CHAGEN: ignore header lines starting with #
    //fgets(line, SZ_LONG_BUF, fp );
    int j = 0;
    endThreshold[i] = numThreshold - 1;
    for (j = 0; j < numThreshold; j++) {
      fgets(line, SZ_LONG_BUF, fp );
      char* token = NULL;
      char* next;
      token = strtok(line,delims);
      double threshold = strtod(token, &next);
      token = strtok(NULL,delims);
      double fpr = strtod(token, &next);
      double ratio = 0;
      mafFPR[i][j] = fpr;
      token = strtok(NULL,delims);
      if (strcmp(token,"NA") != 0) {
        ratio = strtod(token, &next);
      }
      else {
      	endThreshold[i] = j - 1;
      	mafRatio[i][j] = 0;
      	break;
      }
      if (ratio < 0) {
        endThreshold[i] = j - 1;
        mafRatio[i][j] = 0;
        break;
      }
      else {
        mafRatio[i][j] = ratio;
      }
    }
    fclose(fp);
    //printf("%d %d\n",i,endThreshold[i]);
  }
}

int compare (const void * a, const void * b)
{
  if (*(double*)a - *(double*)b < 0) {
    return -1;
  }
  return 1;
}

void sampleMVN(double maxBetaObserved,double minPvalueObserved, double maxGrdnObserved, int seed, double** corrWindow, int numSample, int* mafInt, double** mafFPR, double** mafRatio, int numThreshold, int* endThreshold, char* outf, int noPvalCor, int numWindow, int* numSNPInWindow, int doExPvalCor, double corFactor, float *snp_weights, int optimize) {

  int info = 0;
  int i = 0;
  int pp = 0;
  int j = 0;
  int m = 0;
  int n = 0;
  for (j = 0; j < numWindow; j++) {
    dpotrf("U",&numSNPInWindow[j],corrWindow[j],&numSNPInWindow[j],&info);
  }
  //save_matrix(corrWindow[0], numSNPInWindow[0], numSNPInWindow[0], 7, "temp/cor.matrix.txt");
  //exit(-1);

  VSLStreamStatePtr stream;
  vslNewStream(&stream, VSL_BRNG_MT19937, time(NULL)+seed); //time(NULL)+
  //vslNewStream(&stream, VSL_BRNG_MT19937, seed);

  int numSamplePerBatch = 1000;
  //int numSamplePerBatch = 1;
  int numBatch = numSample / numSamplePerBatch;
  if (numSample < numSamplePerBatch) {
    numBatch = 1;
  }
  if (numSample < numSamplePerBatch) {
    printf("Please run at least %d sampling\n",numSamplePerBatch);
    exit(-1);
  }
  if (numSample % numSamplePerBatch != 0) {
    printf("Please run a multiple of %d sampling\n",numSamplePerBatch);
    exit(-1);
  }
  // create zscore vector for each window
  double** zscoreWindow = (double**)malloc(sizeof(double*)*numWindow);
  double** pvalueWindow = (double**)malloc(sizeof(double)*numWindow);
 
  for (j = 0; j < numWindow; j++) {
    zscoreWindow[j] = (double*)malloc(sizeof(double)*numSNPInWindow[j]);
    pvalueWindow[j] = (double*)malloc(sizeof(double)*numSNPInWindow[j]);
  }
  int numSigSample = 0;
  int numSigSample_mvn = 0;
  //int numSigSample_beta = 0; 
  
  printf("# of batch = %d\n", numBatch);
  printf("# of samples per batch = %d\n", numSamplePerBatch);

  double* allMinPvalue = NULL;
  double* allMaxLR = NULL; 
  
  int sampleIndex = 0;
  //if (doExPvalCor) {
    allMinPvalue = (double*)malloc(sizeof(double)*numSample);
    allMaxLR = (double*)malloc(sizeof(double)*numSample);
  //}

  for (pp = 0; pp < numBatch; pp++) {
    double* minPvalue = (double*)malloc(sizeof(double)*numSamplePerBatch);
    double* maxLR = (double*)malloc(sizeof(double)*numSamplePerBatch);
    //double* maxBeta = (double*)malloc(sizeof(double)*numSamplePerBatch);
    
    for (j = 0; j < numSamplePerBatch; j++) {
      minPvalue[j] = 1;
      maxLR[j] = -INFINITY ; 
      //maxBeta[j] = -INFINITY; 
    }

    for (j = 0; j < numWindow; j++) {
      double* zscoreBatch = (double*)malloc(sizeof(double)*numSamplePerBatch*numSNPInWindow[j]);
      int batchIndex = 0;
      for (i = 0; i < numSamplePerBatch; i++) {
        vdRngGaussian(VSL_METHOD_DGAUSSIAN_BOXMULLER, stream, numSNPInWindow[j], &zscoreBatch[i*numSNPInWindow[j]], 0., 1.);
      }
      dtrmm("L","U","T","N",&numSNPInWindow[j], &numSamplePerBatch, &onef, corrWindow[j], &numSNPInWindow[j], zscoreBatch, &numSNPInWindow[j]);
      batchIndex = 0;
      for (i = 0; i < numSamplePerBatch; i++) {
        for (m = 0; m < numSNPInWindow[j]; m++) {
          zscoreWindow[j][m] = -1.0 * fabs(zscoreBatch[batchIndex]); // abs then negative so all the zscore are in negative 
          // printf("unchanged z score %f\n",zscoreWindow[j][m]);
          batchIndex++;
        }
        vdCdfNorm(numSNPInWindow[j], zscoreWindow[j], pvalueWindow[j]);

        for (m = 0; m < numSNPInWindow[j]; m++) {
          double pvalue = 2.0*pvalueWindow[j][m];
          int snpIndex = m;
          for (n = 0; n < j; n++) {
            snpIndex += numSNPInWindow[n];
          }
          int mafIndex = mafInt[snpIndex] - 1; // MAF 1% becomes index of 0, 2% becomes index of 1, so on
          
          double pvalueCor = pvalue;
          double this_beta = fabs(zscoreWindow [j][m]) ; // zscores are standardized beta. 
          double this_gradient = 1.0 ;
            pvalueCor = correctPvalue(pvalue, mafFPR[mafIndex], mafRatio[mafIndex], endThreshold[mafIndex]);
            double tempps = 1.0-pvalueCor/2.0;
            double tempzscore = 0;
            vdCdfNormInv(1, &tempps, &tempzscore);
            this_beta = tempzscore;
          
          if (pvalueCor < minPvalue[i]) {
            minPvalue[i] = pvalueCor;
          }
          /*if (this_beta > maxBeta[i]) {
            maxBeta[i] = this_beta;
          }*/

          if (optimize==1){
            // sup power. so we need to set the mu as beta_hat 
            this_gradient = computeGradient(this_beta,this_beta,*(snp_weights+snpIndex)) ; // STANDARDIZED BETA, SD=1 
          }
          else if (optimize==0){
            this_gradient = computeGradient(this_beta,MEAN,*(snp_weights+snpIndex)) ; // STANDARDIZED BETA, SD=1 
          }
          else {
            printf("\noptimize indication not recognized.\n");
          }
          if ( this_gradient > maxLR[i]){
            maxLR[i] = this_gradient; 
          }
        }
      }
        free(zscoreBatch);
      }

      for (j = 0; j < numSamplePerBatch; j++) {
        //if (doExPvalCor) {
          allMinPvalue[sampleIndex] = minPvalue[j];
          allMaxLR[sampleIndex] = maxLR[j]; // max llh 
          sampleIndex++; 
        //}
        if (maxLR[j] > maxGrdnObserved) {
          numSigSample_mvn++;
        }
        
        if (minPvalue[j] < minPvalueObserved) {
          numSigSample++;
        }
        /*if (maxBeta[j] > maxBetaObserved) {
          numSigSample_beta++;
        }*/
      }

      if (numBatch >= 10) {
        if (pp % (numBatch/10) == 0) {
          printf("%d batch done\n",pp);
        }
      }
      free(minPvalue);
  }

  qsort(allMinPvalue, numSample, sizeof(double), compare);
  int rank = (int)round((double)numSample*minPvalueObserved);
  if (rank == 0) {
    rank = 1;
  }
  double samplePvalue = allMinPvalue[rank-1];
  double samplePvalue_uncor = (double)numSigSample/(double)numSample;
  FILE* tempfp = fopen(outf, "w");
  fprintf(tempfp, "using_minp_uncorrected_pval %f corrected_pval %f\n", samplePvalue_uncor, samplePvalue_uncor*corFactor);
	// write to console 
  printf("using_minp_uncorrected_pval %-.*lg\t corrected_pval %f\n", DEFAULT_NDIGITS,samplePvalue_uncor, samplePvalue_uncor*corFactor);

  
  double samplePvalue_mvn_uncor = (double)numSigSample_mvn/(double)numSample;
  fprintf(tempfp, "using_lr_uncorrected_pval %f corrected_pval %f\n",samplePvalue_mvn_uncor, samplePvalue_mvn_uncor*corFactor);
  fclose(tempfp);
	// write to console 
  printf("using_lr_uncorrected_pval %-.*lg\t corrected_pval %f\n", DEFAULT_NDIGITS,samplePvalue_mvn_uncor, samplePvalue_mvn_uncor*corFactor);
  
  for (j = 0; j < numWindow; j++) {
    free(zscoreWindow[j]);
    free(pvalueWindow[j]);
  }
  free(zscoreWindow);
  free(pvalueWindow);
  free(allMinPvalue);
}


double correctPvalue(double pvalue, double* fpr, double* ratio, int endThreshold) {
  int closestIndex = binarySearch(pvalue, fpr, endThreshold);
  double correctionFactor = ratio[closestIndex];
  if (correctionFactor == 0) {
    printf("ratio is 0\n");
    exit(-1);
  }
  double pvalueCor = pvalue * correctionFactor;
  //printf("%f %f %f %f\n", pvalue, fpr[closestIndex], ratio[closestIndex], correctionFactor);
  return pvalueCor;
} 


int minDiff(double value, double* A, int* candidates, int numCandidate) {
  int minIndex = -1;
  double minDiffValue = 1;
  int i = 0;
  for (i = 0; i < numCandidate; i++) {
    double diff = fabs(A[candidates[i]]-value);
    if (diff < minDiffValue) {
      minDiffValue = diff;
      minIndex = candidates[i];
    }
  }
  return minIndex;
}

// not sure what is this binary search ?? is this simple search on a vector ? 
int binarySearch(double value, double* A, int imax) {
  int length = imax+1;
  int imin = 0;
  while (imax > imin) {
    int imid = (imax+imin)/2;
    if (imid == imin) {
      break;
    }
    if (value < A[imid]) {
      imin = imid + 1;
    }
    else if (value > A[imid]) {
      imax = imid - 1;
    }
    else {
      return imid;
    }
  }
  if (imax != imin && imax != (imin+1) ) {
    printf("imax != imin && imax != (imin+1)\n");
    exit(-1);
  }
  int candidates[4];
  int numCandidate = 0;
  if (imin != 0) {
    candidates[numCandidate] = imin-1;
    numCandidate++;
  }
  candidates[numCandidate] = imin;
  numCandidate++;
  if (imax == imin && imin != length-1) {
    candidates[numCandidate] = imin+1;
    numCandidate++;
  }
  else {
    candidates[numCandidate] = imax;
    numCandidate++;
    if (imax != length-1) {
      candidates[numCandidate] = imax+1;
      numCandidate++;
    }
  }
  int minIndex = minDiff(value, A, candidates, numCandidate);
  return minIndex;
}

void save_matrix(double* K, int nr, int nc, int ndigits, char* filename) {
   FILE* fp = fopen(filename, "wb");
  int k,l;
  for (k = 0; k < nr; ++k) {
    for (l = 0; l < nc; ++l) {
      if (l == nc - 1) {
	fprintf(fp,"%-.*lg\n",ndigits,K[k+l*nr]);
      }
      else {
	fprintf(fp,"%-.*lg\t",ndigits,K[k+l*nr]);
      }
    }
  }
  fclose(fp);
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

double gammln(double xx) {
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

double tcdf (double t, double nu) {
  if ( isnan(t) ) return 1.;
  else return betai(nu/2.,0.5,nu/(nu+t*t));
}

// for matrix S = I-11'/n, and K
// compute centered trace tr(SKS)
// this is equivalent to sum(diag(K))-sum(K)/nrow(K)
double centered_trace(double* kin, int n) {
  int i,j;
  double dsum,asum;
  dsum = asum = 0;
  for(i=0; i < n; ++i) {
    for(j=0; j < n; ++j) {
      asum += kin[i+n*j];
      if ( i == j ) {
	dsum += kin[i+n*j];
      }
    }
  }
  return (dsum - asum/(double)n);
}

void print_help(void) {
  fprintf(stderr,"\nPlease visit github.com/datduong/FUNC-eGene\n");
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
  //info = clapack_dgetri (CblasColMajor, n, Y, n, ipiv);
  dgetri(&n,Y,&n,ipiv,work,&lwork,&info);
  
  /*  Cleanup and exit  */
  free(ipiv);
  free(work);
  return info;
}

double computeGradient (double testStatistics, double mean, double wts ) {
    double part1 = 0; 
    part1 = log( 0.5 * ( exp (-1.0*testStatistics*mean) + exp(testStatistics*mean) ) ) - 0.5*mean*mean;
    return ( log(wts) + part1 );
}

// return the pointer to the array ? yes. NO OTHER WAY, see the c-guide . 
float * read_wts ( char* file, float* arr, int stop ){
	FILE *fp; 
	if((fp=fopen(file,"r"))==NULL){
			printf("cannot open the file");
			exit(1);
	}
	int i = 0; 
	for (i=0; i<stop; i++){
		fscanf(fp,"%f", &arr[i]); // the value will be store to array. 
	}
	return arr; 
}



double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / 
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    /*if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p 
           << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );
    }
    */
    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}

