/* dlasr.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include <omp.h>
#include <stdio.h>
#include<stdlib.h>
#include <time.h>
#include<string.h>
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
int numberthreads;

//生成对应的参数
float* ag;
float* atg;
float* c__g;
float* sg;


int lsame_(char *a, char *b){
	if(*a==*b) return 1;
	else return 0;
}

void printmatrix(float *t){
	for(int i=0;i<5;i++){
		for(int j=0;j<5;j++){
			printf("%f,",t[j+i*5]);
		}
		printf("\n");
	}
}

//测试用
void matGene(float* A, int lda , int m , int n) {
    srand(time(NULL));
// #pragma omp parallel for num_threads(10)
    for(int j = 0; j < lda; j++) {
        for (int i = 0; i < n; i++) {
            float temp=rand()%100;//产生0-RAND_MAX的数
            float f=(float)rand()/RAND_MAX;//产生0-1的数
            A[i * lda + j] = temp*f; //产生A[j][i]
        }
    }
}

void vecGene(float* A, int size) {
        srand(time(NULL));
        //猜测生成0-1的数
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < size; i++) {
                A[i] = (float)rand()/RAND_MAX; //A[i]
        }
}

void matShow(float* A, int size) {
        for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                        printf("%f,",A[i * size + j]);
                }
                printf("\n");
        }
}

void vecShow(float* A, int size) {
        for (int i = 0; i < size; i++) {
                printf("%f,",A[i]); //A[i]
        }
}
//测试用结束

/* Subroutine */ int dlasr_P(char *side, char *pivot, char *direct, int *m, 
	 int *n, float *c__, float *s, float *a, int *
	lda)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, info;
    float temp;
    float ctemp, stemp;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASR applies a sequence of plane rotations to a real matrix A, */
/*  from either the left or the right. */

/*  When SIDE = 'L', the transformation takes the form */

/*     A := P*A */

/*  and when SIDE = 'R', the transformation takes the form */

/*     A := A*P**T */

/*  where P is an orthogonal matrix consisting of a sequence of z plane */
/*  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R', */
/*  and P**T is the transpose of P. */

/*  When DIRECT = 'F' (Forward sequence), then */

/*     P = P(z-1) * ... * P(2) * P(1) */

/*  and when DIRECT = 'B' (Backward sequence), then */

/*     P = P(1) * P(2) * ... * P(z-1) */

/*  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation */

/*     R(k) = (  c(k)  s(k) ) */
/*          = ( -s(k)  c(k) ). */

/*  When PIVOT = 'V' (Variable pivot), the rotation is performed */
/*  for the plane (k,k+1), i.e., P(k) has the form */

/*     P(k) = (  1                                            ) */
/*            (       ...                                     ) */
/*            (              1                                ) */
/*            (                   c(k)  s(k)                  ) */
/*            (                  -s(k)  c(k)                  ) */
/*            (                                1              ) */
/*            (                                     ...       ) */
/*            (                                            1  ) */

/*  where R(k) appears as a rank-2 modification to the identity matrix in */
/*  rows and columns k and k+1. */

/*  When PIVOT = 'T' (Top pivot), the rotation is performed for the */
/*  plane (1,k+1), so P(k) has the form */

/*     P(k) = (  c(k)                    s(k)                 ) */
/*            (         1                                     ) */
/*            (              ...                              ) */
/*            (                     1                         ) */
/*            ( -s(k)                    c(k)                 ) */
/*            (                                 1             ) */
/*            (                                      ...      ) */
/*            (                                             1 ) */

/*  where R(k) appears in rows and columns 1 and k+1. */

/*  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is */
/*  performed for the plane (k,z), giving P(k) the form */

/*     P(k) = ( 1                                             ) */
/*            (      ...                                      ) */
/*            (             1                                 ) */
/*            (                  c(k)                    s(k) ) */
/*            (                         1                     ) */
/*            (                              ...              ) */
/*            (                                     1         ) */
/*            (                 -s(k)                    c(k) ) */

/*  where R(k) appears in rows and columns k and z.  The rotations are */
/*  performed without ever forming P(k) explicitly. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          Specifies whether the plane rotation matrix P is applied to */
/*          A on the left or the right. */
/*          = 'L':  Left, compute A := P*A */
/*          = 'R':  Right, compute A:= A*P**T */

/*  PIVOT   (input) CHARACTER*1 */
/*          Specifies the plane for which P(k) is a plane rotation */
/*          matrix. */
/*          = 'V':  Variable pivot, the plane (k,k+1) */
/*          = 'T':  Top pivot, the plane (1,k+1) */
/*          = 'B':  Bottom pivot, the plane (k,z) */

/*  DIRECT  (input) CHARACTER*1 */
/*          Specifies whether P is a forward or backward sequence of */
/*          plane rotations. */
/*          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1) */
/*          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1) */

/*  M       (input) int */
/*          The number of rows of the matrix A.  If m <= 1, an immediate */
/*          return is effected. */

/*  N       (input) int */
/*          The number of columns of the matrix A.  If n <= 1, an */
/*          immediate return is effected. */

/*  C       (input) DOUBLE PRECISION array, dimension */
/*                  (M-1) if SIDE = 'L' */
/*                  (N-1) if SIDE = 'R' */
/*          The cosines c(k) of the plane rotations. */

/*  S       (input) DOUBLE PRECISION array, dimension */
/*                  (M-1) if SIDE = 'L' */
/*                  (N-1) if SIDE = 'R' */
/*          The sines s(k) of the plane rotations.  The 2-by-2 plane */
/*          rotation part of the matrix P(k), R(k), has the form */
/*          R(k) = (  c(k)  s(k) ) */
/*                 ( -s(k)  c(k) ). */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The M-by-N matrix A.  On exit, A is overwritten by P*A if */
/*          SIDE = 'R' or by A*P**T if SIDE = 'L'. */

/*  LDA     (input) int */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
// //get suitable number of threads for rotations
    int num_threads = numberthreads; 
    /* Function Body */
    info = 0;
    // if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	// info = 1;
    // } else if (! (lsame_(pivot, "V") || lsame_(pivot, 
	//     "T") || lsame_(pivot, "B"))) {
	// info = 2;
    // } else if (! (lsame_(direct, "F") || lsame_(direct, 
	//     "B"))) {
	// info = 3;
    // } else if (*m < 0) {
	// info = 4;
    // } else if (*n < 0) {
	// info = 5;
    // } else if (*lda < max(1,*m)) {
	// info = 9;
    // }
    // if (info != 0) {
	// xerbla_("DLASR ", &info);
	// return 0;
    // }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
	//a_dim1>=m>=1
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			//a_dim1==1时，无法并行，因为a_dim1>=m>=1，若a_dim1==1则m==1,则无法执行到此，故a_dim1>1
			//可以完全并行
#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__2; ++i__) {
				temp = a[j + 1 + i__ * a_dim1];
				//temp=a[j+1][i];
				a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
					a[j + i__ * a_dim1];

				a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
					+ i__ * a_dim1];
/* L10: */
			}
            }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			//同上一情况
#pragma omp parallel for num_threads(num_threads) private(temp) 
			for (i__ = 1; i__ <= i__1; ++i__) {
				temp = a[j + 1 + i__ * a_dim1];
				a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
					a[j + i__ * a_dim1];
				a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
					+ i__ * a_dim1];
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			//a_dim1>j-1时可以完全并行，又因为a_dim1>=m>=1，j<=m，则a_dim1>=j,则a_dim1>j，则可以完全并行
			//当a_dim1<=j-1时，若(j-1)%ad!=0,则可以完全并行
			//否则串行执行（此处应该还可以进一步改进）
#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__2; ++i__) {
				temp = a[j + i__ * a_dim1];
				a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
					i__ * a_dim1 + 1];
				a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
					i__ * a_dim1 + 1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			//a_dim1>j-1时可以完全并行，又因为a_dim1>=m>=1，j<=m，则a_dim1>=j,则a_dim1>j，则可以完全并行
			//当a_dim1<=j-1时，若(j-1)%ad!=0,则可以完全并行
			//否则串行执行（此处应该还可以进一步改进）
			#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__1; ++i__) {
				temp = a[j + i__ * a_dim1];
				a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
					i__ * a_dim1 + 1];
				a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
					i__ * a_dim1 + 1];
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			//当j+ad>m的时候，可以完全并行，因为a_dim1>=m、j>=1，则j+a_dim1>m一定成立
			//当j+ad<=m且(m-j)%ad!=0时，可以完全并行
			//否则串行执行（此处应该还可以进一步改进）
			#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__2; ++i__) {
				temp = a[j + i__ * a_dim1];
				a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
					+ ctemp * temp;
				a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
					a_dim1] - stemp * temp;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
				i__1 = *n;
				//当j+ad>m的时候，可以完全并行，因为a_dim1>=m、j>=1，则j+a_dim1>m一定成立
				//当j+ad<=m且(m-j)%ad!=0时，可以完全并行
				//否则串行执行（此处应该还可以进一步改进）
#pragma omp parallel for num_threads(num_threads) private(temp)
				for (i__ = 1; i__ <= i__1; ++i__) {
					temp = a[j + i__ * a_dim1];
					a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
						+ ctemp * temp;
					a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
						a_dim1] - stemp * temp;
	/* L110: */
				}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			//理论最大可能线程数取决于a_dim1值线程数=a_dim1
			// num_threads=min(num_threads,a_dim1);
			// int tid;
			// #pragma omp parallel num_threads(num_threads) private(tid,temp,i__)
			// {
			// 	tid=omp_get_thread_num();
			// 	int i__start=tid+1;
			// 	while(i__start<(1+(j+1)*a_dim1)){
			// 		for(i__=i__start;i__<=i__2;i__=i__+a_dim1){
			// 			temp = a[i__ + (j + 1) * a_dim1];
			// 			a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
			// 				a[i__ + j * a_dim1];
			// 			a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
			// 				i__ + j * a_dim1];
			// 		}
			// 		i__start+=num_threads;
			// 	}
			// }
#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			//线程数取决于a_dim1值,同上
// 			num_threads=min(num_threads,a_dim1);
// 			int tid;
// #pragma omp parallel num_threads(num_threads) private(tid,temp,i__)
// 			{
// 				tid=omp_get_thread_num();
// 				int i__start=tid+1;
// 				while(i__start<(1+(j+1)*a_dim1)){
// 					for(i__=i__start;i__<=i__1;i__=i__+a_dim1){
// 						temp = a[i__ + (j + 1) * a_dim1];
// 						a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
// 							a[i__ + j * a_dim1];
// 						a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
// 							i__ + j * a_dim1];
// 					}
// 					i__start+=num_threads;
// 				}
// 			}

#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			//理论上对打可以同时启动(j-1)*a_dim1个线程
// 			num_threads=min(num_threads,(j-1)*a_dim1);
// 			int tid;
// #pragma omp parallel num_threads(num_threads) private(tid,temp,i__)
// 			{
// 				tid=omp_get_thread_num();
// 				int i__start=tid+1;
// 				while(i__start<(1+j*a_dim1)){
// 					for(i__=i__start;i__<=i__2;i__=i__+(j-1)*a_dim1){
// 						temp = a[i__ + j * a_dim1];
// 						a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
// 							i__ + a_dim1];
// 						a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
// 							a_dim1];
// 					}
// 					i__start+=num_threads;
// 				}
// 			}
#pragma omp parallel for num_threads(num_threads) private(temp)
		for (i__ = 1; i__ <= i__2; ++i__) {
			temp = a[i__ + j * a_dim1];
			a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				i__ + a_dim1];
			a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				a_dim1];
/* L170: */
		}
		    }
/* L180: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			//理论上对打可以同时启动(j-1)*a_dim1个线程
			// num_threads=min(num_threads,(j-1)*a_dim1);
			// int tid;
// #pragma omp parallel num_threads(num_threads) private(tid,temp,i__)
// 			{
// 				tid=omp_get_thread_num();
// 				int i__start=tid+1;
// 				while(i__start<(1+j*a_dim1)){
// 					for(i__=i__start;i__<=i__2;i__=i__+(j-1)*a_dim1){
// 						temp = a[i__ + j * a_dim1];
// 						a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
// 							i__ + a_dim1];
// 						a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
// 							a_dim1];
// 					}
// 					i__start+=num_threads;
// 				}
// 			}
#pragma omp parallel for num_threads(num_threads) private(temp)
		for (i__ = 1; i__ <= i__1; ++i__) {
			temp = a[i__ + j * a_dim1];
			a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				i__ + a_dim1];
			a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				a_dim1];
/* L170: */
		}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			//可同时启动(n-j)*ad个线程
// 			num_threads=min(num_threads,(*n-j)*a_dim1);
// 			int tid;
// #pragma omp parallel num_threads(num_threads) private(tid,temp,i__)
// 			{
// 				tid=omp_get_thread_num();
// 				int i__start=tid+1;
// 				while(i__start<(1+(*n)*a_dim1)){
// 					for(i__=i__start;i__<=i__2;i__=i__+(*n-j)*a_dim1){
// 						temp = a[i__ + j * a_dim1];
// 						a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
// 							+ ctemp * temp;
// 						a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
// 							a_dim1] - stemp * temp;
// 					}
// 					i__start+=num_threads;
// 				}
// 			}
#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			//可同时启动(n-j)*ad个线程
			// num_threads=min(num_threads,(*n-j)*a_dim1);
			// int tid;
// #pragma omp parallel num_threads(num_threads) private(tid,temp,i__)
// 			{
// 				tid=omp_get_thread_num();
// 				int i__start=tid+1;
// 				while(i__start<(1+(*n)*a_dim1)){
// 					for(i__=i__start;i__<=i__1;i__=i__+(*n-j)*a_dim1){
// 						temp = a[i__ + j * a_dim1];
// 						a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
// 							+ ctemp * temp;
// 						a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
// 							a_dim1] - stemp * temp;
// 					}
// 					i__start+=num_threads;
// 				}
// 			}
#pragma omp parallel for num_threads(num_threads) private(temp)
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of DLASR */

} /* dlasr_ */



//原函数
/* Subroutine */ int dlasr_(char *side, char *pivot, char *direct, int *m, 
	 int *n, float *c__, float *s, float *a, int *
	lda)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, info;
    float temp;
    float ctemp, stemp;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DLASR applies a sequence of plane rotations to a real matrix A, */
/*  from either the left or the right. */

/*  When SIDE = 'L', the transformation takes the form */

/*     A := P*A */

/*  and when SIDE = 'R', the transformation takes the form */

/*     A := A*P**T */

/*  where P is an orthogonal matrix consisting of a sequence of z plane */
/*  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R', */
/*  and P**T is the transpose of P. */

/*  When DIRECT = 'F' (Forward sequence), then */

/*     P = P(z-1) * ... * P(2) * P(1) */

/*  and when DIRECT = 'B' (Backward sequence), then */

/*     P = P(1) * P(2) * ... * P(z-1) */

/*  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation */

/*     R(k) = (  c(k)  s(k) ) */
/*          = ( -s(k)  c(k) ). */

/*  When PIVOT = 'V' (Variable pivot), the rotation is performed */
/*  for the plane (k,k+1), i.e., P(k) has the form */

/*     P(k) = (  1                                            ) */
/*            (       ...                                     ) */
/*            (              1                                ) */
/*            (                   c(k)  s(k)                  ) */
/*            (                  -s(k)  c(k)                  ) */
/*            (                                1              ) */
/*            (                                     ...       ) */
/*            (                                            1  ) */

/*  where R(k) appears as a rank-2 modification to the identity matrix in */
/*  rows and columns k and k+1. */

/*  When PIVOT = 'T' (Top pivot), the rotation is performed for the */
/*  plane (1,k+1), so P(k) has the form */

/*     P(k) = (  c(k)                    s(k)                 ) */
/*            (         1                                     ) */
/*            (              ...                              ) */
/*            (                     1                         ) */
/*            ( -s(k)                    c(k)                 ) */
/*            (                                 1             ) */
/*            (                                      ...      ) */
/*            (                                             1 ) */

/*  where R(k) appears in rows and columns 1 and k+1. */

/*  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is */
/*  performed for the plane (k,z), giving P(k) the form */

/*     P(k) = ( 1                                             ) */
/*            (      ...                                      ) */
/*            (             1                                 ) */
/*            (                  c(k)                    s(k) ) */
/*            (                         1                     ) */
/*            (                              ...              ) */
/*            (                                     1         ) */
/*            (                 -s(k)                    c(k) ) */

/*  where R(k) appears in rows and columns k and z.  The rotations are */
/*  performed without ever forming P(k) explicitly. */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          Specifies whether the plane rotation matrix P is applied to */
/*          A on the left or the right. */
/*          = 'L':  Left, compute A := P*A */
/*          = 'R':  Right, compute A:= A*P**T */

/*  PIVOT   (input) CHARACTER*1 */
/*          Specifies the plane for which P(k) is a plane rotation */
/*          matrix. */
/*          = 'V':  Variable pivot, the plane (k,k+1) */
/*          = 'T':  Top pivot, the plane (1,k+1) */
/*          = 'B':  Bottom pivot, the plane (k,z) */

/*  DIRECT  (input) CHARACTER*1 */
/*          Specifies whether P is a forward or backward sequence of */
/*          plane rotations. */
/*          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1) */
/*          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1) */

/*  M       (input) int */
/*          The number of rows of the matrix A.  If m <= 1, an immediate */
/*          return is effected. */

/*  N       (input) int */
/*          The number of columns of the matrix A.  If n <= 1, an */
/*          immediate return is effected. */

/*  C       (input) DOUBLE PRECISION array, dimension */
/*                  (M-1) if SIDE = 'L' */
/*                  (N-1) if SIDE = 'R' */
/*          The cosines c(k) of the plane rotations. */

/*  S       (input) DOUBLE PRECISION array, dimension */
/*                  (M-1) if SIDE = 'L' */
/*                  (N-1) if SIDE = 'R' */
/*          The sines s(k) of the plane rotations.  The 2-by-2 plane */
/*          rotation part of the matrix P(k), R(k), has the form */
/*          R(k) = (  c(k)  s(k) ) */
/*                 ( -s(k)  c(k) ). */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          The M-by-N matrix A.  On exit, A is overwritten by P*A if */
/*          SIDE = 'R' or by A*P**T if SIDE = 'L'. */

/*  LDA     (input) int */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    // if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	// info = 1;
    // } else if (! (lsame_(pivot, "V") || lsame_(pivot, 
	//     "T") || lsame_(pivot, "B"))) {
	// info = 2;
    // } else if (! (lsame_(direct, "F") || lsame_(direct, 
	//     "B"))) {
	// info = 3;
    // } else if (*m < 0) {
	// info = 4;
    // } else if (*n < 0) {
	// info = 5;
    // } else if (*lda < max(1,*m)) {
	// info = 9;
    // }
    // if (info != 0) {
	// xerbla_("DLASR ", &info);
	// return 0;
    // }

/*     Quick return if possible */
    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + 1 + i__ * a_dim1];
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + 1 + i__ * a_dim1];
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
/* L170: */
			}
		    }
/* L180: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (ctemp != 1. || stemp != 0.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of DLASR */

} /* dlasr_ */



//测试函数
void test(char *side, char *pivot, char *direct, int *m, 
	int *n, int *lda){
		float* a =(float *)malloc(sizeof(float)*((*n)*(*lda)+2));
		float* at=(float *)malloc(sizeof(float)*((*n)*(*lda)+2));
        float* c__ =(float *)malloc(sizeof(float)*((*n)+2));
        float* s =(float *)malloc(sizeof(float)*((*n)+2));

        //generate
        vecGene(c__,(*n));
        vecGene(s,(*n));
		//生成m*n的矩阵
        matGene(a,(*lda),(*m),(*n));

		for(int j=0;j<*m;j++){
			for(int i=0;i<*n;i++){
				at[j+i*(*lda)]=a[j+i*(*lda)];
			}
		}
		printf("%c,%c,%c:\n",*side,*pivot,*direct);
		float start[2],stop[2];
		float time[2];
		
		start[1]=omp_get_wtime();
		dlasr_(side, pivot, direct, m, n, c__, s, at, lda);
		stop[1]=omp_get_wtime();
		time[1]=stop[1]-start[1];

		start[0]=omp_get_wtime();
	    dlasr_P(side, pivot, direct, m, n, c__, s, a, lda);
		stop[0]=omp_get_wtime();
		time[0]=stop[0]-start[0];
		
		float imp=time[1]/time[0];
        // matShow(a,(*n));
		// printf("----------------------------------------------\n");
		// matShow(at,(*n));
		//计算最大误差
		int errorNums=0;
		for(int j=0;j<*m;j++){
			for(int i=0;i<*n;i++){
				if((at[j+i*(*lda)] - a[j+i*(*lda)])!=0){
					errorNums++;
				}
				//maxNums=max(maxError,fabs(at[j+i*(*n)] - a[j+i*(*n)]) / fabs(at[j+i*(*n)]));
			}
		}
		printf("误差数：%d\n",errorNums);
		printf("并行消耗时间：%f;串行消耗时间：%f;提升了%f倍\n",time[0],time[1],imp);
		printf("----------------------------------------------\n");
}


//测试函数
float test1(char *side, char *pivot, char *direct, int *m, int *n){
		// printf("%c,%c,%c:\n",*side,*pivot,*direct);
		float start[2],stop[2];
		float time[2];
		start[1]=omp_get_wtime();
		dlasr_(side, pivot, direct, m, n, c__g, sg, atg, m);
		stop[1]=omp_get_wtime();
		time[1]=stop[1]-start[1];

		start[0]=omp_get_wtime();
	    dlasr_P(side, pivot, direct, m, n, c__g, sg, ag, m);
		stop[0]=omp_get_wtime();
		time[0]=stop[0]-start[0];
		
		float imp=time[1]/time[0];
        // matShow(a,(*n));
		// printf("----------------------------------------------\n");
		// matShow(at,(*n));
		//计算最大误差
		// int errorNums=0;
		// for(int j=0;j<*m;j++){
		// 	for(int i=0;i<*n;i++){
		// 		if((at[j+i*(*lda)] - a[j+i*(*lda)])!=0){
		// 			errorNums++;
		// 		}
		// 		//maxNums=max(maxError,fabs(at[j+i*(*n)] - a[j+i*(*n)]) / fabs(at[j+i*(*n)]));
		// 	}
		// }
		// printf("误差数：%d\n",errorNums);
		// printf("线程数为：%d;并行消耗时间：%f;串行消耗时间：%f;提升了%f倍\n",numberthreads,time[0],time[1],imp);
		// printf("----------------------------------------------\n");
		return imp;
}

int main(int argc, char* argv[]){
	printf("请输入maxn的大小:\n");
	int maxn;
	//scanf("%d",&maxn);
	maxn=31000;
	printf("开始生成各项参数:\n");
	//生成对应的参数
	ag =(float *)malloc(sizeof(float)*(maxn*maxn+2));
	atg=(float *)malloc(sizeof(float)*(maxn*maxn+2));
	c__g =(float *)malloc(sizeof(float)*(maxn+2));
	sg =(float *)malloc(sizeof(float)*(maxn+2));

	//generate
	//printf("开始生成第一个一维……\n");
	vecGene(c__g,maxn);
	//printf("开始生成第二个一维……\n");
	vecGene(sg,maxn);
	//生成maxn*maxn的矩阵
	//printf("开始生成二维矩阵……\n");
	matGene(ag,maxn,maxn,maxn);
	memcpy(atg,ag,(maxn*maxn)*sizeof(int));
	
	char sidet[2]={'L','R'};
	char pivott[3]={'V','T','B'};
	char directt[2]={'F','B'};
	int mh[9]={10,100,1000,3000,6000 ,10000 ,15000 ,20000 ,30000};
	int ml[9]={10,100,1000,3000,6000 ,10000 ,15000 ,20000 ,30000};
	// int numberThreads[12]={2 ,6 ,10 ,14 ,18 ,22 ,26 ,30 ,32 ,34 ,36 ,40};
	//int numberThreads[12]={13,14,15,16,17,18,19};//缩小范围在16线程数周围进一步测试
	int procsNum=omp_get_num_procs();
	FILE *fp = NULL ;
	fp = fopen("testData.xls","w") ;
    fprintf(fp,"%s\t%s\t%s\t%s\t%s\t%s\t%s\n","row: m","col: n","number of threads","side","pivot","direct","speedup" ) ;

	//printf("开始测试:\n");
	//遍历测试测试情况
	for(int i=0;i<sizeof(sidet)/sizeof(char);i++){
		for(int j=0;j<sizeof(pivott)/sizeof(char);j++){
				for(int p=0;p<sizeof(mh)/sizeof(int);p++){
					for(int q=0;q<sizeof(ml)/sizeof(int);q++){
						//改造线程数，取1、2、4、8、……核数、核数*2
						for(int nt=1;nt<=procsNum*2;nt=nt*2)
						{
							numberthreads=nt;
							float t=test1(&sidet[i],&pivott[j],&directt[0],&mh[p],&ml[q]);
							fprintf(fp,"%d\t%d\t%d\t%c\t%c\t%c\t%f\n",mh[p],ml[q],nt,sidet[i],pivott[j],directt[0],t) ;
						}
					}
				}
			//}
		}
	}
	fclose(fp);
	return 0;
}
