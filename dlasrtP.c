/* dlasrt.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include <omp.h>

int lsame_(char *a, char *b){
	if(*a==*b) return 1;
	else return 0;
}

/* Subroutine */ int dlasrt_(char *id, int *n, double *d__, int *
	info)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    int i__, j;
    double d1, d2, d3;
    int dir;
    double tmp;
    int endd;
    // extern logical lsame_(char *, char *);
    int stack[64]	/* was [2][32] */;
    double dmnmx;
    int start;
    // extern /* Subroutine */ int xerbla_(char *, int *);
    int stkpnt;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  Sort the numbers in D in increasing order (if ID = 'I') or */
/*  in decreasing order (if ID = 'D' ). */

/*  Use Quick Sort, reverting to Insertion sort on arrays of */
/*  size <= 20. Dimension of STACK limits N to about 2**32. */

/*  Arguments */
/*  ========= */

/*  ID      (input) CHARACTER*1 */
/*          = 'I': sort D in increasing order; */
/*          = 'D': sort D in decreasing order. */

/*  N       (input) int */
/*          The length of the array D. */

/*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
/*          On entry, the array to be sorted. */
/*          On exit, D has been sorted into increasing order */
/*          (D(1) <= ... <= D(N) ) or into decreasing order */
/*          (D(1) >= ... >= D(N) ), depending on ID. */

/*  INFO    (output) int */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input paramters. */

    /* Parameter adjustments */
    --d__;

    /* Function Body */
    *info = 0;
    dir = -1;
    if (lsame_(id, "D")) {
	dir = 0;
    } else if (lsame_(id, "I")) {
	dir = 1;
    }
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	//xerbla_("DLASRT", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 1) {
	return 0;
    }

    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {

/*        Do Insertion sort on D( START:ENDD ) */

	if (dir == 0) {

/*           Sort into decreasing order */

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] > d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L30;
		    }
/* L20: */
		}
L30:
		;
	    }

	} else {

/*           Sort into increasing order */

	    i__1 = endd;
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
		i__2 = start + 1;
		for (j = i__; j >= i__2; --j) {
		    if (d__[j] < d__[j - 1]) {
			dmnmx = d__[j];
			d__[j] = d__[j - 1];
			d__[j - 1] = dmnmx;
		    } else {
			goto L50;
		    }
/* L40: */
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

/*        Partition D( START:ENDD ) and stack parts, largest one first */

/*        Choose partition entry as median of 3 */

	d1 = d__[start];
	d2 = d__[endd];
	i__ = (start + endd) / 2;
	d3 = d__[i__];
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

/*           Sort into decreasing order */

	    i__ = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (d__[j] < dmnmx) {
		goto L70;
	    }
L80:
	    ++i__;
	    if (d__[i__] > dmnmx) {
		goto L80;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	} else {

/*           Sort into increasing order */

	    i__ = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (d__[j] > dmnmx) {
		goto L100;
	    }
L110:
	    ++i__;
	    if (d__[i__] < dmnmx) {
		goto L110;
	    }
	    if (i__ < j) {
		tmp = d__[i__];
		d__[i__] = d__[j];
		d__[j] = tmp;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
	    } else {
		++stkpnt;
		stack[(stkpnt << 1) - 2] = j + 1;
		stack[(stkpnt << 1) - 1] = endd;
		++stkpnt;
		stack[(stkpnt << 1) - 2] = start;
		stack[(stkpnt << 1) - 1] = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }
    return 0;

/*     End of DLASRT */

} /* dlasrt_ */

void vecGene(double* A, int size) {
        srand(time(NULL));
        //猜测生成0-1的数
#pragma omp parallel for num_threads(2)
        for (int i = 0; i < size; i++) {
			double temp=rand()%100;//产生0-RAND_MAX的数
            double f=(double)rand()/RAND_MAX;//产生0-1的数
            A[i] = temp*f; //产生A[j][i]
        }
}


//合并两个区间
void merge(int l1, int r1, int r2, double* data, int* temp) {
    int top = l1, p = l1, q = r1;
    while (p < r1 || q < r2) {
        if (q >= r2 || (p < r1 && data[p] <= data[q])) {
            temp[top++] = data[p++];
        }
        else {
            temp[top++] = data[q++];
        }
    }
    for (top = l1; top < r2; top++) {
        data[top] = temp[top];
    }
}

void merge_sort(int l, int r, double* data, int N) {
    int i, j, t, *temp;
    temp = (int*)malloc(N * sizeof(int));
    //这里做了一些优化，预处理合并了单个的区间，略微提高的速度
    #pragma omp parallel for private(i, t) shared(N, data)
    for (i = 0; i < N/2; i++)
        if (data[i*2] > data[i*2+1]) {
            t = data[i*2];
            data[i*2] = data[i*2+1];
            data[i*2+1] = t;
        }

    //i代表每次归并的区间长度，j代表需要归并的两个区间中最小的下标
    for (i = 2; i < r; i *= 2) {
        #pragma omp parallel for private(j) shared(r, i)
        for (j = 0; j < r-i; j += i*2) {
            merge(j, j+i, (j+i*2 < r ? j+i*2 : r), data, temp);
        }
    }
}

void vecShow(double* A, int size) {
        for (int i = 0; i < size; i++) {
                printf("%f,",A[i]); //A[i]
        }
		printf("\n");
}


int main(){
	char *id;
	int *n;
	int *info;
	id=(char *)malloc(sizeof(char)*(1));
	n=(int *)malloc(sizeof(int)*(1));
	info=(int *)malloc(sizeof(int)*(1));
	*id='I';
	double *d__1;
	double *d__2;
	//测量时间的参数
	double start[2],stop[2];
	double time[2];

	printf("请输入带排序数组大小n：");
	scanf("%d",n);
	//生成随机数组
	d__1 =(double *)malloc(sizeof(double)*(*n+2));
	d__2 =(double *)malloc(sizeof(double)*(*n+2));
	vecGene(d__1,*n);
	for(int i=0;i<*n;i++){
		d__2[i]=d__1[i];
	}

	//原函数
	start[0]=omp_get_wtime();
	dlasrt_(id, n, d__1,info);
	stop[0]=omp_get_wtime();
	time[0]=stop[0]-start[0];

	start[1]=omp_get_wtime();
	merge_sort(0,*n,d__2,*n);
	stop[1]=omp_get_wtime();
	time[1]=stop[1]-start[1];
	vecShow(d__1,*n);
	vecShow(d__2,*n);
	printf("原函数时间：%f;并行归并函数时间：%f",time[0],time[1]);

	return 0;
}