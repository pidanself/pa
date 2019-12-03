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
int numProcs;

int lsame_(char *a, char *b){
	if(*a==*b) return 1;
	else return 0;
}


//生成随机数
double* vecGene(int size) {
	double* A=(double *)malloc(size*sizeof(double));
	srand(time(NULL));
	//生成0～size-1的随机数
	for (int i = 0; i < size; i++) {
		double temp=rand();//产生0-RAND_MAX的数
		double f=(double)rand()/RAND_MAX;//产生0-1的数
		A[i] = temp*f; //产生A[j][i]
	}
	return A;
}

void vecShow(double* A, int size) {
        for (int i = 0; i < size; i++) {
                printf("%f,",A[i]); //A[i]
        }
		printf("\n");
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

//原快排
void sort(int *a, int left, int right)
{
    if(left >= right)/*如果左边索引大于或者等于右边的索引就代表已经整理完成一个组了*/
    {
        return ;
    }
    int i = left;
    int j = right;
    int key = a[left];
     
    while(i < j)                               /*控制在当组内寻找一遍*/
    {
        while(i < j && key <= a[j])
        /*而寻找结束的条件就是，1，找到一个小于或者大于key的数（大于或小于取决于你想升
        序还是降序）2，没有符合条件1的，并且i与j的大小没有反转*/ 
        {
            j--;/*向前寻找*/
        }
         
        a[i] = a[j];
        /*找到一个这样的数后就把它赋给前面的被拿走的i的值（如果第一次循环且key是
        a[left]，那么就是给key）*/
         
        while(i < j && key >= a[i])
        /*这是i在当组内向前寻找，同上，不过注意与key的大小关系停止循环和上面相反，
        因为排序思想是把数往两边扔，所以左右两边的数大小与key的关系相反*/
        {
            i++;
        }
         
        a[j] = a[i];
    }
     
    a[i] = key;/*当在当组内找完一遍以后就把中间数key回归*/
    sort(a, left, i - 1);/*最后用同样的方式对分出来的左边的小组进行同上的做法*/
    sort(a, i + 1, right);/*用同样的方式对分出来的右边的小组进行同上的做法*/
                       /*当然最后可能会出现很多分左右，直到每一组的i = j 为止*/
}


void InsertSort(double *p, int low, int high)//指定区间插入排序，即对数组p的指定位置进行插入排序
{
	double temp;
	for (int i = low+1; i <= high; i++) {
		for (int j = i; (j > low) && (p[j] < p[j - 1]); j--) {
			temp = p[j];
			p[j] = p[j - 1];
			p[j - 1] = temp;
		}
	}
}


/*函数作用：取待排序序列中low、mid、high三个位置上数据，选取他们中间的那个数据作为枢轴*/
double SelectPivotMedianOfThree(double *arr, int low, int high)//三数取中
{
	double temp;
	int mid = low + ((high - low) >> 1);//计算数组中间的元素的下标  
 
										//使用三数取中法选择枢轴  
	if (arr[mid] > arr[high])//目标: arr[mid] <= arr[high]  
	{
		//swap(arr[mid], arr[high]);
		temp = arr[mid];
		arr[mid] = arr[high];
		arr[high] = temp;
	}
	if (arr[low] > arr[high])//目标: arr[low] <= arr[high]  
	{
		//swap(arr[low], arr[high]);
		temp = arr[low];
		arr[low] = arr[high];
		arr[high] = temp;
	}
	if (arr[mid] > arr[low]) //目标: arr[low] >= arr[mid]  
	{
		//swap(arr[mid], arr[low]);
		temp = arr[mid];
		arr[mid] = arr[low];
		arr[low] = temp;
	}
	//此时，arr[mid] <= arr[low] <= arr[high]  
	return arr[low];
	//low的位置上保存这三个位置中间的值  
	//分割时可以直接使用low位置的元素作为枢轴，而不用改变分割函数了  
}

void QuickSortAverage(double *p, int low, int high)//快排+三数取中+插入
{
	if (high - low + 1 < 20)
	{
		InsertSort(p, low, high);
		return;
	}//else时，正常执行快排
	int first = low;
	int last = high;
	//double key = p[first];/*用字表的第一个记录作为枢轴*/
	double key = SelectPivotMedianOfThree(p, low, high);
 
	while (first < last)
	{
		while (first < last && p[last] >= key)
		{
			--last;
		}
 
		p[first] = p[last];/*将比第一个小的移到低端*/
 
		while (first < last && p[first] <= key)
		{
			++first;
		}
 
		p[last] = p[first];
		/*将比第一个大的移到高端*/
	}
	p[first] = key;/*枢轴记录到位*/
	QuickSortAverage(p, low, first - 1);
	QuickSortAverage(p, first + 1, high);
}

int Partition(double * a, int low, int high)//分隔
{
	double pivotkey = SelectPivotMedianOfThree(a, low, high);
	while (low<high)
	{
		while (low<high && a[high] >= pivotkey)
			--high;
		a[low] = a[high];
		while (low<high && a[low] <= pivotkey)
			++low;
		a[high] = a[low];
	}
	//此时low==high 
	a[low] = pivotkey;
	return low;
}

//尝试将原函数进行并行
void QuickSortParallel(double *p, int low, int high)//2核快排
{
	//p[0] = BOUNDARY / 2;
	/*for (int i = low; i <= high; i++)
	{
		if (abs(p[i] - BOUNDARY / 2) < 10)
		{
			int temp = p[i];
			p[i] = p[0];
			p[0] = temp;
			break;
		}
	}*/
	int mid = Partition(p, low, high);
	#pragma omp parallel
	{
	#pragma omp sections
	{
	#pragma omp section
	{
		QuickSortAverage(p, low, mid-1);
	}
	#pragma omp section
	{
		QuickSortAverage(p, mid+1, high);
	}
	}
	}
}

//四核
void QuickSortParallel4Core(double *p, int low, int high)//4核快排
{
	int quarter1,quarter2;
	//p[0] = BOUNDARY / 2;
	/*for (int i = low; i <= high; i++)
	{
		if (abs(p[i] - BOUNDARY / 2) < 10)
		{
			int temp = p[i];
			p[i] = p[0];
			p[0] = temp;
			break;
		}
	}*/
	int mid = Partition(p, low, high);
	//p[low] = BOUNDARY / 4;
	//p[mid + 1] = BOUNDARY / 4 * 3;
	
	
#pragma omp parallel
	{

#pragma omp sections
	{
#pragma omp section
	{
		// printf("-1");
		quarter1 = Partition(p, low, mid - 1);
	}
#pragma omp section
	{
		// printf("0");
		quarter2 = Partition(p, mid + 1, high);
	}
	}
		
#pragma omp sections
	{
#pragma omp section
	{
		// printf("1");
		//double start1 = omp_get_wtime();
		QuickSortAverage(p, low, quarter1-1);
		//double end1 = omp_get_wtime();
		//printf("%lf\n", end1 - start1);
	}
#pragma omp section
	{
		// printf("2");
		//double start2 = omp_get_wtime();
		QuickSortAverage(p, quarter1 + 1, mid-1);
		//double end2 = omp_get_wtime();
		//printf("%lf\n", end2 - start2);
	}
#pragma omp section
	{
		// printf("3");
		//double start3 = omp_get_wtime();
		QuickSortAverage(p, mid+1, quarter2-1);
		//double end3 = omp_get_wtime();
		//printf("%lf\n", end3 - start3);
	}
#pragma omp section
	{
		// printf("4");
		//double start4 = omp_get_wtime();
		QuickSortAverage(p, quarter2+1, high);
		//double end4 = omp_get_wtime();
		//printf("%lf\n", end4 - start4);
	}
	}
	}

}



//合并两个区间
void merge(int l1, int r1, int r2, double* data, double* temp) {
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
    int i, j;
	double t;
	double *temp;
    temp = (double*)malloc(N * sizeof(double));
    //这里做了一些优化，预处理合并了单个的区间，略微提高的速度
    #pragma omp parallel for private(i, t) shared(N, data) num_threads(2*numProcs-1)
    for (i = 0; i < N/2; i++)
        if (data[i*2] > data[i*2+1]) {
            t = data[i*2];
            data[i*2] = data[i*2+1];
            data[i*2+1] = t;
        }

    //i代表每次归并的区间长度，j代表需要归并的两个区间中最小的下标
    for (i = 2; i < r; i *= 2) {
        #pragma omp parallel for private(j) shared(r, i) num_threads(2*numProcs-1)
        for (j = 0; j < r-i; j += i*2) {
            merge(j, j+i, (j+i*2 < r ? j+i*2 : r), data, temp);
        }
    }
}


double *test(int N){
	double *time=(double *)malloc(sizeof(double)*(4));
	char *id;
	int *n;
	int *info;
	id=(char *)malloc(sizeof(char)*(1));
	n=(int *)malloc(sizeof(int)*(1));
	info=(int *)malloc(sizeof(int)*(1));
	*id='I';
	//测量时间的参数
	double start[4],stop[4];
	*n=N;
	double *d__1;
	double *d__2;
	double *d__3;
	double *d__4;
	//生成随机数组
	d__2 =(double *)malloc(sizeof(double)*(*n+2));
	d__3=(double *)malloc(sizeof(double)*(*n+2));
	d__4=(double *)malloc(sizeof(double)*(*n+2));
	d__1=vecGene(*n);
	for(int i=0;i<*n;i++){
		d__2[i]=d__1[i];
	}
	for(int i=0;i<*n;i++){
		d__3[i]=d__1[i];
	}
	for(int i=0;i<*n;i++){
		d__4[i]=d__1[i];
	}
	//原函数
	start[0]=omp_get_wtime();
	dlasrt_(id, n, d__1,info);
	stop[0]=omp_get_wtime();
	time[0]=stop[0]-start[0];

	//并行归并排序
	start[1]=omp_get_wtime();
	merge_sort(0,*n,d__2,*n);
	stop[1]=omp_get_wtime();
	time[1]=stop[1]-start[1];
	// vecShow(d__1,*n);
	// vecShow(d__2,*n);

	//原函数并行（2线程）
	start[2]=omp_get_wtime();
	merge_sort(0,*n,d__2,*n);
	QuickSortParallel(d__3,0,*n-1);
	stop[2]=omp_get_wtime();
	time[2]=stop[2]-start[2];

	//原函数并行（4线程）
	start[3]=omp_get_wtime();
	QuickSortParallel4Core(d__4,0,*n-1);
	stop[3]=omp_get_wtime();
	time[3]=stop[3]-start[3];
	return time;
}

int main(){
	numProcs=omp_get_num_procs();
	double time[4];
	//打开xls文件
	FILE *fp = NULL ;
	fp = fopen("sortData.xls","w") ;
    fprintf(fp,"n           kind\t0\t1\t2\t3\n") ;
	//printf("原函数时间：%f;并行归并函数时间：%f;原函数并行（2线程）：%f;原函数并行（4线程）：%f\n",time[0],time[1],time[2],time[3]);
	int N[]={10,20,40,80,160,320,640,1000,2000,4000,8000,16000};//,20000,40000,80000,100000,500000,1000000,2000000,6000000,10000000,30000000,60000000,90000000,100000000,500000000};
	
	for(int i=0;i<sizeof(N)/sizeof(int);i++){
		fprintf(fp,"%d\t",N[i]) ;
		//初始化time
		for(int j=0;j<(sizeof(time)/sizeof(double));j++){
			time[j]=0;
		}
		//做三次计算
		for(int j=0;j<3;j++){
			double *temp=test(N[i]);
			for(int jj=0;jj<(sizeof(time)/sizeof(double));jj++){
				time[jj]+=temp[jj];
			}
		}
		//输出到xls表格
		fprintf(fp,"%f\t%f\t%f\t%f\n",time[0]/3,time[1]/3,time[2]/3,time[3]/3);
	}
	fclose(fp);
	return 0;
}