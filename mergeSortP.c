#include <omp.h>
#include <stdio.h>

//合并两个区间
void merge(int l1, int r1, int r2, int* data, int* temp) {
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

void merge_sort(int l, int r, int* data, int N) {
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