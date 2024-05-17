#include <algorithm>
#include <vector>
#include <climits>
// #include <cstdio>
#define MIN_SORT 16
#define STEP 8

using namespace std;

typedef unsigned long long data_t;

void merge(data_t *data, int start, int end, int tabs)
{
    /* for (int t = 0; t < tabs; t++)
        printf("\t");
    printf("merge: %d %d\n", start, end);
    for (int t = 0; t < tabs; t++)
        printf("\t");
    printf("data: ");
    for (int t = start; t < end; t++)
        printf("%llu ", data[t]);
    printf("\n"); */
    vector<data_t> temp(end - start);
    vector<int> indexes(STEP, 0);
    int step = (end - start) / STEP + 1;
    for (int j = 0; j < end - start; j++)
    {
        int min_index = 0;
        data_t min_value = ULLONG_MAX;
        for (int i = 0; i < STEP; i++)
        {
            if (indexes[i] >= step || start + indexes[i] + i * step >= end)
                continue;
            if (data[start + indexes[i] + i * step] >= min_value)
                continue;
            min_index = i;
            min_value = data[start + indexes[i] + i * step];
        }
        temp[j] = data[start + indexes[min_index] + min_index * step];
        indexes[min_index]++;
    }
    for (int j = start; j < end; j++)
    {
        data[j] = temp[j - start];
    }
    /* for (int t = 0; t < tabs; t++)
        printf("\t");
    printf("data: ");
    for (int t = start; t < end; t++)
        printf("%llu ", data[t]);
    printf("\n"); */
}

void merge_sort(data_t *data, int start, int end, int tabs)
{
    int step = (end - start) / STEP + 1;
    /* for (int t = 0; t < tabs; t++)
        printf("\t");
    printf("step: %d\n", step); */
    if (end - start < MIN_SORT || step == 0)
    {
        std::sort(data + start, data + end);
        return;
    }

#pragma omp parallel for
    for (int i = start; i < end; i += step)
    {
        /*for (int t = 0; t < tabs; t++)
                 printf("\t");
            printf("i: %d, end: %d\n", i, min(i + step, end)); */
        merge_sort(data, i, min(i + step, end), tabs + 1);
    }
    merge(data, start, end, tabs);
}

void psort(int n, data_t *data)
{
    // FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
    // using the basic idea of merge sort.
    // std::sort(data, data + n);
    /* printf("n: %d\n", n);
    printf("data: ");
    for (int t = 0; t < n; t++)
        printf("%llu ", data[t]);
    printf("\n"); */
    merge_sort(data, 0, n, 0);
}
