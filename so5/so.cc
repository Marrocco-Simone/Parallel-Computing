#include <algorithm>
#include <vector>
#include <cstdio>
#define MIN_SORT 2048

using namespace std;

typedef unsigned long long data_t;

int partition(data_t *data, int low, int high)
{
    data_t pivot = data[low];
    int i = low - 1, j = high + 1;

    while (true)
    {
        do
        {
            i++;
        } while (data[i] < pivot);

        do
        {
            j--;
        } while (data[j] > pivot);

        if (i >= j)
            return j;

        swap(data[i], data[j]);
    }
}

void quick_sort(data_t *data, int low, int high, int min_sort)
{
    if (high - low < min_sort)
    {
        std::sort(data + low, data + high + 1);
        return;
    }

    int pi = partition(data, low, high);

#pragma omp task
    quick_sort(data, low, pi, min_sort);
#pragma omp task
    quick_sort(data, pi + 1, high, min_sort);
}

void psort(int n, data_t *data)
{
    // FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
    // using the basic idea of quicksort.
    // std::sort(data, data + n);
    int min_sort = n / MIN_SORT;
    if (min_sort < 2)
        min_sort = 2;

    int low = 0;
    int high = n - 1;

#pragma omp parallel
#pragma omp single
    quick_sort(data, low, high, min_sort);
}
