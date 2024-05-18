#include <algorithm>
#include <vector>
#define MIN_SORT 1024

using namespace std;

typedef unsigned long long data_t;

void merge_sort(data_t *data, int start, int end, int min_sort)
{
    if (end - start < min_sort)
    {
        std::sort(data + start, data + end);
        return;
    }
    int mid = (start + end) / 2;
#pragma omp task
    merge_sort(data, start, mid, min_sort);
#pragma omp task
    merge_sort(data, mid, end, min_sort);
#pragma omp taskwait
    std::inplace_merge(data + start, data + mid, data + end);
}

void psort(int n, data_t *data)
{
    // FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
    // using the basic idea of merge sort.
    // std::sort(data, data + n);
    int min_sort = n / MIN_SORT;
    if (min_sort < 2)
        min_sort = 2;

#pragma omp parallel
#pragma omp single
    merge_sort(data, 0, n, min_sort);
}
