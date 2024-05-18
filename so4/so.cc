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
    {
        vector<data_t> tmp;
        tmp.reserve(end - start);
        std::merge(data + start, data + mid, data + mid, data + end, tmp.begin());

        for (int i = start; i < end; i++)
            data[i] = tmp[i - start];
    }
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
