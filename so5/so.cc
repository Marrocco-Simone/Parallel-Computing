#include <algorithm>
#include <vector>
#define MIN_SORT 2048

using namespace std;

typedef unsigned long long data_t;

int partition(data_t *data, int low, int high)
{
    data_t pivot = data[high];
    int i = (low - 1);

    for (int j = low; j <= high; j++)
    {
        if (data[j] < pivot)
        {
            i++;
            swap(data[i], data[j]);
        }
    }
    swap(data[i + 1], data[high]);
    return (i + 1);
}

void quick_sort(data_t *data, int low, int high, int min_sort)
{
    if (low >= high)
        return;
    if (high - low < min_sort)
    {
        std::sort(data + low, data + high + 1);
        return;
    }

    int pi = partition(data, low, high);

#pragma omp task
    quick_sort(data, low, pi - 1, min_sort);
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
    if (n <= min_sort)
        min_sort = 1;

    int low = 0;
    int high = n - 1;

#pragma omp parallel
#pragma omp single
    quick_sort(data, low, high, min_sort);
}
