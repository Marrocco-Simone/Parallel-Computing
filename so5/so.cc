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

void quick_sort(data_t *data, int low, int high)
{
    if (low >= high)
        return;

    int pi = partition(data, low, high);

    quick_sort(data, low, pi - 1);
    quick_sort(data, pi + 1, high);
}

void psort(int n, data_t *data)
{
    // FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
    // using the basic idea of quicksort.
    // std::sort(data, data + n);
    quick_sort(data, 0, n - 1);
}
