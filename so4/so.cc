// #include <algorithm>
#include <vector>

using namespace std;

typedef unsigned long long data_t;

void merge(data_t *data, int start, int end, int mid)
{
    int i = start;
    int j = mid + 1;
    int k = start;
    vector<data_t> temp(end + 1);
    while (i <= mid && j <= end)
    {
        if (data[i] < data[j])
        {
            temp[k] = data[i];
            i++;
        }
        else
        {
            temp[k] = data[j];
            j++;
        }
        k++;
    }
    j = end;
    for (int h = mid; h >= i; h--, j--)
    {
        data[j] = data[h];
    }
    for (j = start; j < k; j++)
    {
        data[j] = temp[j];
    }
}

void merge_sort(data_t *data, int start, int end)
{
    if (start >= end)
        return;

    int mid = (start + end) / 2;
    merge_sort(data, start, mid);
    merge_sort(data, mid + 1, end);
    merge(data, start, end, mid);
}

void psort(int n, data_t *data)
{
    // FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
    // using the basic idea of merge sort.
    // std::sort(data, data + n);
    merge_sort(data, 0, n - 1);
}
