#include <algorithm>
#include <vector>
#include <climits>
#define MIN_SORT 16
#define STEP 8

using namespace std;

typedef unsigned long long data_t;

void merge(data_t *data, int start, int end)
{
    vector<data_t> temp(end - start);
    vector<int> indexes(STEP, 0);
    int step = (end - start) / STEP + 1;
    for (int j = 0; j < end - start; j++)
    {
        int min_index = 0;
        data_t min_value = ULLONG_MAX;
        for (int i = 0; i < STEP; i++)
        {
            if (indexes[i] >= step || start + indexes[i] + i * step >= end || data[start + indexes[i] + i * step] >= min_value)
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
}

void merge_sort(data_t *data, int start, int end)
{
    int step = (end - start) / STEP + 1;
    if (end - start < MIN_SORT || step == 0)
    {
        std::sort(data + start, data + end);
        return;
    }

#pragma omp parallel for
    for (int i = start; i < end; i += step)
    {
        merge_sort(data, i, min(i + step, end));
    }
    merge(data, start, end);
}

void psort(int n, data_t *data)
{
    // FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
    // using the basic idea of merge sort.
    // std::sort(data, data + n);
    merge_sort(data, 0, n);
}
