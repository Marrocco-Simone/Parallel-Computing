#include <algorithm>
#include <vector>
#include <climits>
#define MIN_SORT 128
#define STEP 10

using namespace std;

typedef unsigned long long data_t;

void merge(data_t *data, int start, int end)
{
    vector<data_t> temp(end - start);
    int indexes[STEP] = {0};
    int step = (end - start) / STEP + 1;
    for (int j = 0; j < end - start; j++)
    {
        int min_index = 0;
        data_t min_value = ULLONG_MAX;
        for (int i = 0; i < STEP; i++)
        {
            if (indexes[i] >= step)
                continue;
            if (start + indexes[i] + i * step >= end)
                continue;
            if (data[start + indexes[i] + i * step] >= min_value)
                continue;
            min_index = i;
            min_value = data[start + indexes[i] + i * step];
        }
        temp[j] = min_value;
        indexes[min_index]++;
    }
    for (int j = 0; j < end - start; j++)
    {
        data[j + start] = temp[j];
    }
}

void merge_sort(data_t *data, int start, int end)
{
    int step = (end - start) / STEP + 1;
    if (end - start < MIN_SORT)
    {
        std::sort(data + start, data + end);
        return;
    }

    for (int i = start; i < end; i += step)
    {
#pragma omp task
        merge_sort(data, i, min(i + step, end));
    }
#pragma omp taskwait
    merge(data, start, end);
}

void psort(int n, data_t *data)
{
// FIXME: Implement a more efficient parallel sorting algorithm for the CPU,
// using the basic idea of merge sort.
// std::sort(data, data + n);
#pragma omp parallel
#pragma omp single
    {
        merge_sort(data, 0, n);
    }
}
