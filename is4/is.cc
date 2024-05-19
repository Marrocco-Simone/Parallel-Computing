#include <vector>
#define T 4
typedef double double4_t __attribute__((vector_size(T * sizeof(double))));
using namespace std;

constexpr double infty = 1.7e308;
constexpr int C = 3;

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[C];
    float inner[C];
};

int id4_t(int x, int y, int nx)
{
    return x + nx * y;
}

double color(int x, int y, int c, int nx, const float *data)
{
    return data[c + C * x + C * nx * y];
}

/** O(n^2)*/
void calculate_total_avg(int ny, int nx, const float *data, double4_t &total_avg)
{
#pragma omp parallel for
    for (int c = 0; c < C; c++)
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
                total_avg[c] += color(x, y, c, nx, data);
}

/** O(n^2) */
void calculate_sum_from_zero(int ny, int nx, const float *data, vector<double4_t> &sum_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
        {
            double4_t point_block = {color(x, y, 0, nx, data), color(x, y, 1, nx, data), color(x, y, 2, nx, data), 0.0};
            double4_t prev_left_block = {0.0};
            double4_t prev_up_block = {0.0};
            double4_t prev_left_up_block = {0.0};

            if (x != 0)
                prev_left_block = sum_from_zero[id4_t(x - 1, y, nx)];
            if (y != 0)
                prev_up_block = sum_from_zero[id4_t(x, y - 1, nx)];
            if (x != 0 && y != 0)
                prev_left_up_block = sum_from_zero[id4_t(x - 1, y - 1, nx)];

            sum_from_zero[id4_t(x, y, nx)] = (prev_left_block + prev_up_block - prev_left_up_block + point_block);
        }
}

/** O(n^2) */
void calculate_sum_squared_from_zero(int ny, int nx, const float *data, vector<double4_t> &sum_squared_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
        {
            double4_t point_block = {color(x, y, 0, nx, data) * color(x, y, 0, nx, data), color(x, y, 1, nx, data) * color(x, y, 1, nx, data), color(x, y, 2, nx, data) * color(x, y, 2, nx, data), 0.0};
            double4_t prev_left_block = {0.0};
            double4_t prev_up_block = {0.0};
            double4_t prev_left_up_block = {0.0};

            if (x != 0)
                prev_left_block = sum_squared_from_zero[id4_t(x - 1, y, nx)];
            if (y != 0)
                prev_up_block = sum_squared_from_zero[id4_t(x, y - 1, nx)];
            if (x != 0 && y != 0)
                prev_left_up_block = sum_squared_from_zero[id4_t(x - 1, y - 1, nx)];

            sum_squared_from_zero[id4_t(x, y, nx)] = (prev_left_block + prev_up_block - prev_left_up_block + point_block);
        }
}

/** O(1) */
int calculate_in_points(int x0, int x1, int y0, int y1)
{
    return (y1 - y0 + 1) * (x1 - x0 + 1);
}

/** O(1) - access all combinations of x1 / y1 / x0-1 / y0-1 */
void calculate_avg_in_color(int in_points, int x0, int x1, int y0, int y1, int nx, vector<double4_t> const &sum_from_zero, double4_t &inner)
{
    double4_t in_points4_t;
    for (int t = 0; t < T; t++)
    {
        in_points4_t[t] = in_points;
    }

    double4_t point_block = sum_from_zero[id4_t(x1, y1, nx)];
    double4_t prev_left_block = {0.0};
    double4_t prev_up_block = {0.0};
    double4_t prev_left_up_block = {0.0};

    if (x0 != 0 && x1 != 0)
        prev_left_block = sum_from_zero[id4_t(x0 - 1, y1, nx)];
    if (y0 != 0 && y1 != 0)
        prev_up_block = sum_from_zero[id4_t(x1, y0 - 1, nx)];
    if (x0 != 0 && x1 != 0 && y0 != 0 && y1 != 0)
        prev_left_up_block = sum_from_zero[id4_t(x0 - 1, y0 - 1, nx)];

    inner = (prev_left_up_block + point_block - prev_up_block - prev_left_block) / in_points4_t;
}

/** O(1) */
void calculate_avg_out_color(int in_points, int full_points, double4_t &inner, double4_t &total_avg, double4_t &outer)
{
    if (in_points == full_points)
    {
        return;
    }
    double4_t in_points4_t, out_points4_t;
    for (int t = 0; t < T; t++)
    {
        in_points4_t[t] = in_points;
        out_points4_t[t] = full_points - in_points;
    }
    double4_t in_block = inner * in_points4_t;
    outer = (total_avg - in_block) / out_points4_t;
}

/** O(1) - access all combinations of x1 / y1 / x0-1 / y0-1 */
void calculate_in_squared_sum(int x0, int x1, int y0, int y1, int nx, vector<double4_t> const &sum_squared_from_zero, double4_t &inner)
{
    double4_t point_block = sum_squared_from_zero[id4_t(x1, y1, nx)];
    double4_t prev_left_block = {0.0};
    double4_t prev_up_block = {0.0};
    double4_t prev_left_up_block = {0.0};

    if (x0 != 0 && x1 != 0)
        prev_left_block = sum_squared_from_zero[id4_t(x0 - 1, y1, nx)];
    if (y0 != 0 && y1 != 0)
        prev_up_block = sum_squared_from_zero[id4_t(x1, y0 - 1, nx)];
    if (x0 != 0 && x1 != 0 && y0 != 0 && y1 != 0)
        prev_left_up_block = sum_squared_from_zero[id4_t(x0 - 1, y0 - 1, nx)];

    inner = (prev_left_up_block + point_block - prev_up_block - prev_left_block);
}

/** O(1) */
void calculate_out_squared_sum(double4_t &inner, double4_t &end_sum_squared, double4_t &outer)
{
    outer = end_sum_squared - inner;
}

/** O(1) */
double calculate_in_error(double4_t &inner, double4_t &in_squared_sum, int in_points)
{
    double sum0 = in_squared_sum[0] - in_points * inner[0] * inner[0];
    double sum1 = in_squared_sum[1] - in_points * inner[1] * inner[1];
    double sum2 = in_squared_sum[2] - in_points * inner[2] * inner[2];
    return sum0 + sum1 + sum2;
}

/** O(1) */
double calculate_out_error(double4_t &outer, double4_t &out_squared_sum, int in_points, int full_points)
{
    int out_points = full_points - in_points;
    return calculate_in_error(outer, out_squared_sum, out_points);
}

/** O(1) */
void set_result(int x0, int x1, int y0, int y1, int nx, int ny, const std::vector<double4_t> &sum_from_zero, double4_t &total_avg, Result &result)
{
    double4_t outer = {0.0};
    double4_t inner = {0.0};

    int in_points = calculate_in_points(x0, x1, y0, y1);
    calculate_avg_in_color(in_points, x0, x1, y0, y1, nx, sum_from_zero, inner);
    calculate_avg_out_color(in_points, nx * ny, inner, total_avg, outer);

    result.x0 = x0;
    result.x1 = x1 + 1;
    result.y0 = y0;
    result.y1 = y1 + 1;
    for (int c = 0; c < C; c++)
    {
        result.inner[c] = inner[c];
        result.outer[c] = outer[c];
    }
}

double calculate_sq_error(int x0, int x1, int y0, int y1, int nx, int ny, double4_t &total_avg, double4_t &end_sum_squared, vector<double4_t> &sum_from_zero, vector<double4_t> &sum_squared_from_zero)
{
    double4_t outer = {0.0};
    double4_t inner = {0.0};

    int in_points = calculate_in_points(x0, x1, y0, y1);
    calculate_avg_in_color(in_points, x0, x1, y0, y1, nx, sum_from_zero, inner);
    calculate_avg_out_color(in_points, nx * ny, inner, total_avg, outer);

    double4_t outer_squared_sums = {0.0};
    double4_t inner_squared_sums = {0.0};
    calculate_in_squared_sum(x0, x1, y0, y1, nx, sum_squared_from_zero, inner_squared_sums);
    calculate_out_squared_sum(inner_squared_sums, end_sum_squared, outer_squared_sums);
    double sq_err = calculate_in_error(inner, inner_squared_sums, in_points) + calculate_out_error(outer, outer_squared_sums, in_points, nx * ny);

    return sq_err;
}

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/
Result segment(int ny, int nx, const float *data)
{
    Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
    double4_t total_avg = {0.0};
    calculate_total_avg(ny, nx, data, total_avg);
    vector<double4_t> sum_from_zero(nx * ny);
    calculate_sum_from_zero(ny, nx, data, sum_from_zero);
    vector<double4_t> sum_squared_from_zero(nx * ny);
    calculate_sum_squared_from_zero(ny, nx, data, sum_squared_from_zero);
    double4_t end_sum_squared = sum_squared_from_zero[id4_t(nx - 1, ny - 1, nx)];

    vector<int> best_solutions_coord(nx * ny * 4);
    vector<double> best_solutions(nx * ny);

#pragma omp parallel for schedule(dynamic, 2)
    for (int x0y0 = 0; x0y0 < ny * nx; x0y0++)
    {
        double min_err = infty;
        int i = x0y0;
        int y0 = x0y0 / nx;
        int x0 = x0y0 % nx;

        for (int y1 = y0; y1 < ny; y1++)
        {
            for (int x1 = x0; x1 < nx; x1++)
            {
                double sq_err = calculate_sq_error(x0, x1, y0, y1, nx, ny, total_avg, end_sum_squared, sum_from_zero, sum_squared_from_zero);

                if (sq_err < min_err)
                {
                    min_err = sq_err;
                    best_solutions[i] = sq_err;
                    best_solutions_coord[i * 4 + 0] = x0;
                    best_solutions_coord[i * 4 + 1] = x1;
                    best_solutions_coord[i * 4 + 2] = y0;
                    best_solutions_coord[i * 4 + 3] = y1;
                }
            }
        }
    }

    int min_i = 0;
    for (int i = 0; i < nx * ny; i++)
        if (best_solutions[i] < best_solutions[min_i])
            min_i = i;

    int x0 = best_solutions_coord[min_i * 4 + 0];
    int x1 = best_solutions_coord[min_i * 4 + 1];
    int y0 = best_solutions_coord[min_i * 4 + 2];
    int y1 = best_solutions_coord[min_i * 4 + 3];
    set_result(x0, x1, y0, y1, nx, ny, sum_from_zero, total_avg, result);

    return result;
}