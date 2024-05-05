// #include <limits>
// TODO delete
#include <cstdio>
#include <vector>
using namespace std;

// constexpr double infty = numeric_limits<double>::infinity();
constexpr double infty = 1.7e308;
constexpr int C = 3;
// TODO delete
constexpr double epsilon = 1e-10;

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[C];
    float inner[C];
};

int id(int x, int y, int c, int nx)
{
    return c + C * x + C * nx * y;
}

double color(int x, int y, int c, int nx, const float *data)
{
    return data[id(x, y, c, nx)];
}

double error(int x, int y, int c, int nx, const float *data, double *area)
{
    return area[c] - color(x, y, c, nx, data);
}

bool is_inside(int x, int y, int x0, int x1, int y0, int y1)
{
    return x >= x0 && x < x1 && y >= y0 && y < y1;
}

double squared_error(int x0, int x1, int y0, int y1, int nx, int ny, const float *data, double *outer, double *inner)
{
    double sq_err = 0.0;
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            double *area;
            if (is_inside(x, y, x0, x1, y0, y1))
                area = inner;
            else
                area = outer;
            for (int c = 0; c < C; c++)
            {
                double err = error(x, y, c, nx, data, area);
                sq_err += err * err;
            }
        }
    }
    return sq_err;
}

void avg_colors(int x0, int x1, int y0, int y1, int nx, int ny, const float *data, double *outer, double *inner)
{
    for (int c = 0; c < C; c++)
        outer[c] = inner[c] = 0;

    int inner_points = (x1 - x0) * (y1 - y0);
    int outer_points = nx * ny - inner_points;
    // printf("inner points: %d * %d = %d\n", x1 - x0, y1 - y0, inner_points);

    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            bool inside = is_inside(x, y, x0, x1, y0, y1);
            // printf("(%d,%d) is inside: %d\n", x, y, inside);
            for (int c = 0; c < C; c++)
            {
                // printf("color %d: %f,\t", c, color(x, y, c, nx, data));
                if (inside)
                    inner[c] += color(x, y, c, nx, data) / inner_points;
                else
                    outer[c] += color(x, y, c, nx, data) / outer_points;
            }
            // printf("\n");
        }
    }
}

/** O(n^2)*/
void calculate_total_avg(int ny, int nx, const float *data, double *total_avg)
{
    int tot = ny * nx;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int c = 0; c < C; c++)
                total_avg[c] += color(x, y, c, nx, data) / tot;
}

/** O(n^2) */
void calculate_avg_from_zero(int ny, int nx, const float *data, vector<double> &avg_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int c = 0; c < C; c++)
            {
                if (y == 0)
                {
                    if (x == 0)
                    {
                        // * set first element
                        avg_from_zero[id(x, y, c, nx)] = color(x, y, c, nx, data);
                    }
                    else
                    {
                        // * set first row
                        double prev_left_block = avg_from_zero[id(x - 1, y, c, nx)] * ((y + 1) * x);
                        avg_from_zero[id(x, y, c, nx)] = (prev_left_block + color(x, y, c, nx, data)) / ((x + 1) * (y + 1));
                    }
                }
                else if (x == 0)
                {
                    // * set first column
                    double prev_up_block = avg_from_zero[id(x, y - 1, c, nx)] * (y * (x + 1));
                    avg_from_zero[id(x, y, c, nx)] = (prev_up_block + color(x, y, c, nx, data)) / ((x + 1) * (y + 1));
                }
                else
                {
                    double prev_left_block = avg_from_zero[id(x - 1, y, c, nx)] * ((y + 1) * x);
                    double prev_up_block = avg_from_zero[id(x, y - 1, c, nx)] * (y * (x + 1));
                    double prev_left_up_block = avg_from_zero[id(x - 1, y - 1, c, nx)] * (y * x);

                    avg_from_zero[id(x, y, c, nx)] = (prev_left_block + prev_up_block - prev_left_up_block + color(x, y, c, nx, data)) / ((x + 1) * (y + 1));
                }
            }
}

/** O(n^2) */
void calculate_sum_squared_from_zero(int ny, int nx, const float *data, vector<double> &sum_squared_from_zero)
{
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int c = 0; c < C; c++)
            {
                double squared_color = color(x, y, c, nx, data) * color(x, y, c, nx, data);
                if (y == 0)
                {
                    if (x == 0)
                    {
                        // * set first element
                        sum_squared_from_zero[id(x, y, c, nx)] = squared_color;
                    }
                    else
                    {
                        // * set first row
                        double prev_left_block = sum_squared_from_zero[id(x - 1, y, c, nx)];
                        sum_squared_from_zero[id(x, y, c, nx)] = (prev_left_block + squared_color);
                    }
                }
                else if (x == 0)
                {
                    // * set first column
                    double prev_up_block = sum_squared_from_zero[id(x, y - 1, c, nx)];
                    sum_squared_from_zero[id(x, y, c, nx)] = (prev_up_block + squared_color);
                }
                else
                {
                    double prev_left_block = sum_squared_from_zero[id(x - 1, y, c, nx)];
                    double prev_up_block = sum_squared_from_zero[id(x, y - 1, c, nx)];
                    double prev_left_up_block = sum_squared_from_zero[id(x - 1, y - 1, c, nx)];

                    sum_squared_from_zero[id(x, y, c, nx)] = (prev_left_block + prev_up_block - prev_left_up_block + squared_color);
                }
            }
}

int calculate_in_points(int x0, int x1, int y0, int y1)
{
    return (y1 - y0 + 1) * (x1 - x0 + 1);
}

/** O(1) */
void calculate_avg_in_color(int x0, int x1, int y0, int y1, int nx, vector<double> const &avg_from_zero, double *in)
{
    int in_points = calculate_in_points(x0, x1, y0, y1);
    for (int c = 0; c < C; c++)
    {
        if (y0 == 0)
        {
            if (x0 == 0)
            {
                in[c] = avg_from_zero[id(x1, y1, c, nx)];
            }
            else
            {
                double point_block = avg_from_zero[id(x1, y1, c, nx)] * ((y1 + 1) * (x1 + 1));
                double prev_left_block = avg_from_zero[id(x0 - 1, y1, c, nx)] * ((y1 + 1) * x0);
                in[c] = (point_block - prev_left_block) / in_points;
            }
        }
        else if (x0 == 0)
        {
            double point_block = avg_from_zero[id(x1, y1, c, nx)] * ((y1 + 1) * (x1 + 1));
            double prev_up_block = avg_from_zero[id(x1, y0 - 1, c, nx)] * (y0 * (x1 + 1));
            in[c] = (point_block - prev_up_block) / in_points;
        }
        else if (y1 == 0)
        {
            {
                // * set first row
                double x0y1_block = avg_from_zero[id(x0 - 1, y1, c, nx)] * x0;
                double x1y1_block = avg_from_zero[id(x1, y1, c, nx)] * (x1 + 1);
                in[c] = (x1y1_block - x0y1_block) / in_points;
            }
        }
        else if (x1 == 0)
        {
            // * set first column
            double x1y0_block = avg_from_zero[id(x1, y0 - 1, c, nx)] * y0;
            double x1y1_block = avg_from_zero[id(x1, y1, c, nx)] * (y1 + 1);
            in[c] = (x1y1_block - x1y0_block) / in_points;
        }
        else
        {
            double x0y0_block = avg_from_zero[id(x0 - 1, y0 - 1, c, nx)] * (y0 * x0);
            double x1y0_block = avg_from_zero[id(x1, y0 - 1, c, nx)] * (y0 * (x1 + 1));
            double x0y1_block = avg_from_zero[id(x0 - 1, y1, c, nx)] * ((y1 + 1) * x0);
            double x1y1_block = avg_from_zero[id(x1, y1, c, nx)] * ((y1 + 1) * (x1 + 1));
            in[c] = (x0y0_block + x1y1_block - x1y0_block - x0y1_block) / in_points;
        }
    }
}

/** O(1) */
void calculate_in_squared_sum(int x0, int x1, int y0, int y1, int nx, vector<double> const &sum_squared_from_zero, double *in)
{
    for (int c = 0; c < C; c++)
    {
        if (y0 == 0)
        {
            if (x0 == 0)
            {
                in[c] = sum_squared_from_zero[id(x1, y1, c, nx)];
            }
            else
            {
                double point_block = sum_squared_from_zero[id(x1, y1, c, nx)];
                double prev_left_block = sum_squared_from_zero[id(x0 - 1, y1, c, nx)];
                in[c] = point_block - prev_left_block;
            }
        }
        else if (x0 == 0)
        {
            double point_block = sum_squared_from_zero[id(x1, y1, c, nx)];
            double prev_up_block = sum_squared_from_zero[id(x1, y0 - 1, c, nx)];
            in[c] = point_block - prev_up_block;
        }
        else if (y1 == 0)
        {
            {
                // * set first row
                double x0y1_block = sum_squared_from_zero[id(x0 - 1, y1, c, nx)];
                double x1y1_block = sum_squared_from_zero[id(x1, y1, c, nx)];
                in[c] = x1y1_block - x0y1_block;
            }
        }
        else if (x1 == 0)
        {
            // * set first column
            double x1y0_block = sum_squared_from_zero[id(x1, y0 - 1, c, nx)];
            double x1y1_block = sum_squared_from_zero[id(x1, y1, c, nx)];
            in[c] = x1y1_block - x1y0_block;
        }
        else
        {
            double x0y0_block = sum_squared_from_zero[id(x0 - 1, y0 - 1, c, nx)];
            double x1y0_block = sum_squared_from_zero[id(x1, y0 - 1, c, nx)];
            double x0y1_block = sum_squared_from_zero[id(x0 - 1, y1, c, nx)];
            double x1y1_block = sum_squared_from_zero[id(x1, y1, c, nx)];
            in[c] = (x0y0_block + x1y1_block - x1y0_block - x0y1_block);
        }
    }
}

/** O(1) */
void calculate_avg_out_color(int x0, int x1, int y0, int y1, int nx, int ny, double *in, double *total_avg, double *out)
{
    int in_points = calculate_in_points(x0, x1, y0, y1);
    int full_points = (nx * ny);
    if (in_points == full_points)
    {
        return;
    }
    for (int c = 0; c < C; c++)
    {
        double in_block = in[c] * in_points;
        double full_block = total_avg[c] * full_points;
        out[c] = (full_block - in_block) / (full_points - in_points);
    }
}

void calculate_out_squared_sum(int nx, int ny, double *in, vector<double> const &sum_squared_from_zero, double *out)
{
    for (int c = 0; c < C; c++)
    {
        out[c] = sum_squared_from_zero[id(nx - 1, ny - 1, c, nx)] - in[c];
    }
}

double
calculate_in_error(double *in, double *in_squared_sum, int in_points)
{
    double sum0 = in_squared_sum[0] - in_points * in[0] * in[0];
    double sum1 = in_squared_sum[1] - in_points * in[1] * in[1];
    double sum2 = in_squared_sum[2] - in_points * in[2] * in[2];
    // printf("%f + %f + %f = %f\n", sum0, sum1, sum2, sum0 + sum1 + sum2);
    return sum0 + sum1 + sum2;
}

double calculate_out_error(double *out, double *out_squared_sum, int in_points, int full_points)
{
    int out_points = full_points - in_points;
    return calculate_in_error(out, out_squared_sum, out_points);
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
    double min_err = infty;
    double total_avg[C] = {0.0};
    calculate_total_avg(ny, nx, data, total_avg);
    vector<double> avg_from_zero(C * nx * ny, 0.0);
    calculate_avg_from_zero(ny, nx, data, avg_from_zero);
    vector<double> sum_squared_from_zero(C * nx * ny, 0.0);
    calculate_sum_squared_from_zero(ny, nx, data, sum_squared_from_zero);

    // TODO delete
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            printf("%f %f %f\n", avg_from_zero[id(x, y, 0, nx)], avg_from_zero[id(x, y, 1, nx)], avg_from_zero[id(x, y, 2, nx)]);
    printf("\n");
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            printf("%f %f %f\n", sum_squared_from_zero[id(x, y, 0, nx)], sum_squared_from_zero[id(x, y, 1, nx)], sum_squared_from_zero[id(x, y, 2, nx)]);
    printf("\n");

    for (int xdim = 1; xdim <= nx; xdim++)
    {
        for (int ydim = 1; ydim <= ny; ydim++)
        {
            for (int x0 = 0; x0 <= nx - xdim; x0++)
            {
                for (int y0 = 0; y0 <= ny - ydim; y0++)
                {
                    int x1 = x0 + xdim;
                    int y1 = y0 + ydim;
                    double outer[C] = {0.0};
                    double inner[C] = {0.0};

                    // printf("x0: %d, x1: %d, y0: %d, y1: %d\n", x0, x1, y0, y1);
                    avg_colors(x0, x1, y0, y1, nx, ny, data, outer, inner);
                    double sq_err = squared_error(x0, x1, y0, y1, nx, ny, data, outer, inner);
                    // printf("inner: %f,%f,%f\n", inner[0], inner[1], inner[2]);
                    // printf("outer: %f,%f,%f\n", outer[0], outer[1], outer[2]);
                    // printf("error = %f\n", sq_err);
                    // printf("\n");
                    if (sq_err < min_err)
                    {
                        min_err = sq_err;
                        result.x0 = x0;
                        result.x1 = x1;
                        result.y0 = y0;
                        result.y1 = y1;
                        for (int c = 0; c < C; c++)
                        {
                            result.inner[c] = inner[c];
                            result.outer[c] = outer[c];
                        }
                    }

                    // TODO delete
                    printf("x0: %d, x1: %d, y0: %d, y1: %d\n", x0, x1, y0, y1);

                    double outer2[C] = {0.0};
                    double inner2[C] = {0.0};
                    calculate_avg_in_color(x0, x1 - 1, y0, y1 - 1, nx, avg_from_zero, inner2);
                    calculate_avg_out_color(x0, x1 - 1, y0, y1 - 1, nx, ny, inner2, total_avg, outer2);

                    double outer_squared_sums[C] = {0.0};
                    double inner_squared_sums[C] = {0.0};
                    calculate_in_squared_sum(x0, x1 - 1, y0, y1 - 1, nx, sum_squared_from_zero, inner_squared_sums);
                    calculate_out_squared_sum(nx, ny, inner_squared_sums, sum_squared_from_zero, outer_squared_sums);
                    int in_points = calculate_in_points(x0, x1 - 1, y0, y1 - 1);
                    double sq_err2 = calculate_in_error(inner2, inner_squared_sums, in_points) + calculate_out_error(outer2, outer_squared_sums, in_points, nx * ny);

                    for (int c = 0; c < C; c++)
                        printf("%d: %f-%f = %f-%f\n", c, inner[c], inner2[c], outer[c], outer2[c]);
                    printf("sq_err: %f - %f\n", sq_err, sq_err2);
                    for (int c = 0; c < C; c++)
                    {
                        if ((inner2[c] - inner[c]) * (inner2[c] - inner[c]) > epsilon || (outer2[c] - outer[c]) * (outer2[c] - outer[c]) > epsilon)
                        {
                            printf("ERROR in %d\n", c);
                            Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
                            return result;
                        }
                    }
                    if ((sq_err - sq_err2) * (sq_err - sq_err2) > epsilon)
                    {
                        {
                            printf("ERROR in sq_err\n");
                            Result result{0, 0, 0, 0, {0, 0, 0}, {0, 0, 0}};
                            return result;
                        }
                    }

                    printf("ok\n\n");
                }
            }
        }
    }

    return result;
}