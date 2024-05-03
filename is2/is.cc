#include <limits>
#include <algorithm>
using namespace std;

struct Result
{
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[3];
    float inner[3];
};

constexpr double infty = numeric_limits<double>::infinity();
constexpr int C = 3;

double color(int x, int y, int c, int nx, const float *data)
{
    return data[c + C * x + C * nx * y];
}

double error(int x, int y, int c, int nx, const float *data, double *area)
{
    return area[c] - color(x, y, c, nx, data);
}

double squared_error(int x0, int x1, int y0, int y1, int nx, int ny, const float *data, double *outer, double *inner)
{
    double sq_err = 0.0;
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            double *area;
            if (x >= x0 && x <= x1 && y >= y0 && y <= y1)
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

    return result;
}
