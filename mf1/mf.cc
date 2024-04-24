#include <algorithm>
#include <vector>

using namespace std;

float median(vector<float> &v)
{
  int n = v.size();
  if (n % 2 == 0)
  {
    int k = n / 2;
    nth_element(v.begin(), v.begin() + k, v.end());
    nth_element(v.begin(), v.begin() + k - 1, v.end());
    return (v[k] + v[k - 1]) / 2;
  }
  else
  {
    int k = (n - 1) / 2;
    nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
  }
}

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in in[x + y*nx]
- for each pixel (x, y), store the median of the pixels (a, b) which satisfy
  max(x-hx, 0) <= a < min(x+hx+1, nx), max(y-hy, 0) <= b < min(y+hy+1, ny)
  in out[x + y*nx].
*/
void mf(int ny, int nx, int hy, int hx, const float *in, float *out)
{
  for (int x = 0; x < nx; x++)
  {
    for (int y = 0; y < ny; y++)
    {
      vector<float> values = {};
      for (int i = max(0, x - hx); i < min(nx, x + hx + 1); i++)
      {
        for (int j = max(0, y - hy); j < min(ny, y + hy + 1); j++)
        {
          values.push_back(in[i + nx * j]);
        }
      }
      out[x + nx * y] = median(values);
    }
  }
}
