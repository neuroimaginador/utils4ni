#include <Rcpp.h>
using namespace Rcpp;

# define CUBE(x)   ((x) * (x) * (x))
# define SQR(x)    ((x) * (x))

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


int getIndex(int* size, int x, int y, int z) {

  int size_x = size[0], size_y = size[1], size_z = size[2];

  if (x >= size_x) x = size_x - 1;
  if (y >= size_y) y = size_y - 1;
  if (z >= size_z) z = size_z - 1;

  if (x < 0) x = 0;
  if (y < 0) y = 0;
  if (z < 0) z = 0;

  return z * size_x * size_y + y * size_x + x;

}

int get_index(int* size, int x, int y, int z) {

  int size_x = size[0], size_y = size[1], size_z = size[2];

  if (x >= size_x) return -1;
  if (y >= size_y) return -1;
  if (z >= size_z) return -1;

  if (x < 0) return -1;
  if (y < 0) return -1;
  if (z < 0) return -1;

  return z * size_x * size_y + y * size_x + x;

}

float evaluate_volume(double* vol, int* size, int x, int y, int z) {

  int index = get_index(size, x, y, z);
  float value = 0;

  if (index >= 0) value = vol[index];

  return value;

}

// https://github.com/Washington-University/gradunwarp/blob/master/gradunwarp/core/interp3_ext.c
float TriCubic (float px, float py, float pz, double* volume, int xDim, int yDim, int zDim)
{
  int             x, y, z;
  int    i, j, k;
  float           dx, dy, dz;
  // float *pv;
  int index;
  float           u[4], v[4], w[4];
  float           r[4], q[4];
  float           vox = 0;
  int             xyDim;

  xyDim = xDim * yDim;

  x = (int) px, y = (int) py, z = (int) pz;
  // necessary evil truncating at dim-2 because tricubic needs 2 more values
  // which is criminal near edges
  // future work includes doing trilinear for edge cases
  // range checking is extremely important here
  if (x < 2 || x > xDim-3 || y < 2 || y > yDim-3 || z < 2 || z > zDim-3)
    return (0);

  dx = px - (float) x, dy = py - (float) y, dz = pz - (float) z;
  index = (x - 1) + (y - 1) * xDim + (z - 1) * xyDim;

  /* factors for Catmull-Rom interpolation */

  u[0] = -0.5 * CUBE (dx) + SQR (dx) - 0.5 * dx;
  u[1] = 1.5 * CUBE (dx) - 2.5 * SQR (dx) + 1;
  u[2] = -1.5 * CUBE (dx) + 2 * SQR (dx) + 0.5 * dx;
  u[3] = 0.5 * CUBE (dx) - 0.5 * SQR (dx);

  v[0] = -0.5 * CUBE (dy) + SQR (dy) - 0.5 * dy;
  v[1] = 1.5 * CUBE (dy) - 2.5 * SQR (dy) + 1;
  v[2] = -1.5 * CUBE (dy) + 2 * SQR (dy) + 0.5 * dy;
  v[3] = 0.5 * CUBE (dy) - 0.5 * SQR (dy);

  w[0] = -0.5 * CUBE (dz) + SQR (dz) - 0.5 * dz;
  w[1] = 1.5 * CUBE (dz) - 2.5 * SQR (dz) + 1;
  w[2] = -1.5 * CUBE (dz) + 2 * SQR (dz) + 0.5 * dz;
  w[3] = 0.5 * CUBE (dz) - 0.5 * SQR (dz);

  for (k = 0; k < 4; k++)
  {
    q[k] = 0;
    for (j = 0; j < 4; j++)
    {
      r[j] = 0;
      for (i = 0; i < 4; i++)
      {
        r[j] += u[i] * volume[index];
        index++;
      }
      q[k] += v[j] * r[j];
      index += xDim - 4;
    }
    vox += w[k] * q[k];
    index += xyDim - 4 * xDim;
  }
  return vox;
}

//Trilinear interpolation kernel, wrapped from FLIRT
float trilinear_interpolation(float v000, float v001, float v010,
                              float v011, float v100, float v101,
                              float v110, float v111,
                              float dx, float dy, float dz) {

  float temp1, temp2, temp3, temp4, temp5, temp6;

  temp1 = (v100 - v000) * dx + v000;
  temp2 = (v101 - v001) * dx + v001;
  temp3 = (v110 - v010) * dx + v010;
  temp4 = (v111 - v011) * dx + v011;

  // second order terms
  temp5 = (temp3 - temp1) * dy + temp1;
  temp6 = (temp4 - temp2) * dy + temp2;

  // final third order term
  float result = (temp6 - temp5) * dz + temp5;

  return result;

}

float nearest_interpolation(float v000, float v001, float v010,
                            float v011, float v100, float v101,
                            float v110, float v111,
                            float dx, float dy, float dz) {

  float result = 0.0;

  if ((dx <= 0.5) & (dy <= 0.5) & (dz <= 0.5)) result = v000;
  if ((dx <= 0.5) & (dy <= 0.5) & (dz >= 0.5)) result = v001;
  if ((dx <= 0.5) & (dy >= 0.5) & (dz <= 0.5)) result = v010;
  if ((dx <= 0.5) & (dy >= 0.5) & (dz >= 0.5)) result = v011;

  if ((dx >= 0.5) & (dy <= 0.5) & (dz <= 0.5)) result = v100;
  if ((dx >= 0.5) & (dy <= 0.5) & (dz >= 0.5)) result = v101;
  if ((dx >= 0.5) & (dy >= 0.5) & (dz <= 0.5)) result = v110;
  if ((dx >= 0.5) & (dy >= 0.5) & (dz >= 0.5)) result = v111;

  return result;

}


//Get the interpolated value, wrapped from FLIRT
float interpolate(double* vol, int* vol_size,
                  float point_x, float point_y, float point_z,
                  int method) {

  if (method == 3) {

    float result = TriCubic (point_x, point_y, point_z, vol, vol_size[0], vol_size[1], vol_size[2]);

    return result;

  }

  int ix, iy, iz;
  ix = (int) floor(point_x);
  iy = (int) floor(point_y);
  iz = (int) floor(point_z);

  //left-top-front
  float dx = point_x - ix, dy = point_y - iy, dz = point_z - iz;

  float v000, v001, v010, v011, v100, v101, v110, v111;

  v000 = evaluate_volume(vol, vol_size, ix, iy, iz);
  v001 = evaluate_volume(vol, vol_size, ix, iy, iz + 1);
  v010 = evaluate_volume(vol, vol_size, ix, iy + 1, iz);
  v011 = evaluate_volume(vol, vol_size, ix, iy + 1, iz + 1);
  v100 = evaluate_volume(vol, vol_size, ix + 1, iy, iz);
  v101 = evaluate_volume(vol, vol_size, ix + 1, iy, iz + 1);
  v110 = evaluate_volume(vol, vol_size, ix + 1, iy + 1, iz);
  v111 = evaluate_volume(vol, vol_size, ix + 1, iy + 1, iz + 1);

  float result;
  if (method == 1) {

    //Use trilinear interpolation
    result = trilinear_interpolation(v000, v001, v010, v011, v100, v101, v110, v111, dx, dy, dz);

  } else {

    result = nearest_interpolation(v000, v001, v010, v011, v100, v101, v110, v111, dx, dy, dz);

  }

  return result;

}

void transform_volume(double* source_vol, double* target_vol, int* source_dims, int* target_dims, double* matrix, int method) {

  for (int x = 0; x < target_dims[0]; x++) {

    for (int y = 0; y < target_dims[1]; y++) {

      for (int z = 0; z < target_dims[2]; z++) {

        float value;

        int index = getIndex(target_dims, x, y, z);

        float px = x - target_dims[0] / 2;
        float py = y - target_dims[1] / 2;
        float pz = z - target_dims[2] / 2;

        float source_x = matrix[0] * px + matrix[4] * py + matrix[8] * pz + matrix[12];
        float source_y = matrix[1] * px + matrix[5] * py + matrix[9] * pz + matrix[13];
        float source_z = matrix[2] * px + matrix[6] * py + matrix[10] * pz + matrix[14];

        source_x += source_dims[0] / 2;
        source_y += source_dims[1] / 2;
        source_z += source_dims[2] / 2;

        value = interpolate(source_vol, source_dims, source_x, source_y, source_z, method);

        target_vol[index] = value;

      }

    }

  }

}

// [[Rcpp::export]]
NumericVector transform_volume(NumericVector V, NumericVector M, IntegerVector target_dims, int method) {

  IntegerVector source_dims = V.attr("dim");

  NumericVector transformed_volume(target_dims[0] * target_dims[1] * target_dims[2]);

  transform_volume(V.begin(),
                   transformed_volume.begin(),
                   source_dims.begin(),
                   target_dims.begin(),
                   M.begin(),
                   method);

  transformed_volume.attr("dim") = target_dims;

  return(transformed_volume);

}


void deform_volume(double* source_vol, double* target_vol, int* source_dims, int* target_dims, double* dx, double* dy, double* dz, int method) {

  for (int x = 0; x < target_dims[0]; x++) {

    for (int y = 0; y < target_dims[1]; y++) {

      for (int z = 0; z < target_dims[2]; z++) {

        float value;

        int index = getIndex(target_dims, x, y, z);

        float source_x = x + dx[index];
        float source_y = y + dy[index];
        float source_z = z + dz[index];

        value = interpolate(source_vol, source_dims, source_x, source_y, source_z, method);

        target_vol[index] = value;

      }

    }

  }

}

void deform_volume(double* source_vol, double* target_vol, int* source_dims, int* target_dims, int* tx, int* ty, int* tz, int method) {

  for (int x = 0; x < target_dims[0]; x++) {

    for (int y = 0; y < target_dims[1]; y++) {

      for (int z = 0; z < target_dims[2]; z++) {

        float value;

        int index = getIndex(target_dims, x, y, z);

        float source_x = x + tx[index];
        float source_y = y + ty[index];
        float source_z = z + tz[index];

        value = interpolate(source_vol, source_dims, source_x, source_y, source_z, method);

        target_vol[index] = value;

      }

    }

  }

}

// [[Rcpp::export]]
NumericVector deform_volume(NumericVector V, NumericVector Dx, NumericVector Dy, NumericVector Dz, IntegerVector target_dims, int method) {

  IntegerVector source_dims = V.attr("dim");

  NumericVector transformed_volume(target_dims[0] * target_dims[1] * target_dims[2]);

  deform_volume(V.begin(),
                   transformed_volume.begin(),
                   source_dims.begin(),
                   target_dims.begin(),
                   Dx.begin(),
                   Dy.begin(),
                   Dz.begin(),
                   method);

  transformed_volume.attr("dim") = target_dims;

  return(transformed_volume);

}

// [[Rcpp::export]]
NumericVector deform_volume_candidates(NumericVector V, IntegerVector Dx, IntegerVector Dy, IntegerVector Dz, IntegerVector target_dims, int method) {

  IntegerVector source_dims = V.attr("dim");

  NumericVector transformed_volume(target_dims[0] * target_dims[1] * target_dims[2]);

  deform_volume(V.begin(),
                transformed_volume.begin(),
                source_dims.begin(),
                target_dims.begin(),
                Dx.begin(),
                Dy.begin(),
                Dz.begin(),
                method);

  transformed_volume.attr("dim") = target_dims;

  return(transformed_volume);

}
