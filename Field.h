#pragma once
#include <vector>
#include <stdexcept>
#include "Params.h"
#include "Index.h"

// Multi-variable spatial field.
// Layout: data[v * (nx*ny*nz) + spatial_flat_index]
// Each variable's spatial data is contiguous — good for stencil sweeps.
//
// Access:
//   1D  field(v, i)
//   2D  field(v, i, j)
//   3D  field(v, i, j, k)

struct Field
{
    int dim  = 1;
    int nvar = 1;
    int nx   = 1;
    int ny   = 1;
    int nz   = 1;

    std::vector<double> data;

    Field() = default;

    explicit Field(const Params& p)
        : dim(p.dim), nvar(p.nvar), nx(p.nx), ny(p.ny), nz(p.nz)
    {
        data.assign(nvar * nx * ny * nz, 0.0);
    }

    int spatial_size() const { return nx * ny * nz; }
    int total_size()   const { return nvar * nx * ny * nz; }

    // --- 1D access ---
    double& operator()(int v, int i)
    {
        return data[v * nx + idx1d(i, nx)];
    }
    const double& operator()(int v, int i) const
    {
        return data[v * nx + idx1d(i, nx)];
    }

    // --- 2D access ---
    double& operator()(int v, int i, int j)
    {
        return data[v * nx * ny + idx2d(i, j, nx, ny)];
    }
    const double& operator()(int v, int i, int j) const
    {
        return data[v * nx * ny + idx2d(i, j, nx, ny)];
    }

    // --- 3D access ---
    double& operator()(int v, int i, int j, int k)
    {
        return data[v * nx * ny * nz + idx3d(i, j, k, nx, ny, nz)];
    }
    const double& operator()(int v, int i, int j, int k) const
    {
        return data[v * nx * ny * nz + idx3d(i, j, k, nx, ny, nz)];
    }

    // Fill one variable.
    void fill(int v, double value)
    {
        int base = v * spatial_size();
        std::fill(data.begin() + base, data.begin() + base + spatial_size(), value);
    }

    // Fill all variables.
    void fill(double value)
    {
        std::fill(data.begin(), data.end(), value);
    }
};
