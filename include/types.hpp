#if !defined(NUMALALGOPDES_TYPES_H)
#define NUMALALGOPDES_TYPES_H

#include <Eigen/Core>
#include <complex>

namespace numal {
/*short hand for std double complex type */
using dcomplex = std::complex<double>;

/* double type Eigen matrix sized 3x1 */
using vec3 = Eigen::Matrix<double, 3, 1>;

/* double complex type Eigen matrix sized 3x1 */
using cvec3 = Eigen::Matrix<dcomplex, 3, 1>;

/* double type NdArray */
template <int dim> using Array = Eigen::Matrix<double, dim, 1>;

/* NdArray with date type defined by user */
template <typename T, int dim> using NdArray = Eigen::Matrix<T, dim, 1>;

} // namespace numal

#endif // NUMALALGOPDES_TYPES_H
