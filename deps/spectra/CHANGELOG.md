## [Unreleased]
### Added
- Added an `Uplo` template parameter to the `DenseSymShiftSolve` class
- Added the generalized eigen solver `SymGEigsSolver` in Cholesky mode
- Added the wrapper classes `DenseCholesky` and `SparseCholesky` that can be
  used in `SymGEigsSolver`
- Added test code for generalized eigen solver

### Changed
- Allowing basic math functions such as `abs()` and `sqrt()` to be overloaded
  (avoid using `std::abs` and `std::sqrt` directly), thanks to
  [@jdbancal](https://github.com/jdbancal). This makes it possible to use
  user-defined float number types with Spectra
- Replaced other `std` functions by their Eigen counterparts, for example using
  `Eigen::NumTraits<Scalar>::epsilon()` to substitute
  `std::numeric_limits<Scalar>::epsilon()`
- Fixed an out-of-bound [bug](https://github.com/yixuan/spectra/issues/14)
  detected by [@jdbancal](https://github.com/jdbancal)
- Improved numerical stability, e.g. the function `hypot(x, y)` is used to
  compute `sqrt(x^2 + y^2)`
- More careful in using "approximate zero" constants
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.5.7
- Improved documentation
- Updated Travis CI script


## [0.3.0] - 2016-07-03
### Added
- Added the wrapper classes `SparseSymMatProd` and `SparseSymShiftSolve`
  for sparse symmetric matrices
- Added the wrapper class `SparseGenRealShiftSolve` for general sparse matrices
- Added tests for sparse matrices
- Using Travis CI for automatic unit test

### Changed
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.5.6
- **API change**: Each eigen solver was moved to its own header file.
  For example to use `SymEigsShiftSolver` one needs to include
  `<SymEigsShiftSolver.h>`
- Header files for internal use were relocated


## [0.2.0] - 2016-02-28
### Added
- Benchmark script now outputs number of matrix operations
- Added this change log
- Added a simple built-in random number generator, so that the algorithm
  was made to be deterministic
- Added the wrapper class `DenseSymMatProd` for symmetric matrices

### Changed
- Improved Arnoldi factorization
  - Iteratively corrects orthogonality
  - Creates new residual vector when invariant subspace is found
  - Stability for matrices with repeated eigenvalues is greatly improved
- Adjusted deflation tolerance in double shift QR
- Updated result analyzer
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.3.4
- Updated copyright information
- **API change**: Default operator of `SymEigsSolver` was changed from
  `DenseGenMatProd` to `DenseSymMatProd`



## [0.1.0] - 2015-12-19
### Added
- Initial release of Spectra
