// The threading runtimes: the OpenMP loops take the thread count we set; MKL (PARDISO) runs on its TBB
// threading layer, limited by tbb::global_control.
#include <catch2/catch_test_macros.hpp>

#ifdef TURBOLENZA_HAVE_OPENMP
#include <omp.h>
#endif
#ifdef TURBOLENZA_HAVE_PARDISO
#include <mkl.h>
#include <oneapi/tbb/global_control.h>
#endif

#include <vector>

TEST_CASE("OpenMP takes the thread count; MKL runs on TBB", "[turbolenza][runtime]")
{
#ifdef TURBOLENZA_HAVE_OPENMP
    for (int n : {1, 2, 3}) {
        omp_set_num_threads(n);
        int seen = 0;
#pragma omp parallel
        {
#pragma omp single
            seen = omp_get_num_threads();
        }
        CHECK(seen == n);
#ifdef TURBOLENZA_HAVE_PARDISO
        const oneapi::tbb::global_control limit(oneapi::tbb::global_control::max_allowed_parallelism,
                                                static_cast<std::size_t>(n));
        CHECK(oneapi::tbb::global_control::active_value(oneapi::tbb::global_control::max_allowed_parallelism) ==
              static_cast<std::size_t>(n));
        MKLVersion v;
        mkl_get_version(&v);
        CHECK(v.MajorVersion >= 2024);
#endif
    }
#else
    SKIP("built without OpenMP");
#endif
}
