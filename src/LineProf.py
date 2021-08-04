from line_profiler import LineProfiler
import scattering_main_routine as barr
import integral_ijkr as ijk
lprofiler = LineProfiler()
lprofiler.add_function(ijk.integral_k_ijkr)
lp_wrapper = lprofiler(barr.main)
lp_wrapper()
lprofiler.print_stats()