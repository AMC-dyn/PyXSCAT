from line_profiler import LineProfiler
import scattering_main_routine as barr
import integral_ijkr_vect as ijk
from integrals_wrapper import integrals_ijkr
lprofiler = LineProfiler()
lprofiler.add_function(integrals_ijkr.tot_integral_k_ijkr)
lp_wrapper = lprofiler(barr.main)
lp_wrapper()
lprofiler.print_stats()