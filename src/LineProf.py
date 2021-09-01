from line_profiler import LineProfiler
import scattering_main_routine as barr
import integral_ijkr_vect as ijk
from integrals_wrapper import integrals_ijkr
import molproinp_out as mp
import twordm as td
import molproinp_out as mp
lprofiler = LineProfiler()
lprofiler.add_function(mp.create_input)
lp_wrapper = lprofiler(barr.main)
lp_wrapper()
lprofiler.print_stats()