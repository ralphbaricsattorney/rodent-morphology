import cleaning
from cleaning import *
from figures import make_effect_variability_scatter_plots


species_cleaners = [x for x in dir() if x.startswith("process")]
for cleaner in species_cleaners:
    e, v = getattr(cleaning, cleaner)()
    make_effect_variability_scatter_plots(e, v, "_".join(cleaner.split("_")[1:-1]))