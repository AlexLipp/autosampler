This repository contains some tools and algorithms which streamline the analysis of point observations (e.g., sample sites, gauging stations, dams etc) alongside drainage basins. There are two python scripts: 

1. `auto_cathments.py` contains an algorithm (`AutoCatchments`) which automatically subdivides a drainage basin into sub-basins of approximately equal area, greater than some given threshold area. 
2. `process_samples.py` contains tools which streamline the analysis of a given set of river sample/observation sites. Specifically, it provides functions for 1) [aligning samples to drainage](https://onlinelibrary.wiley.com/doi/abs/10.1029/2007WR006507), and 2) subdividing a drainage basin into sub-catchments based on provided sample sites.  

### Cite

If you use this please cite it as:
`Lipp, Alex. 2022. AutoCatchments. Software. doi:  10.5281/zenodo.7311352`. 

Please also acknowledge: `Barnhart et al. 2020. doi: 10.5194/esurf-8-379-2020` and `Barnes et al. 2014. doi: 10.1016/j.cageo.2013.04.024`.

### Prerequisites 

These `python` scripts require the `LandLab` surface process modelling package. `LandLab` can be easily [installed](https://landlab.github.io/) using most package management software (e.g., `conda`, `pip`). It also requires standard scientific computing packages.


## 1. AutoCatchments


![Example output from AutoCatchments](https://user-images.githubusercontent.com/10188895/201136932-9fe7db1e-4d4c-4672-ad07-73aca1e5dd23.png)

`AutoCatchments` is an algorithm that, given a DEM and a target sub-catchment area, identifies N sample sites on a drainage network that sub-divide it into N sub-catchments which have an area greater than the given target area. This algorithm could be useful for designing sample-campaigns we ensure equal _area_ coverage across drainage basins for [geochemical exploration](https://doi.org/10.1016/0375-6742(87)90081-1), and [ecological/environmental monitoring](https://www.biorxiv.org/content/10.1101/2022.01.25.475970v1.abstract). 

This algorithm is different to previous catchment delineation approaches as it aims to produce sub-catchments which are approximately equal in size, and no smaller than a specified target area. Automatic catchment delineation as it is generally implemented (e.g., in [ArcGIS](https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-watershed-works.htm) or [QGIS](https://docs.qgis.org/2.8/en/docs/training_manual/processing/hydro.html)) typically works by subdividing a *stream* network into catchments. This can result in sub-basins with a wide degree of different _areas_, depending on the topology of the network. Additionally, they only allow for sample sites to be present at stream junctions. The new algorithm presented here automatically identifies the optimal sampling strategy, optimising for equal coverage, irrespective of stream junctions.

I make no claim that this is novel or even the best way of solving the problem... but hopefully it is useful. Boris Gailleton helped work out the most efficient approach.

### Usage 

The algorithm is designed to be run from the command line. It takes as argument a DEM (in `.asc`) format, and a target area (in the units of `.asc` file):

`python auto_catchments.py [path_to_dem.asc] [target_catchment_area]`

The locations of the sample sites (in units of **model grid cells**) are saved to `sample_sites.csv` and a raster map of the delineated catchments is saved to `area_IDs.asc`. 

#### Example

An example DEM from [North East Scotland, UK ](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GC009838) is provided for testing purposes. The above image shows the output from subdividing this DEM into sub-catchments greater than 500 km2 in size, using the command: 

`python auto_catchments.py cairngorms_topo.asc 5e8` 


### Algorithm Description 

I briefly detail how the algorithm works here but the code is well commented. I don't know if this has been done before (a very brief search found nothing) or whether there is a faster way of doing it (there probably is). 

1. Fill sinks and route drainage across DEM (using `LandLab` functions). 
2. Create a list of 'unvisited' nodes in topological order from downstream to upstream.
3. Identify sink nodes which have area smaller than the target area.
4. Remove these sinks and all their upstream nodes from the unvisited nodes list as these cannot be sampled.
5. Loop through all nodes from upstream to downstream 
6. For each node, calculate the number of unvisited nodes (and thus area) uniquely upstream of that point. 
7. If this area is greater than the threshold we add it to the list of sample sites and remove its upstream basin from the list of unvisited nodes.
8. Repeat until all nodes have been visited. 


For moderately sized DEMs (<100,000) cells the runtime is nearly instantaneous. For a DEM with 1 million cells it took around 20 minutes to complete on a standard desktop computer. Information on the progress of the algorithm is produced to the terminal. 

## 2. Sample Processing 

`process_samples.py` contains (well documented) tools which streamline the analysis of a given set of river sample/observation sites. Specifically, it provides functions for 1) [aligning samples to drainage](https://onlinelibrary.wiley.com/doi/abs/10.1029/2007WR006507), and 2) subdividing a drainage basin into sub-catchments based on provided sample sites.  

### Usage 

If used from command line this script is designed to take, 1) a DEM as a `.asc` 2) coordinates (in the same units as the DEM) of sample sites and 3) A drainage area threshold, which indicates the size of the channels which sample sites should be aligned to. 

Alternatively, the functions can be imported e.g., `import process_samples.py`. Each function is well documented.

### Example 

`python process_samples.py cairngorms_topo.asc noisy_sample_sites.csv 5e8` will align the sample site localities given in `noisy_sample_sites.csv` to channels greater than 500 km2, and then subdivide the DEM into the sub-catchments for each sample-site. It outputs `fitted_localities.csv` which gives the coordinates of the fitted localities, and the corresponding 'Area ID'. This ID number corresponds to a unique sub-catchment in the raster map output `area_IDs.asc`. 
