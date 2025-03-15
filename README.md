athena with James Poisson solver
======

<p align="center">
<img src="https://github.com/user-attachments/assets/d1c7fd6a-39a8-4d73-8596-d8d76374abc5" width="600">
</p>

This version implements the James algorithm based on [Moon, Kim, & Ostriker (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJS..241...24M/abstract) to enable the self-gravity with open boundary conditions in three-dimensional Cartesian and **cylindrical** coordinates.

## Usage
Simply configure the code with `--grav james` option.
```
> ./configure.py --prob poisson --grav james -fft -mpi
> make clean
> make
```
Note that `-mpi` is necessary even when using a single core.

You need to write your own problem generator under the `pgen` directory. In your problem generator, make sure you set the gravitational constant by calling either `Mesh::SetGravitationalConstant` or `Mesh::SetFourPiG`. Examples can be found in, e.g., `pgen/poisson.cpp` and `pgen/disk.cpp`. 

## Restrictions
Because the current parallel FFT interface requires Mesh to be evenly divisible for both the block and pencil decompositions, the solver may not work for certain number of cells or MeshBlock decompositions. In addition, the James algorithm involves FFTs acting only on surfaces, which add complications on the possible decompositions. **The easiest way to meet all the requirements is to set the number of cells and MeshBlocks in each direction as powers of two.**

## Acknowledgement

Due to repository changes, the current repository does not contain full commit histories. I thank **Hans Baehr** and **Jeong Hyeon Ahn** for correcting the errors in the gravitational source term in cylindrical coordinates.

## Citation
To cite Athena++ in your publication, please use the following BibTeX to refer to the code's [method paper](https://ui.adsabs.harvard.edu/abs/2020ApJS..249....4S/abstract):
```
@article{Stone2020,
	doi = {10.3847/1538-4365/ab929b},
	url = {https://doi.org/10.3847%2F1538-4365%2Fab929b},
	year = 2020,
	month = jun,
	publisher = {American Astronomical Society},
	volume = {249},
	number = {1},
	pages = {4},
	author = {James M. Stone and Kengo Tomida and Christopher J. White and Kyle G. Felker},
	title = {The Athena$\mathplus$$\mathplus$ Adaptive Mesh Refinement Framework: Design and Magnetohydrodynamic Solvers},
	journal = {The Astrophysical Journal Supplement Series},
}
```
Additionally, you can add a reference to `https://github.com/PrincetonUniversity/athena` in a footnote.

Finally, we have minted DOIs for each released version of Athena++ on Zenodo. This practice encourages computational reproducibility, since you can specify exactly which version of the code was used to produce the results in your publication. `10.5281/zenodo.4455879` is the DOI which cites _all_ versions of the code; it will always resolve to the latest release. Click on the Zenodo badge above to get access to BibTeX, etc. info related to these DOIs, e.g.:

```
@software{athena,
  author       = {Athena++ development team},
  title        = {PrincetonUniversity/athena: Athena++ v24.0},
  month        = jun,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {24.0},
  doi          = {10.5281/zenodo.11660592},
  url          = {https://doi.org/10.5281/zenodo.11660592}
}
```

In addition to above, please cite our [method paper](https://ui.adsabs.harvard.edu/abs/2019ApJS..241...24M/abstract) for the cylindrical James algorithm  
```
@ARTICLE{Moon2019,
       author = {{Moon}, Sanghyuk and {Kim}, Woong-Tae and {Ostriker}, Eve C.},
        title = "{A Fast Poisson Solver of Second-order Accuracy for Isolated Systems in Three-dimensional Cartesian and Cylindrical Coordinates}",
      journal = {The Astrophysical Journal Supplement Series},
         year = 2019,
        month = apr,
       volume = {241},
       number = {2},
        pages = {24},
          doi = {10.3847/1538-4365/ab09e9},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJS..241...24M},
}
```


