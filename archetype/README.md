# archetype

Fast prototyping in architecture can be of immense help to engineers and designers to understand the solution space. This project aims to show how fast Quality Diversity and Lattice-Boltzmann techniques can be integrated into an architectural design process. The resulting morphological design classes are evolved towards the criterion of minimizing wind discomfort around a building.

# Computer Aided Ideation: Prototype Discovery using Quality-Diversity

This software shows how quality diversity algorithms can be used for (interactive) computer aided ideation. 

It includes: 
- Matlab demo (demo.m)
- Demo GUI (demoGUItabs.mlapp) created with Matlab Appdesigner (2019a)
- Standalone application (needs Matlab Runtime, see demoGUItabs/for_redistribution)

Please include the following references in any publication using this code. For Bibtex please see the end of this file.

Hagg, A., Asteroth, A. and Bäck, T., 2018, September. Prototype discovery using quality-diversity. In International Conference on Parallel Problem Solving from Nature (pp. 500-511). Springer, Cham.
Hagg, A., Asteroth, A. and Bäck, T., 2019, July. Modeling user selection in quality diversity. In Proceedings of the Genetic and Evolutionary Computation Conference (pp. 116-124). ACM.


## Before you begin
- Code was written for Matlab 2019a
- Please compile /similarity/bhtsne. In order to do so, please follow the /classes/bhtsne/README.md


## Examples
- demoGUItabs.mlapp     -   the GUI as a Matlab application
- ... contains a standalone version of the GUI that you can use after installing the Matlab runtime
- demo.m                -   non-GUI code usage example


## Documentation
- More explanation about the function headers is provided in every file, which can be accessed using "help <functionname>"


## How to add a domain
Add a domain folder in domain/. Please use the same structure and naming conventions as in domain/npoly_ffd. The utils folder can contain some domain-specific functions. 

The following functions are expected:

categorize.m        -   phenotypic features that can be extracted from a phenotype in this domain
domain.m            -   parameterization of the domain
getPhenotype.m      -   create phenotype based on genome
showExample.m       -   visualize phenotype
validate.m          -   can contain validation of known (non-user-selection) validity constraints. 

## Code Structure
classification	    - 	Create similarity space with Barnes-Hut t-SNE and k-medoids clustering
			& extract prototypes (medoid, per class)
constraints	    - 	Create a model of user selection in context of similarity space
domain		    - 	Contains domain specifications
models		    - 	Create models and predict with models using GPML library
QD		    - 	Variations of quality diversity algorithms: MAP-Elites (grid) or voronoi-
 			Elites (voronoi)


## todo

1. QD polygon domain with symmetry objective
    1. ~~Simple representation of polygon nodes~~
    2. ~~Features: circumference, area~~
    3. ~~Objective: point symmetry~~
2. QD polygon domain with drag objective
    1. Integrate Lattice Boltzman solver
        1. Check necessity of surrogate model
    2. Features: area, point symmetry
3. QD polygon domain with wind comfort objective
    1. Define measuring points
        1. Fixed grid, remove points covered by building footprints
        2. Measure total wind discomfort in an approximation of the NEN8100 building norm
        3. Single fixed wind direction, or multiple?
    2. Standard environment (find building requirement)
    3. Measure wind flow at n points
4. QD polygon domain with inducing flow in situ
    1. Add buildings in environment
    2. Get data from NYC, Amsterdam Zuidas?
5. VAE representation 
    1. Use polygons from MAP-Elites run
    2. Polygon images
    3. Postprocessing of VAE output
        1. Binary output
        2. Binary operations
    4. Show classes
6. Train VAE on building footprints
    1. Data set: open data NL, NYC? openstreetmap or https://code.waag.org/buildings/


## Notes about NEN8100 and wind discomfort

Measurement: 
- 108 locations around the project at 1.75 meters. Points are usually positioned at pedestrian walkways or any point where a pedestrian might be walking. Points are classified as “fast walking” or “strawling”.
- 12 different wind angles. Determination of windspeedcoefficients c_v (= v_1.75/v_60)
- Critical v_1.75 = 5 and 15
- Total probability of threshold surpassing = sum of surpassing probabilites per angle

Main causes of bad wind climate:
- Wind circulating around building (case 1: 2D)
  - Corners
  - Building mass (area)
- Drop winds (case 2: 3D)
  - Building height
  - Corners

Conditions that are not possible to solve with additional simple measures:
- Building mass

## Bibtex entrees
@inproceedings{hagg2018prototype,
  title={Prototype discovery using quality-diversity},
  author={Hagg, Alexander and Asteroth, Alexander and B{\"a}ck, Thomas},
  booktitle={International Conference on Parallel Problem Solving from Nature},
  pages={500--511},
  year={2018},
  organization={Springer}
}

@inproceedings{hagg2019modeling,
  title={Modeling user selection in quality diversity},
  author={Hagg, Alexander and Asteroth, Alexander and B{\"a}ck, Thomas},
  booktitle={Proceedings of the Genetic and Evolutionary Computation Conference},
  pages={116--124},
  year={2019},
  organization={ACM}
}


Copyright 2019 Alexander Hagg

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
