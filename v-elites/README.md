# Computer Aided Ideation: Prototype Discovery using Quality-Diversity

## Literature references

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