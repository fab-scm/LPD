# Development of a local data point density approach to assess the area of applicability for spatial prediction models
This repository contains the code for the analysis in the thesis "Development of a local data point density approach to assess the area of applicability for spatial prediction models" by Fabian Schumacher.

Following contents can be found in the repository:

* [simulationStudy.R](analysis/simulationStudy.R): Contains the code to run the simulation study for all scenarios individually. Plots are generated step by step, but they do not match the ones in the thesis. Recommended to run, if results want to be reproduced and understood.
* [simulationStudy_withPlotting.R](analysis/simulationStudy_withPlotting.R): Contains the code to run the simulation study for all sampling scenarios simultaneously. At the same time the plots that can be found in the thesis are generated.
* [caseStudy.R](analysis/caseStudy.R): Contains the code to run the case study.
* [data](data/): Contains the predictor data used in the case study, the sampling locations for all studies, the model used in the case study and the modeldomain for the case study.

Please use the [LPD.RProj] to ensure proper relative file paths.