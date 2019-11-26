# MongooseGUI3
-----------------------------------------------------------------
MONGOOSE ( MetabOlic Network GrOwth Optimization Solved Exactly )
-----------------------------------------------------------------

MONGOOSE is a package for structural analysis and refinement of constraint-based metabolic networks.
Unlike other existing softwares, MONGOOSE uses exact rational arithmetic, which makes its results
certifiably accurate. The MetaMerge algorithm (Chindelevitch et al, Genome Biology 2012,
13:r6, http://genomebiology.com/2012/13/1/r6), is based on and fully integrated with MONGOOSE.
MONGOOSE's operation requires the esolver executable from QSOpt_ex, available to download
at http://www.dii.uchile.cl/~daespino/ESolver_doc/main.html to be located in the working directory
from which you run the Python script ModelProcessing.py. To check the results of an external
program's solution to a metabolic problem please use the script ExternalChecker.py.

The GUI can be run using Docker. Please see our Dockerhub: https://hub.docker.com/r/ctlevn/mongoose/

