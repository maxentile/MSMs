These are initial scripts to make MSMs of Src and Abl.

The bash and python scripts do the same thing as the ipynb, only the ipynb's also add ensembler models.

Data used is available on abl in Sonya's dir (eventually we should move this to choderalab dir):

/home/hansons/sims/src/10467/long_sims/
/home/hansons/sims/abl/10468/long_sims/

In order to use ipython notebooks remotely:
On abl:
`ipython notebook --no-browser --port=7000`

On your machine:
`ssh -N -L localhost:8878:localhost:7000 hansons@abl.mskcc.org&`

Then open http://localhost:8878/ in a browser.
