# Genomic Architecture

This is based in principle by code for [Ashraf and Lawson 2020](https://www.biorxiv.org/content/10.1101/2020.08.17.254110v1) (biorxiv link) "Genetic drift from the out-of-Africa bottleneck leads to biased estimation of genetic architecture and selection".

Some improvements were made based on comments by the Stan community

The model can be applied in posteriors from BayesRR-RC under the assumption that the the S parameter is  independant to all other genetic parameters conditioned on the posterior Betas. 
