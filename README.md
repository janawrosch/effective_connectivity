# effective_connectivity
This repository includes all MATLAB code used to generate the data published in Wrosch et al., Rewiring of neuronal networks during synaptic silencing, 2017, doi: 10.1038/s41598-017-11729-5.

Please use the scripts in the following order:

A Fluorescence trace extractation 
  Input: Stacks of multilayered tif images
  Output: Fluorescence traces for each detected region of interest
B Spike estimation
  Input: Fluorescence traces for each detected region of interest
  Output: Binary spike time traces for each region of interest
C Network reconstruction
  Input: Binary spike time traces for each region of interest
  Output: Reconstructed effective network
D Network topology analysis
  Input: Reconstructed effective network
  Output: Network topology parameter
D2 Network activity analysis
  Input: Binary spike time traces for each region of interest
  Output: Network activity parameter
E Compare results of different experimental conditions
  Input: Reconstructed effective network, Binary spike time traces for each region of interest

Best
Jana Wrosch and colleagues
