Copyright 2017 Michael Kellman, Michael Lustig

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This is the dataset for the work titled and authored by:

# Node-Pore Coded Coincidence Correcting Microfluidic Channel Framework: Code Design and Sparse Deconvolution

## Michael Kellman, Francois Rivest, Alina Pechacek, Lydia Sohn, Michael Lustig

We present a novel method to perform individual particle (e.g. cells or viruses) coincidence correction through joint channel design and algorithmic methods. Inspired by multiple-user communication theory, we modulate the channel response, with Node-Pore Sensing, to give each particle a binary Barker code signature. When processed with our modified successive interference cancellation method, this signature enables both the separation of coincidence particles and a high sensitivity to small particles. We identify several sources of modeling error and mitigate most effects using a data-driven self-calibration step and robust regression. Additionally, we provide simulation analysis to highlight our robustness, as well as our limitations, to these sources of stochastic system model error. Finally, we conduct experimental validation of our techniques using several encoded devices to screen a heterogeneous sample of several size particles.

## What is included:

The software is included here for demostrations, experiments with Barker 7, 11, and 13 encoded channels, and simulation for analysis purposes.

* Demo - demonstrates our algorithmic methods on experimental data (see below) for Barker 7, 11, 13 encoded channels
* Analysis - performs and displays pulse height, transit time, and false alarm analysis
* Processing - general processing to perform algorithmic methods on whole datasets
* Util - several useful/required helper functions

## How to use:

This software can be used with these datasets for reproducibility purposes at this DOI:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.860191.svg)](https://doi.org/10.5281/zenodo.860191)

## Citation:

This dataset can be cited with the following DOI:

10.5281/zenodo.846448

## Requirements:

cvx, free optimization software, is required for calibration and some of the simulations.
Link: http://cvxr.com/cvx/download/





