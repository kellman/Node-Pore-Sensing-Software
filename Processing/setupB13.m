% Copyright 2017 Michael Kellman, Michael Lustig
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% file paths
filepath = '../Data/';
filename = 'B13_colloidmix_newprotein_42116_1PSI';
infile = [filepath filename '.mat'];
outfile = ['../Results/' filename '_output.mat'];

% MB13 code structure
code = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
rcode = [1 1 1 1 1 1 1 1 1 2 1 1 2 1 1 2 2 2 2 1];
rcode = rcode'./sum(rcode);
len = 26;

% fitting parameters
epsilon = 1e-6; % reweighting epsilon
maxRWLSIter = 5; % iterations (5 is default) but 1 is Least Squares
lambda = 7.5e3; % tradeoff between data consistency and smooth baseline

% stopping criteria
maxSICIter = 6; % maximum number of iterations
thr = 0.9e-4; % minimum particle size

% detection parameters
beta = 7;
