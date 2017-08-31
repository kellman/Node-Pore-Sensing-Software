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
clear all; close all; clc; warning off;
addpath '../';
addpath '../../Util';
addpath '../../Processing';
plotting = 0;
%% Summary
% Runs the transit time estimation analysis highlighted in Figure 6 and
% laid out in Section 6.3.
%% setup
Geometry = 'MB11';

if(strcmp(Geometry,'MB13'))
    setupB13;
elseif(strcmp(Geometry,'MB11'))
    setupB11;
elseif(strcmp(Geometry,'MB7'))
    setupB7;
end
    
transitTime = 0.15; % seconds
ph = 0.004;
SNR = 30;
sizeDevPer = 1/2;
edgeWidth = 1;
fs = 50000;
M = 15;
% [output1,output2] = simulation(Geometry,SNR,ph,transitTime,sizeDevPer,edgeWidth,0);

% setup filterbank
minRatio = 0.2;
maxRatio = 1.8;
numFilters = 500;
[~,~,~,pFB_lr,uFB_lr,~,ttList] = genFB(rcode,code,fs,M,numFilters,minRatio,maxRatio,transitTime);
figure, imagesc(pFB_lr);

% setup MC simulation
Niter = 1000;
Ntimes = 3;
times = [0.1875 0.15 0.1125]; % pm25% average and average
results = zeros(Ntimes, Niter);

%% simulation
for ii = 1:Ntimes
    for jj = 1:Niter
        disp(['timing: ' num2str(times(ii)) '| Iteration: ' num2str(jj)]);
        [sim,ideal] = simulation(Geometry,SNR,ph,times(ii),sizeDevPer,edgeWidth,plotting);

        R = applyFB(sim,pFB_lr);
        
        [~,mi] = max(R(:));
        [mi,mj] = ind2sub(size(R),mi);

        results(ii,jj) = ttList(mi) - times(ii); 
    end
end
%% saving
outfile = ['../../Results/' Geometry '_transitTime_stats.mat'];
save(outfile);