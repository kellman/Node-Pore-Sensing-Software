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
% Runs the pulse height estimation analysis highlighted in Figure 5 and
% laid out in Section 6.2.
%% setup 
Geometry = 'MB13';
outfile = ['../../Results/' Geometry '_amplitude_stats.mat'];

transitTime = 0.2; % seconds
ph = 0.004;
SNR = 30;
sizeDevPer = 1/2;
edgeWidth = 1;
[output1,output2] = simulation(Geometry,SNR,ph,transitTime,sizeDevPer,edgeWidth,0);

%% simulating
% experimentally estimated noise level (average)
nstd =  1.2418e-04;
sizeDevPer = 1/2;
edgeWidth = 1;

Nsnr = 16;
SNR = linspace(0,31,Nsnr);
amps = sqrt((nstd.^2)*10.^(SNR/10));

Niter = 1000;
results = zeros(Nsnr,Niter,2);
for ii = 1:Nsnr
    for jj = 1:Niter
        disp(['SNR: ' num2str(SNR(ii)) '| iteration: ' num2str(jj)]);
        [sim,ideal,noisy_ideal] = simulation(Geometry,SNR(ii),amps(ii),transitTime,sizeDevPer,edgeWidth,plotting);

        % estimate pulse height
        ampls = ls(ideal,sim);
        amprr = rr(ideal,sim);
        
        
        cleanls = ls(ideal,noisy_ideal);
        cleanrr = rr(ideal,noisy_ideal);
        
        
        % record error
        results(ii,jj,2) = amps(ii) - ampls;
        results(ii,jj,3) = amps(ii) - amprr;
        results(ii,jj,4) = amps(ii) - cleanls;
        results(ii,jj,5) = amps(ii) - cleanrr;
    end
end
%% saving
save(outfile);