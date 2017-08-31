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

function [nonideal, ideal, noisy] = simulation(Geometry,SNR,ph,transitTime,sizeDevPer,edgeWidth,plotting)
    %% setup to get geometry of channel
    if(strcmp(Geometry,'MB13'))
        setupB13;
    elseif(strcmp(Geometry,'MB11'))
        setupB11;
    elseif(strcmp(Geometry,'MB7'))
        setupB7;
    end
    
    %% simulation setup
    fs = 50000;
    N = fs/2;
    Npore = fs*transitTime/len;
    M = 15;
    
    if(edgeWidth == 0)
        edgeImpulse = 1;
    else
        edgeImpulse = hanning(edgeWidth*fs*transitTime/100);
        edgeImpulse = edgeImpulse./sum(edgeImpulse);
    end
    
    %% simulation
    sps = rcode'./rcode(1)*Npore;
    varp = 2*sizeDevPer*rand(length(rcode),1)-sizeDevPer;
%     varp = 1 + varp/100;
    varp = varp/100;
%     Nsam = sps'.*varp;
    Nsam = sps' + sum(sps)*varp;
%     length(Nsam)
    Nsam = round(Nsam);
    
    % non-ideal signal
    test = genFilter(code,Nsam);
%     length(test)
    test2 = zeros(N,1);
%     length(test2)

    start = length(test2)/2 - length(test)/2;
    stop = start + length(test) - 1;
    
    test2(start:stop) = test;
    test3 = conv(test2,edgeImpulse,'same');
    dstest1 = downsample(test3,M);
    
    % scale and add noise
    nstd = 1./sqrt(10.^(SNR/10)./(ph^2));
    n = nstd*randn(length(dstest1),1);
    nonideal = ph*dstest1+n;
%     size(nonideal)
    
    
    N = fs/2;
    % ideal signal
    ord = genFilter(code,sps);
%     length(ord)
    ord2 = zeros(N,1);
%     length(ord2)
    start = length(ord2)/2 - length(ord)/2;
    stop = start + length(ord) - 1;
    ord2(start:stop) = ord;
%     ord3 = conv(ord2,edgeImpulse,'same');
    ord4 = downsample(ord2,M);
    ideal = ord4; % height 1
    
    noisy = ph*ideal + n;
    
    if(plotting)
        figure, subplot 411, plot(nonideal), hold on, plot(ph.*ideal);
        subplot 412, plot(ph*ideal-nonideal);
        subplot 413, plot(noisy), hold on, plot(ph*ideal);
        subplot 414, plot(noisy - ph*ideal);
%         size(ideal)
%         size(ph)
    end
end