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

function [pFB_hr,uFB_hr,indexDict_hr,pFB_lr,uFB_lr,indexDict_lr,ttList,bw] = genFB(rcode,code,fs,M,numFilters,minBaudRatio,maxBaudRatio,baseTime)
%GENFB3 Summary of this function goes here
%   Detailed explanation goes here [pFB,uFB,ldict,su]
    baseSam = baseTime*fs;
    rangeSam = linspace(minBaudRatio*baseSam,maxBaudRatio*baseSam,numFilters);
    ttList = rangeSam/fs;
    % generate high resolution filter banks
    maxLen = sum(ceil(max(rangeSam)*rcode))+1;
    pFB_hr = zeros(numFilters,maxLen);
    uFB_hr = zeros(numFilters,maxLen);
    bw = zeros(numFilters);
    for ii = 1:numFilters
        samps = rangeSam(ii)*rcode;
        utfilt = genFilter(code,samps);
        ptfilt = 2*utfilt - 1;
        bw(ii) = samps(1);
        
        % generate FB
        start = floor((size(pFB_hr,2) - length(ptfilt))/2) + 1;
        stop = start + length(ptfilt) - 1;
        pFB_hr(ii,start:stop) = ptfilt'/norm(ptfilt);
        uFB_hr(ii,start:stop) = utfilt';
        
        % generate index dictionary
        indexDict_hr(ii,:) = [length(utfilt); start; stop];
    end
    
    % generates low resolution filter bank
    pFB_lr = zeros(numFilters,ceil(max(rangeSam)/M)+1);
    uFB_lr = zeros(numFilters,ceil(max(rangeSam)/M)+1);
    for ii = 1:numFilters
%         pFB_lr(ii,:) = resample(pFB_hr(ii,:),1,M);
%         uFB_lr(ii,:) = resample(uFB_hr(ii,:),1,M);
        temp = interp1(1:length(uFB_hr(ii,:)),uFB_hr(ii,:),1:M:length(uFB_hr(ii,:)));
        uFB_lr(ii,1:length(temp)) = temp;%./norm(temp);
        temp = interp1(1:length(pFB_hr(ii,:)),pFB_hr(ii,:),1:M:length(pFB_hr(ii,:)));
        pFB_lr(ii,1:length(temp)) = temp;%./norm(temp);
        

        % generate index dictionary
        start = find(temp~=0,1,'first');
        stop = find(temp~=0,1,'last');
        indexDict_lr(ii,:) = [stop-start,start,stop];
    end
    bw = bw/M;
    bw = 1./bw;
end

