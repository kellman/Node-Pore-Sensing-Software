% coded-NPS-software 
% copyright 2019 Kristen Cotner, Alan Dong

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

clear all; clc; close all; warning off;
addpath '../Util/'
addpath '../Data/'
%% Setup Demo
% loading experimental specs
% demo_setupB13;
demo_setupB11;
% demo_setupB7;

% load data
load(infile);
raw = 1./data';
global_start = 1;
global_stop = length(raw)-1;

% downsampling rate
M = 15; % from 50 kHz to ~3.333 kHz

% filterbank spec
baseTime = 0.2; % seconds
minRatio = 0.05;
maxRatio = 1.5;
numFilters = 100;
transitTime = baseTime*linspace(minRatio,maxRatio,numFilters);

% processing setup
% blockSize = 2500;
blockSize = Nstop - Nstart + 1;
blockOverlap = blockSize/2;

% baseline subtraction
bl_s = 1e6;
bl_p = 0.01;

% option flags
calib = 0; % perform system model calibration
plotting = 1; % plot while processing

%% Channel Calibration & Generate Matched Filter Bank
% (optional calibration)
hdata = double(1./raw);

if(calib)
    disp('Calibrating System Model');
    % perform system model calibration
    [~,~,~,pFB_lr_org,uFB_lr_org,indexDict_lr_org] = genFB(rcode,code,fs,M,numFilters,minRatio,maxRatio,baseTime);
    start_calib = 1;
    stop_calib = length(hdata) - 1;
    [rcode,~,~] = calibFB(hdata(start_calib:stop_calib),M,pFB_lr_org,10000,bl_s,bl_p,Geometry,0);
    rcode = [rcode; mean(rcode(1:5))];
    rcode = rcode./sum(rcode);
end

% generate matched filterbank (ideal/calibrated)
[pFB_hr,uFB_hr,indexDict_hr,pFB_lr,uFB_lr,indexDict_lr] = genFB(rcode,code,fs,M,numFilters,minRatio,maxRatio,baseTime);

if(plotting)
    figure; subplot 121; imagesc(pFB_lr); title('Matched Filterbank');
    subplot 122; imagesc(uFB_lr); title('System Model');
end
%% Downsample Experimental Data
hdata = hdata(global_start:global_stop); % offset to point of interest
ldata = resample(hdata,1,M); % decimate by M
ldata = ldata(50:end-50); % gets rid of transient resampling effects
N = floor(length(ldata)/blockOverlap); % compute number of frames to process

%% Successive Interference Cancellation Setup
% set beginning offset
offset = ceil(indexDict_lr(end,1)/2);

% 2nd Difference matrix
Ibig = eye(blockSize);
for zz = 1:length(Ibig)
    Ibig(zz,:) = conv(Ibig(zz,:),[1 -2 1],'same');
end; Ibig(1,1) = -1; Ibig(end,end) = -1;

% detection list initialization
isaved = [];
detectList = [];

% start index to not go out of bounds
istart = offset+1;

%% Demo Process Begin
gstart = Nstart;
gstop = Nstop;
offset = 1;
blockSize = Nstop - Nstart + 1;
tt = linspace(gstart/(fs/M),gstop/(fs/M),blockSize);
Frame = ldata(gstart:gstop);
bLine = ASLS(Frame,bl_s,bl_p,10);
dFrame = Frame-bLine;
R = [];
detectSICList = []; % detections for frame
%% Successive Iterference Cancellation
for gg = 1:maxSICIter
    disp(['SIC Iteration ' num2str(gg)]);
    %% Applying Matched Filterbank (applying adjoint)
    R(:,:,gg) = applyFB(dFrame,pFB_lr);
    if(plotting)
        f1 = figure; subplot 511; plot(tt,dFrame); axis([min(tt) max(tt) min(dFrame) max(dFrame)]);
        subplot 512; imagesc(tt,transitTime,R(:,:,gg)); axis xy; drawnow();
    end

    %% Adaptive Detection Step
    gdim1 = 5; idim1 = gdim1;
    gdim2 = 2*floor(baseTime*fs/M/len); idim2 = gdim2/2;
    bdim1 = 5*gdim1+idim1;
    bdim2 = gdim2+idim2;
    [detections,PSLR_2d] = adaptiveDetection(R(:,:,gg),indexDict_lr,beta,len);
    %% Selecting candidate detections
    candidateList = []; values = []; values2 = [];

    if(~isempty(detections))
        values = (detections(:,2) > offset) & (detections(:,2) < (blockSize + offset)) & (detections(:,1) > bdim1) & (detections(:,1) <= (numFilters-bdim1));
        detections = detections(values,:);

        p_values = diag(PSLR_2d(detections(:,1),detections(:,2)));
        corr_values = diag(R(detections(:,1),detections(:,2),gg));

        for pp = 1:size(detections,1)
            candidate = struct('time',detections(pp,2),'tti',detections(pp,1),'p_value',p_values(pp),'corr_value',corr_values(pp),'pol',1);
            candidateList = [candidateList candidate];
        end
    end

    if(isempty(detections))
        % stop condition: if there are no more detections
        disp('No more detections'); break;
    end

    if(plotting)
        figure(f1);
        subplot 513;
        tempval = getAll(candidateList,'corr_value')';
        cscale = 1-repmat(tempval,1,3)./max(tempval);
        mergedDetections = [getAll(candidateList,'tti'); getAll(candidateList,'time')]';
        scatter(tt(mergedDetections(:,2)-offset),transitTime(mergedDetections(:,1)),6,cscale);
        axis([min(tt) max(tt) min(transitTime) max(transitTime)]); hold on;
        hold off; drawnow(); shg;
    end
    %% Select Best Detection
    % maximum peak in frame that is not already in list
    [y,si] = sort(abs(getAll(candidateList,'corr_value'))); % sort in order of correlation
    cval = abs(getAll(candidateList,'corr_value')); % sort in order of correlation
    cval = cval(si);
    ss = length(si); % initialize from highest correlation detection
    while(ss > 0)
        currentDetection = struct('pos',candidateList(si(ss)).time,'ttindex',candidateList(si(ss)).tti, ...
            'pol',candidateList(si(ss)).pol,'amp',[],'correng',cval(ss),'pslr',y(ss),'blamp',[], ...
            'spfilt',[],'estBaud',[],'ttest',[],'index',NaN,'sig',[]);
            % removed reference to nonexistant variable `ii`
        if(~isempty(detectSICList))
        if(~isempty(detectSICList))
            gdim = indexDict_lr(currentDetection.ttindex,1)/len; % already in terms of decimated samples
            ddist2 = abs(currentDetection.pos - getAll(detectSICList,'pos')) < 2*gdim;
            if(all(0 == ddist2))
                % add to current frames detection list
                detectSICList = [detectSICList; currentDetection];
                disp('Adding detection to frame list'); break;
            end
        else
            % initialize the currect frame's detection list
            detectSICList = [detectSICList; currentDetection];
            disp('Initalizing frame detection list'); break;
        end
        ss = ss - 1;
    end
    if(ss == 0)
        % Stopping Conditions: couldn't find any unique detections
        disp('Cannot find anymore detections'); break;
    end
    %% Signal Model Regression Step (joint baseline and pulse height estimation)
    % setup pruned system model matrix
    A_RWLS = [];
    for ll = 1:size(detectSICList,1)
        CD = detectSICList(ll);
        TID = indexDict_lr;
        TFB = uFB_lr;

        % get signal to fit
        tempFilt = TFB(CD.ttindex,TID(CD.ttindex,2)+1:TID(CD.ttindex,3));
        dlength = length(tempFilt);
        start = floor(CD.pos) - round(dlength/2);
        stop = start + dlength - 1;
        zpbuffer = zeros(size(Frame));

        % support matrix
        zpbuffer(start:stop) = tempFilt;
        detectSICList(ll).spfilt = zpbuffer;
        detectSICList(ll).estBaud = transitTime(detectSICList(ll).ttindex);
        A_RWLS = [A_RWLS zpbuffer];

        if(plotting)
            % plot current frame's detections
            figure(f1);
            subplot 513; hold on; plot(tt(CD.pos-offset),transitTime(CD.ttindex),'g*'); hold off; drawnow(); shg;
        end
    end

    % fit setup
    A = sparse([A_RWLS eye(length(Frame)); zeros(size(A_RWLS)) lambda*Ibig]);
    y_sp = [Frame; zeros(size(Frame))];

    % Least Squares Regression (LS)
%         x_mod = A \ y_sp;

    % Robust Regression: Iterative Reweighted Least Squares (RWLS)
    [x_mod,~,~] = RWLS(A,y_sp,epsilon,maxRWLSIter,length(Frame));

    % model fitting results
    x = x_mod(1:size(A_RWLS,2)); % detected components amplitudes
    base = x_mod(size(A_RWLS,2)+1:end); % fitted baseline

    for ll = 1:length(x)
        detectSICList(ll).amp = x(ll);
    end

    if(plotting)
        % plot current iteration's fit to frame
        temp1 = A_RWLS*x + base;
        temp2 = base;
%         subplot 514; plot(tt,Frame(offset:end-offset-1)); hold on; plot(tt,temp1(offset:end-offset-1)); plot(tt,temp2(offset:end-offset-1)); hold off;
%         axis([min(tt) max(tt) min(Frame(offset:end-offset-1)) max(Frame(offset:end-offset-1))]); drawnow(); shg;
        subplot 514; plot(tt,Frame); hold on; plot(tt,temp1); plot(tt,temp2); hold off;
        axis([min(tt) max(tt) min(Frame) max(Frame)]); drawnow(); shg;
    end
    %% Interference Cancellation Step
    dFrame = Frame - ( bLine + A_RWLS*x );
    if(plotting)
        figure(f1)
        subplot 515; plot(tt,dFrame);
        axis([min(tt) max(tt) min(dFrame) max(dFrame)]);
        drawnow(); shg;
    end

    if(any(x <= thr))
        % stopping if any of the fit detections are too small to really
        % be detections
        disp('Breaking b/c below noise floor');
        break;
    end
end
%% Final Demo Data System Model fitting
A_RWLS = [];
prunedList = [];

for ss = 1:size(detectSICList,1)
    if(detectSICList(ss).amp > thr && detectSICList(ss).pos > offset+250 && detectSICList(ss).pos < blockSize+offset-250)
        prunedList = [prunedList; detectSICList(ss)];
        A_RWLS = [A_RWLS detectSICList(ss).spfilt];
    end
end
A = sparse([A_RWLS eye(length(Frame)); zeros(size(A_RWLS)) lambda*Ibig]);
y = [Frame; zeros(size(Frame))];

% Iterative Reweighted Least Squares (RWLS)
[x_mod,~,~] = RWLS(A,y,epsilon,maxRWLSIter,length(Frame));
x = x_mod(1:size(A_RWLS,2));
base = x_mod(size(A_RWLS,2)+1:end);

%% Display Model Fitting Results
if(plotting && ~isempty(x)); figure; subplot 211, plot(Frame); hold on; plot(A_RWLS*x+base); title('Model Fit'), subplot 212, plot(Frame - A_RWLS*x - base); title('Residual'), end;
disp('Final Detected Particle Pulse Heights: ');
disp(num2str(x'));
