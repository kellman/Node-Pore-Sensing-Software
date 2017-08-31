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
addpath '../../Util/'
addpath '../../Processing/'
plotting = 0; % plot while processing
%% Generate particle size to pulse height conversion
% this allows us to convert back and forth between pulse height and
% particle diameter for our channels
Geometry = 'B7';
sz = 1; % 1 for 10 micron particles and 2 for 15 micron particles

N = 10000;
maxDiam = 20;
minDiam = 0;
bl = 3;
alpha = linspace(minDiam,maxDiam,N);

d10 = 10;
L = 4000;

% dipi10 = 0.0004249; % B13 (10um cluster)
if(strcmp(Geometry,'B13'))
    dipi10 = 0.0003900; % B13 (10um cluster)
elseif(strcmp(Geometry,'B11'))
    dipi10 = 0.00036; % B11 (10um cluster)
elseif(strcmp(Geometry,'B7'))
    dipi10 =  0.0004336; % B7 (10um cluster)
end

p = [dipi10*L/d10^3 0 -1 -0.8*dipi10*L];
r = roots(p);
D = r(1);

iph = @(alpha) (((alpha./bl)*D^3*L)./(D + 0.8*L*(alpha./bl))).^(1/3);
ph = @(alpha) bl*((alpha).^3/(D^2*L))*(1-(0.8/(D^3)));
clearvars -except iph ph Geometry sz plotting

%% setup & loading
% loading experimental specs
if(strcmp(Geometry,'B13'))
    setupB13;
elseif(strcmp(Geometry,'B11'))
    setupB11;
elseif(strcmp(Geometry,'B7'))
    setupB7;
end

load(['../' outfile]);
load(['../' infile]);
maxRWLSIter = 1;
%% locate non-overlapping high SNR pulse responses
times = getAll(detectList,'pos');
[~,p] = sort(times);
detectList = detectList(p);
times = getAll(detectList,'pos');
tt_samp = (fs/M)*getAll(detectList,'ttest');
lower_bound = times-tt_samp/2/2;
upper_bound = times+tt_samp/2/2;

overLapInd = zeros(size(times));
for ii = 1:length(times)
    % given transit time, find other particles that overlap with current
    current = detectList(ii);
    if(ii == 1 & upper_bound(ii) >= lower_bound(ii+1))
        overLapInd(ii) = 1;
    elseif(ii == length(times) & lower_bound(ii) <= upper_bound(ii-1))
        overLapInd(ii) = 1;
    elseif(ii~=1 & ii~=length(times) & (upper_bound(ii) >= lower_bound(ii+1) || lower_bound(ii) <= upper_bound(ii-1)))
        overLapInd(ii) = 1;
    end
end

np = find(~overLapInd);

if(strcmp(Geometry,'B13') & sz == 2)
    % these indices are relative to np
%     ikeep = [1:9 11 12 14 16];
    ikeep = [1 3 4 10 11 14 16 17 24 33 38 41 42 47 50 52 56 59 61 67 69 76 85 87 94 97:98 101 102];
    % length is 29
    istart = 5000;
    istop = 15000;
    dthr_min = 0.002;
    dthr_max = inf;
    fromAll = 1;

elseif(strcmp(Geometry,'B11') & sz == 2)
    % these indices are relative to np
%     ikeep = [1:12 14:16 19 20];
    ikeep = [1:4 6:10 12 13 15:18 20:25 29 30 33:36 43 45 48:50 59 61 70 73 79];
    % length is 37
    istart = 5000;
    istop = 15000;
    dthr_min = 0.0015;
    dthr_max = inf;
    fromAll = 1;

elseif(strcmp(Geometry,'B7') & sz == 2)
    % these indices are relative to detectList
    ikeep = [1 13 16 19 22 23:28 32 37 40]; % 14
    istart = 5000;
    istop = 15000;
    dthr_min = 0.0012;
    dthr_max = inf;
    fromAll = 1;

elseif(strcmp(Geometry,'B13') & sz == 1)
    % these indices are relative to np
%     ikeep = [1:4 6:9 11:16 19:24];
    ikeep = [1 6 7 14 17 21 22 23 27 32 35:39 41 43 46 47 56 60 61 62 63 66 69 75]; % 27
    istart = 5000;
    istop = 15000;
    dthr_min = 0.0011;
    dthr_max = 0.0015;
    fromAll = 1;

elseif(strcmp(Geometry,'B11') & sz == 1)
    % these indices are relative to np
%     ikeep = [1:9 11:29 31:32 35:50 62:73 75:79 81:95];
%     ikeep = [2:11 18 20:22 25 26 30 33 43 45 58 60 67 69 71 75 78:80 86 87];
    ikeep = [1 5:10 12:15 23 26 31 32 33 34 38 39 44 49 50 61 63 64 65 88 91 100 102]; % 30
    istart = 5000;
    istop = 15000;
    dthr_min = 0.0010;
    dthr_max = 0.0014;
    fromAll = 1;
    
elseif(strcmp(Geometry,'B7') & sz == 1)
    % these indices are relative to np
    ikeep = [1 2 4:13 15 18:23 26:36 40 41]; % ikeep is 32 in length
    istart = 5000;
    istop = 15000;
    dthr_min = 0.0007;
    dthr_max = 0.0015;
    fromAll = 0;

end

% manually look for non-coincidence setting frames from all
% figure
% gg = 1;
% for ii = 1:length(detectList)
%     if(detectList(ii).amp > dthr_min & detectList(ii).amp < dthr_max)
%         temp = detectList(ii).sig;
%         plot(temp(istart:istop));
% %         temp = resample(temp,1,M);
% %         temp = temp(50:end-50);
%         title(num2str(gg));
%         gg = gg + 1;
%         waitforbuttonpress();
%     end
% end

% manually check for non-coincidence setting frames from auto
% figure
% gg = 1;
% for ii = 1:length(np)
%     if(detectList(np(ii)).amp > dthr_min & detectList(np(ii)).amp < dthr_max)
%         temp = detectList(np(ii)).sig;%(length(temp)/2-5000:length(temp)/2+5000);
%         plot(temp(length(temp)/2-5000:length(temp)/2+5000));
%         temp = resample(temp,1,M);
%         temp = temp(50:end-50);
% %         plot(temp);
% %         drawnow limitrate;
%         title(num2str(gg));
%         gg = gg + 1;
%         waitforbuttonpress();
%     end
% end

if(fromAll)
    pnp = find(getAll(detectList,'amp') > dthr_min & getAll(detectList,'amp') < dthr_max);
else
    pnp = np(getAll(detectList(np),'amp') > dthr_min & getAll(detectList(np),'amp') < dthr_max);
end
pnp = pnp(ikeep);
disp(ikeep)
disp(pnp)

% figure
% for ii = 1:length(pnp)
%     temp = detectList(pnp(ii)).sig;
%     plot(temp(istart:istop))
%     title(num2str(ii))
%     waitforbuttonpress();
% end 

%% Processing Setup
raw = 1./data';
global_start = 1;
global_stop = length(raw)-1;

% down sampling rate
M = 15; % from 50 kHz to ~3.333 kHz

% filterbank spec
baseTime = 0.2; % seconds
minRatio = 0.05;
maxRatio = 1.5;
numFilters = 100;
transitTime = baseTime*linspace(minRatio,maxRatio,numFilters);

% processing setup
blockSize = 2500;
blockOverlap = blockSize/2;

% baseline subtraction
bl_s = 1e6;
bl_p = 0.01;

% setup flags
% calib = 1; % perform system model calibration

% hdata = double(1./raw);
% 
% if(calib)
%     % perform system model calibration
%     [~,~,~,pFB_lr_org,uFB_lr_org,indexDict_lr_org] = genFB(rcode,code,fs,M,numFilters,minRatio,maxRatio,baseTime);
%     start_calib = 1;
%     stop_calib = length(hdata) - 1;
%     [rcode,~,~] = calibFB(hdata(start_calib:stop_calib),M,pFB_lr_org,10000,bl_s,bl_p,Geometry,0);
%     rcode = [rcode; mean(rcode(1:5))];
%     rcode = rcode./sum(rcode);
% end
% 
% % generate matched filterbank (ideal/calibrated)
% use saved calibrated matched filter bank from processing saved in results
[~,~,~,~,~,~,ttList] = genFB(rcode,code,fs,M,numFilters,minRatio,maxRatio,baseTime);

if(plotting)
    figure; subplot 121; imagesc(pFB_lr); title('Matched Filterbank');
    subplot 122; imagesc(uFB_lr); title('System Model');
end
%% Single Frame Setup
Ntrials = 1;

parameters = 2; %linspace(2,8,3);
Nparam = length(parameters);
particleRange = [1]; %[5 6 7 8 9 10 11 12 13 14 15];
Nsizes = length(particleRange);

beta = ph(particleRange); % multiply signal by this divided by true amplitude
Results = zeros(length(pnp),Ntrials,Nparam,Nsizes,4);

% alter processing to analyze false alarms
thr = 0; % processes everything it can detect.
offset = 1; % no need only a single frame (no overlap).
maxSICIter = 6; % temporary
plotting = 0;

% total number of detections | true positives | missed detections (false negatives) | false positives
rspecs = [];
for ff = 1:Nparam
    for kk = 1:length(pnp)
        for tt = 1:Ntrials
%             disp(num2str(detectList(kk).pos));
            temp = detectList(pnp(kk)).sig;
            dtemp = resample(temp(istart:istop),1,M);
            dtemp = dtemp(50:end-50);            
            for jj = 1:Nsizes
                %% setup data
%                 frame = boa(jj)*dtemp + sqrt(nadd(jj))*noise;
                frame = dtemp;
                bLine = ASLS(frame,bl_s,bl_p,10);
                stopping = 0;
                R = [];
                detectSICList = [];
                dFrame = frame - bLine;
                for gg = 1:maxSICIter
                    %% Matched Filterbank (adjoint)
                    R(:,:,gg) = applyFB(dFrame,pFB_lr);
                    if(plotting)
                        f1 = figure; subplot 511; plot(dFrame);% axis([min(tt) max(tt) min(dFrame) max(dFrame)]);
                        subplot 512; imagesc(R(:,:,gg)); axis xy; drawnow();
                    end

                    %% Detection
                    gdim1 = 5; idim1 = gdim1;
                    gdim2 = 2*floor(baseTime*fs/M/len); idim2 = gdim2/2;
                    bdim1 = 5*gdim1+idim1;
                    bdim2 = gdim2+idim2;
                    [detections,PSLR_2d] = adaptiveDetection(R(:,:,gg),indexDict_lr,parameters(ff),len);
                    %% selecting candidate detections
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
                        % stop if there are no more detections
                        disp('No more detections'); break;
                    end

                    if(plotting)
                        figure(f1);
                        subplot 513;
                        tempval = getAll(candidateList,'corr_value')';
                        cscale = 1-repmat(tempval,1,3)./max(tempval);
                        mergedDetections = [getAll(candidateList,'tti'); getAll(candidateList,'time')]';
                        scatter(mergedDetections(:,2)-offset,transitTime(mergedDetections(:,1)),6,cscale);
                        axis([1 length(dFrame) min(transitTime) max(transitTime)]);
%                         axis([min(tt) max(tt) min(transitTime) max(transitTime)]); hold on;
                        hold off; drawnow(); shg;
                    end
                    %% Select best detection
                    % maximum peak in frame that is not already in list
                    [y,si] = sort(abs(getAll(candidateList,'corr_value'))); % sort in order of correlation
                    cval = abs(getAll(candidateList,'corr_value')); % sort in order of correlation
                    cval = cval(si);
                    ss = length(si); % initialize from highest correlation detection
                    while(ss > 0)
                        currentDetection = struct('pos',candidateList(si(ss)).time,'ttindex',candidateList(si(ss)).tti,'pol',candidateList(si(ss)).pol,'amp',[],'correng',cval(ss),'pslr',y(ss),'blamp',[],'spfilt',[],'estBaud',[],'ttest',[],'index',ii,'sig',[]);
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
                        % couldn't find any unique detections
                        disp('Cannot find anymore detections'); break;
                    end
                    %% Signal Fit (joint baseline and pulse height)
                    % setup pruned system model
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
                        zpbuffer = zeros(size(dFrame));

                        % support matrix
                        zpbuffer(start:stop) = tempFilt;
                        detectSICList(ll).spfilt = zpbuffer;
                        detectSICList(ll).estBaud = transitTime(detectSICList(ll).ttindex);
                        A_RWLS = [A_RWLS zpbuffer];

                        if(plotting)
                            % plot current frame's detections
                            figure(f1);
                            subplot 513; hold on; plot(CD.pos-offset,transitTime(CD.ttindex),'g*'); hold off; drawnow(); shg;
                        end
                    end
                    
                    % Laplacian/Difference matrix
                    Ibig = eye(length(dFrame));
                    for zz = 1:length(Ibig)
                        Ibig(zz,:) = conv(Ibig(zz,:),[1 -2 1],'same');
                    end; Ibig(1,1) = -1; Ibig(end,end) = -1;

                    
                    % fit setup
                    A = sparse([A_RWLS eye(length(dFrame)); zeros(size(A_RWLS)) lambda*Ibig]);
                    y_sp = [frame; zeros(size(frame))];

                    % Least Squares regression (LS)
            %         x_mod = A \ y_sp;

                    % Iterative Reweighted Least Squares (RWLS)
                    [x_mod,~,~] = RWLS(A,y_sp,epsilon,maxRWLSIter,length(dFrame));

                    % results
                    x = x_mod(1:size(A_RWLS,2)); % detected components amplitudes
                    base = x_mod(size(A_RWLS,2)+1:end); % fitted baseline

                    disp(num2str(iph(x)));
                    for ll = 1:length(x)
                        detectSICList(ll).amp = x(ll);
                    end

                    if(plotting)
                        % plot current iteration's fit to frame
                        temp1 = A_RWLS*x + base;
                        temp2 = base;
                        subplot 514; plot(frame(offset:end-offset-1)); hold on; plot(temp1(offset:end-offset-1)); plot(tt,temp2(offset:end-offset-1)); hold off;
%                         axis([min(tt) max(tt) min(Frame(offset:end-offset-1)) max(Frame(offset:end-offset-1))]); drawnow(); shg;
                    end
                    %% Interference Cancellation
                    dFrame = frame - ( bLine + A_RWLS*x );
                    if(plotting)
                        figure(f1)
                        subplot 515; plot(dFrame(offset:end-offset-1));
%                         axis([min(tt) max(tt) min(dFrame(offset:end-offset-1)) max(dFrame(offset:end-offset-1))]);
                        drawnow(); shg;
                    end

                    if(any(x <= thr))
                        % stopping if any of the fit detections are too small to really
                        % be detections
                        disp('Breaking b/c below noise floor');
                        break;
                    end
                end
                %% results
                detectFinal = [];
                for ii = 1:length(detectSICList)
                    if(detectSICList(ii).amp > thr)
                        detectFinal = [detectFinal detectSICList(ii)];
                    end
                end
%                 detectFinal
                Results(kk,tt,ff,jj,1) = length(detectFinal);
                for ee = 1:length(detectFinal)
                    if(abs(detectFinal(ee).pos-length(frame)/2) < 15)
                        Results(kk,tt,ff,jj,2) = 1;
                    end
                end
                
                if(Results(kk,tt,ff,jj,2)==0), Results(kk,tt,ff,jj,3) = 1; end
                Results(kk,tt,ff,jj,4) = Results(kk,tt,ff,jj,1) - Results(kk,tt,ff,jj,2);
                
                rtemp = struct('list',detectFinal,'Nparam',ff,'indxdata',kk,'trial',tt,'size',jj,'amps',getAll(detectFinal,'amp'),'ttimes',ttList(getAll(detectFinal,'ttindex')));
                rspecs = [rspecs; rtemp];
            end
        end
    end
end

%% saving
% temporary
outfile = ['../../Results/' Geometry '_detection_analysis_single_size' num2str(sz) '_3.mat'];
save(outfile);
