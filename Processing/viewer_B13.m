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

% clear all; clc; close all; warning off;
%% loading results
filepath = '../Results/';
Geometry = 'B13';
GeoIndex = 3;

filename = [Geometry '_colloidmix_newprotein_42116_1PSI_output'];
infile = [filepath filename '.mat'];
load(infile);
load(['../Processing/thrlist_real.mat']);

% temporary
M = 15;

% sort based on starting time
times = getAll(detectList,'pos');
[~,p] = sort(times);
detectList = detectList(p);
pos = getAll(detectList,'pos');
detectList = detectList(pos<length(raw)/M/2);
%%  Transit time - pulse height scatters
% figure, subplot 211
amp = getAll(detectList,'amp');
tt = getAll(detectList,'ttest');
% 
% scatter(amp,tt,3,'b');
% title('Transit time vs. dR'); xlabel('dR (Ohms)'); ylabel('Transit time (sec)');
% 
% subplot 212
bl = abs(getAll(detectList,'blamp'));
% scatter(amp./bl,tt,3,'b');
% title('Transit time vs. dR/R'); xlabel('dR/R (unitless)'); ylabel('Transit time (sec)');

%% Convert pulse height to particle diameter
d10 = 10;
d5 = 5;
L = 4000;
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
dipi = amp./bl;
colloidDiam = ((dipi*D^3*L)./(D + 0.8*L*dipi)).^(1/3);

%% Transit time - particle Diameter scatter
% figure;
% scatter(colloidDiam,tt,3,'b');
% title('Transit time vs. Diameter (um)'); xlabel('Diameter (um)'); ylabel('Transit time (sec)');
% axis([3 16 0.05 0.25]);

%% Separate the Coincidence Detections
times = getAll(detectList,'pos');
tt_samp = (fs/M)*getAll(detectList,'ttest');
lower_bound = times-tt_samp/2;
upper_bound = times+tt_samp/2;

overLapInd = zeros(size(times));
for ii = 1:length(times)
    disp(num2str(times(ii)));
    
    % given transit time, find other particles that overlap with current
    current = detectList(ii);
    if(ii == 1 && upper_bound(ii) >= lower_bound(ii+1))
        overLapInd(ii) = 1;
    elseif(ii == length(times) && lower_bound(ii) <= upper_bound(ii-1))
        overLapInd(ii) = 1;
    elseif(ii~=1 && ii~=length(times) &&(upper_bound(ii) >= lower_bound(ii+1) || lower_bound(ii) <= upper_bound(ii-1)))
        overLapInd(ii) = 1;
    end
end
np = find(~overLapInd);
p = find(overLapInd);

% temporary to look at large non-overlapping detections
% for ii = 1:length(p)
%     if(detectList(p(ii)).amp > 0.0011)
%         figure();
%         plot(detectList(p(ii)).sig)
%     end
% end



% This is individual plots
% figure;
% subplot 311;
% scatter(colloidDiam, tt, 3, 'b');
% title('All Particles'); xlabel('Diameter (um)'); ylabel('Transit time (sec)');
% axis([3 16 0.05 0.25]);
% subplot 312;
% scatter(colloidDiam(np), tt(np), 3, 'b');
% title('Single Particles'); xlabel('Diameter (um)'); ylabel('Transit time (sec)');
% axis([3 16 0.05 0.25]);
% subplot 313;
% scatter(colloidDiam(p), tt(p), 3, 'b');
% title('Coincidence Events'); xlabel('Diameter (um)'); ylabel('Transit time (sec)');
% axis([3 16 0.05 0.25]);
% set(gcf,'color','white');

%% Coincidence Setting Pruning
% this is based off of simulation on experimental data
load('phmap.mat');
particleRange = thrlist(:,1);
times = getAll(detectList,'pos');
tt_samp = (fs/M)*getAll(detectList,'ttest');
lower_bound = times-tt_samp/2;
upper_bound = times+tt_samp/2;

overLapInd = zeros(size(times));
kept = zeros(size(times));

for ii = 1:length(times)
    overList = [ii find((upper_bound > lower_bound(ii)) & ii>(1:length(times))) , find(lower_bound < upper_bound(ii) & ii<(1:length(times)))];
%     disp(num2str(overList));
    % given transit time, find other particles that overlap with current
    current = detectList(ii);
    [~,mi] = min(abs(iph(max(getAll(detectList(overList),'amp')))-[10 15]));
    if(mi == 1)
        dthr = thrlist(GeoIndex,1,1) + 4*thrlist(GeoIndex,1,2);
    elseif(mi == 2)
        dthr = thrlist(GeoIndex,2,1) + 2*thrlist(GeoIndex,2,2);
    end



disp(num2str(dthr))
    
    if(length(overList) == 1)
        kept(ii) = 1;
    elseif(length(overList)>1)
        
        overLapAmps = iph(getAll(detectList(overList),'amp'));
        overPruneList = overLapAmps > dthr;
        
        kept(ii) = overPruneList(1);
        overLapInd(ii) = any(overPruneList(2:end));
 
    end
end
np = find(~overLapInd); % all non-overlapping detections are kept
p = find(overLapInd); % all overlapping detections
pkept = find(overLapInd & kept); % all kept overlapping detections

% This is the collated plots
% figure;
% scatter(colloidDiam(np), tt(np),'b.');
% hold on, scatter(colloidDiam(p), tt(p),'r.');
% hold on, scatter(colloidDiam(pkept),tt(pkept),'g.');
% title('Particle Statistics'); xlabel('Diameter (um)'); ylabel('Transit time (sec)');
% axis([3 16 0.05 0.25]);
% legend('Singleton Particles','Coincidence Particles','Corrected Coincidence Particles');

% This is individual plots

% subplot 411;
% scatter(colloidDiam, tt, 3, 'b');
% title('All Detections'); xlabel('Diameter (um)'); ylabel('Transit time (sec)');
% axis([4.5 16 0.05 0.25]);

if(exist('f10')==0)
    f10 = figure();
else
    figure(f10);
end

subplot 333;
scatter(colloidDiam(np), tt(np), 3, 'g');
set(gca, 'FontSize', 18)
% ftitle('Non-Overlapping Detections',28);
flabel('Diameter (\mum)','Transit time (s)',20);
axis([4 16 0.09 0.21]);
set(gca,'ytick',[0.10 0.15 0.20]);
set(gca,'yticklabel',{'0.10' '0.15' '0.20'});
set(gca,'xtick',[4 5 6 7 8 9 10 11 12 13 14 15 16]);
set(gca,'xticklabel',{'4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16'});

subplot 336;
scatter(colloidDiam(p), tt(p), 3, 'r');
set(gca, 'FontSize', 18)
% ftitle('Coincidence Event Detections',28); 
flabel('Diameter (\mum)','Transit time (s)',20);
axis([4 16 0.09 0.21]);
set(gca,'ytick',[0.10 0.15 0.20]);
set(gca,'yticklabel',{'0.10' '0.15' '0.20'});
set(gca,'xtick',[4 5 6 7 8 9 10 11 12 13 14 15 16]);
set(gca,'xticklabel',{'4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16'});

subplot 339;
scatter(colloidDiam(pkept), tt(pkept), 3, 'b');
set(gca, 'FontSize', 18)
% ftitle('Coincidence Event Detections with Pruning' ,28);
flabel('Diameter (\mum)','Transit time (s)',20);
axis([4 16 0.09 0.21]);
set(gca,'ytick',[0.10 0.15 0.20]);
set(gca,'yticklabel',{'0.10' '0.15' '0.20'});
set(gca,'xtick',[4 5 6 7 8 9 10 11 12 13 14 15 16]);
set(gca,'xticklabel',{'4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16'});

set(gcf,'color','white');

%% detection counts
% 15 micron detections
dc15 = sum([colloidDiam(np)>12 colloidDiam(pkept)>12]);
disp(num2str(dc15));

% 10 micron detections
dc10 = sum([(colloidDiam(np)>9.5 & colloidDiam(np)<11) (colloidDiam(pkept)>9.5 & colloidDiam(pkept)<11)]);
disp(num2str(dc10));

% 5 micron detections
dc5 = sum([(colloidDiam(np)<6) (colloidDiam(pkept)<6)]);
disp(num2str(dc5));
