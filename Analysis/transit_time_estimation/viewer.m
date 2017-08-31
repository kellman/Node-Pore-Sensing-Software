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
Geometry = 'MB11';
infile = ['../../Results/' Geometry '_transitTime_stats.mat'];
load(infile);

%% percent error of transit time

bnwd = 1.1; % binwidth is ___% error
figure
for ii = 1:3
    subplot(1,3,ii);
    histogram(100*results(ii,:,1)/times(ii),'binWidth',bnwd,'Normalization','Probability');
%     hold on, plot(100*tdf/times(ii)*ones(2,1),[0 0.5],'r--');
%     plot(-100*tdf/times(ii)*ones(2,1),[0 0.5],'r--');
    axis([-10 10 0 0.5]);
    flabel('Percent Transit Time Error','Normalized Histogram');
    ftitle(['Transit Time: ' num2str(times(ii)) ' sec']);
end
set(gcf,'color','white');
% saveas(gcf,'ttanal','svg');

% note: B11 has bais idk why yet...

%% supplemental figure
% all 3 sets of 3
m = zeros(3,3);
v = zeros(3,3);

figure
for mm = 1:3
    if mm == 1
        Geometry = 'MB7';
    elseif mm == 2
        Geometry = 'MB11';
    elseif mm == 3
        Geometry = 'MB13';
    end
    infile = ['../../Results/' Geometry '_transitTime_stats.mat'];
    load(infile);

    %% percent error of transit time

    bnwd = 1.2; % binwidth is ___% error
    for ii = 1:3
        subplot(3,3,(mm-1)*3 + ii);
        histogram(100*results(ii,:,1)/times(ii),'binWidth',bnwd,'Normalization','Probability');
    %     hold on, plot(100*tdf/times(ii)*ones(2,1),[0 0.5],'r--');
    %     plot(-100*tdf/times(ii)*ones(2,1),[0 0.5],'r--');
        axis([-10 10 0 0.75]);
        set(gca, 'FontSize', 18)
        flabel('Percent Transit Time Error','Normalized Histogram',22);
        ftitle(['Transit Time: ' num2str(times(ii),'%1.4f') ' sec'],28);
        
%         print('\mu');
        m(mm,ii) = mean(100*results(ii,:,1)/times(ii));
        v(mm,ii) = var(100*results(ii,:,1)/times(ii));
%         message = sprintf('Line #1\nThe second line.\nAnd finally a third line.');

        comment = sprintf('%s = %1.2f %% \n%s^{2} = %1.2f %%','\mu', m(mm,ii),'\sigma',v(mm,ii));
        text(4,0.6,comment,'FontSize', 18);
    end
end
set(gcf,'color','white');
