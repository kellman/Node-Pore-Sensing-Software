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

clear all; close all; clc;
%% setup and load
Geometry = 'MB13';
infile = ['../../Results/' Geometry '_amplitude_stats.mat'];
load(infile);
%% display
lsresm = [];
rrresm = [];
lsresv = [];
rrresv = [];

clsresm = [];
crrresm = [];
clsresv = [];
crrresv = [];

% results = abs(results);

for ss = 1:Nsnr
    lsresm = [lsresm mean(results(ss,:,2))];
    rrresm = [rrresm mean(results(ss,:,3))];
    lsresv = [lsresv var(results(ss,:,2))];
    rrresv = [rrresv var(results(ss,:,3))];
    
    clsresm = [clsresm mean(results(ss,:,4))];
    crrresm = [crrresm mean(results(ss,:,5))];
    clsresv = [clsresv var(results(ss,:,4))];
    crrresv = [crrresv var(results(ss,:,5))];
end


lsresm = 100*lsresm./amps;
rrresm = 100*rrresm./amps;
lsresv = 10000*lsresv./(amps.^2);
rrresv = 10000*rrresv./(amps.^2);

clsresm = 100*clsresm./amps;
crrresm = 100*crrresm./amps;
clsresv = 10000*clsresv./(amps.^2);
crrresv = 10000*crrresv./(amps.^2);

resfig = figure;
errorbar(SNR,clsresm,sqrt(clsresv));
hold on,
errorbar(SNR,lsresm,sqrt(lsresv));
errorbar(SNR,rrresm,sqrt(rrresv));


axis([0 33 -7 20]);
legend('Least-Squares Regression without Imperfections','Least-Squares Regression with Imperfections','Robust Regression with Imperfections');
flabel('SNR(dB)','Percent Amplitude Error of True Amplitude',12);
ftitle([Geometry ': Amplitude Estimation'],14);

%% display all
figure;
for ll = 1:3
    %% loading
    if ll == 3
        Geometry = 'MB13';
    elseif ll == 2
        Geometry = 'MB11';
    elseif ll == 1
        Geometry = 'MB7';
    end
    infile = ['../../Results/' Geometry '_amplitude_stats.mat'];
    load(infile);
    
    %% display
    lsresm = [];
    rrresm = [];
    lsresv = [];
    rrresv = [];
    clsresm = [];
    crrresm = [];
    clsresv = [];
    crrresv = [];
    
    for ss = 1:Nsnr
        lsresm = [lsresm mean(results(ss,:,2))];
        rrresm = [rrresm mean(results(ss,:,3))];
        lsresv = [lsresv var(results(ss,:,2))];
        rrresv = [rrresv var(results(ss,:,3))];

        clsresm = [clsresm mean(results(ss,:,4))];
        crrresm = [crrresm mean(results(ss,:,5))];
        clsresv = [clsresv var(results(ss,:,4))];
        crrresv = [crrresv var(results(ss,:,5))];
    end


    lsresm = 100*lsresm./amps;
    rrresm = 100*rrresm./amps;
    lsresv = 10000*lsresv./(amps.^2);
    rrresv = 10000*rrresv./(amps.^2);
    clsresm = 100*clsresm./amps;
    crrresm = 100*crrresm./amps;
    clsresv = 10000*clsresv./(amps.^2);
    crrresv = 10000*crrresv./(amps.^2);

    subplot(1,3,ll)
    % hold on,
    errorbar(SNR,clsresm,sqrt(clsresv),'LineWidth',3.5);
    hold on,
    errorbar(SNR,lsresm,sqrt(lsresv),'LineWidth',2.75);
    errorbar(SNR,rrresm,sqrt(rrresv),'LineWidth',1.75);
    set(gca,'FontSize',14)
    axis([0 33 -7 20]);
    axis square
    flabel('SNR(dB)','Mean Percent Amplitude Error of True Amplitude',14);
    ftitle([Geometry ': Amplitude Estimation'],16);
end

subplot 131
flabel('SNR(dB)','Mean Percent Amplitude Error of True Amplitude',16);
ftitle(['MB 7: Amplitude Estimation'],18);

subplot 132
flabel('SNR(dB)','Mean Percent Amplitude Error of True Amplitude',16);
ftitle(['MB 11: Amplitude Estimation'],18);

subplot 133
flabel('SNR(dB)','Mean Percent Amplitude Error of True Amplitude',16);
ftitle(['MB 13: Amplitude Estimation'],18);

% legend only on single subplot
legend('Least-Squares Regression without Imperfections','Least-Squares Regression with Imperfections','Robust Regression with Imperfections');

set(gcf,'color','white');

% %% error figure
% figure;
% for ll = 1:3
%     %% loading
%     if ll == 3
%         Geometry = 'MB13';
%     elseif ll == 2
%         Geometry = 'MB11';
%     elseif ll == 1
%         Geometry = 'MB7';
%     end
%     infile = ['../../Results/' Geometry '_amplitude_stats.mat'];
%     load(infile);
%     
%     %% display
%     lsresm = [];
%     rrresm = [];
%     lsresv = [];
%     rrresv = [];
%     for ss = 1:Nsnr
%         lsresm = [lsresm mean(results(ss,:,2))];
%         rrresm = [rrresm mean(results(ss,:,3))];
%         temp = results(ss,:,2);
%         lsresv = [lsresv var(temp)];
%         temp = results(ss,:,3);
%         rrresv = [rrresv var(temp)];
%     end
%     
%     lsresm = lsresm;
%     rrresm = rrresm;
%     lsresv = lsresv;
%     rrresv = rrresv;
% 
%     subplot(1,3,ll)
%     errorbar(SNR,lsresm,sqrt(lsresv));
%     hold on,
%     errorbar(SNR,rrresm,sqrt(rrresv));
% %     axis([0 33 0 20]);
%     % xticklabels({'5', '10', '15', '20', '25', '30'})
%     % yticklabels({'0', '0.02',  '0.04',  '0.06',  '0.08',  '0.1', '0.12', '0.14', '0.16'})
%     legend('Least-Squares Regression','Robust Regression');
%     flabel('SNR(dB)','Percent Amplitude Error of True Amplitude',12);
%     ftitle([Geometry ': Amplitude Estimation'],14);
% end
