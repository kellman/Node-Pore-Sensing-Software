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

clear all; clc; warning off; close all;
GeoList = [7 11 13];
SizeList = [1 2];
f10 = figure();
thrlist = zeros(3,2,3);
%%
for vv = 1:length(GeoList)
    for ww = 1:length(SizeList)
        disp(length(SizeList)*(vv-1) + ww)
        Geometry = ['B' int2str(GeoList(vv))];
        sz = ww;
        load(['../../Results/' Geometry '_detection_analysis_single_size' num2str(sz) '_2.mat']);
        addpath '../';
        %% setup
        N = 10000;
        maxDiam = 20;
        minDiam = 0;
        bl = 3;
        alpha = linspace(minDiam,maxDiam,N);

        d10 = 10;
        L = 4000;

        % dipi10 = 0.0004249; % B13 (10um cluster)
        if(strcmp(Geometry,'B13'))
            dipi10 = 0.00034; % B13 (10um cluster)
        elseif(strcmp(Geometry,'B11'))
            dipi10 = 0.00034; % B11 (10um cluster)
        elseif(strcmp(Geometry,'B7'))
            dipi10 =  0.00035; % B7 (10um cluster)
        end
%         if(strcmp(Geometry,'B13'))
%             dipi10 = 0.0003900; % B13 (10um cluster)
%         elseif(strcmp(Geometry,'B11'))
%             dipi10 = 0.00036; % B11 (10um cluster)
%         elseif(strcmp(Geometry,'B7'))
%             dipi10 =  0.0004336; % B7 (10um cluster)
%         end
        
        
        p = [dipi10*L/d10^3 0 -1 -0.8*dipi10*L];
        r = roots(p);
        D = r(1);

        iph = @(alpha) (((alpha./bl)*D^3*L)./(D + 0.8*L*(alpha./bl))).^(1/3);
        ph = @(alpha) bl*((alpha).^3/(D^2*L))*(1-(0.8/(D^3)));
        %% setup
        % new format: length(pnp),Ntrials,Nparam,Nsizes,4
        % old format: length(pnp),Nparam,Nsizes,4,Ntrials
        % total number of detections | true positives | missed detections (false negatives) | false positives
        Results = permute(Results,[1 3 4 5 2]);
        %% 
        nstd = 0;
        fathr = 0.75;
        ss = ph(particleRange);
        cc = 1;
        PP = zeros(Nparam,Nsizes,3);
        for ii = 1:Nsizes
            for jj = 1:Nparam
                sfilt = find(getAll(rspecs,'size')==ii & getAll(rspecs,'Nparam')==jj);
        %         sfilt
                spec = rspecs(sfilt);
                Ntd = length(spec);
                fa = 0; td = 0;
                fa_amp_list = [];
                td_amp_list = [];
                for kk = 1:length(spec)
                    sl = spec(kk).list;
                    % tally false alarms and true detections for this trial
                    % expiment. If just one in middle then true positive, if any
                    % others than false alarms. Also if true detection is not 
%                     if(kk ~= 27 && kk ~= 11)
                        for ll = 1:length(sl)
                            curr = sl(ll);
                            if(abs(curr.pos - length(curr.spfilt)/2) < 20)
                                td = td + 1;
                                td_amp_list = [td_amp_list, curr.amp];
                            else
                                fa = fa + 1;
                                fa_amp_list = [fa_amp_list, curr.amp];
%                                 disp([num2str(kk) '_' num2str(ll) '_' num2str(iph(curr.amp))])
                            end
                        end
%                         disp('\n')
%                     end
                end
                PP(jj,ii,1) = td;
                PP(jj,ii,2) = Ntd - td;
                PP(jj,ii,3) = fa;
                P(jj,ii) = struct('fa',fa,'td',td,'md',Ntd-td,'falist',fa_amp_list,'tdlist',td_amp_list);
            end
        end

        %% display
        cc = 1;
        bw = 0.4;
        Nplist = [1];
        Nslist = [1]; %1:1:Nsizes;
        figure(f10);
        subplot(3,2,length(SizeList)*(vv-1) + ww)

        for ii = Nplist
            for jj = Nslist
%                 subplot(length(Nplist),length(Nslist),cc);
                cc = cc + 1;
                tdlist = real(P(ii,jj).tdlist);
                falist = real(P(ii,jj).falist);

                % plot histogram of false alarms
                histogram(iph(tdlist),'BinWidth',bw,'Normalization','Probability');
                hold on,
                H = histogram(iph(falist),'BinWidth',bw,'Normalization','Probability');
                ma = max(iph(falist));
                mfa = mean(iph(falist));
                sfa = std(iph(falist));
                disp([num2str(mfa) '_' num2str(sfa)]);
                plot((mfa+2*sfa)*ones(2,1),[0 1],'k--','LineWidth',2);

%                 thrlist = [thrlist; [particleRange(jj) ma mfa sfa]];
                thrlist(vv,ww,:) = [mfa sfa mfa+2*sfa];
                
                if(vv == 1 && ww == 1)
                    legend('True Detection Histogram','False Alarm Histogram','Cutoff Threshold');
                end
                comment = sprintf('%s = %1.2f \n%s = %1.2f','\mu', mfa,'\sigma',sfa);
                text(0.5,0.7,comment,'FontSize', 20);

                axis([0 17.5 0 0.8]); axis square;
                xticks([0 2.5 5 7.5 10 12.5 15]);

                set(gca, 'FontSize', 18)
                flabel('Detected Particle Sizes', 'Normalized Counts',22);
            end
        end
        set(gcf,'color','white');
    end
end

%% saving threshold list
save(['../../Processing/thrlist_new.mat'],'thrlist');
