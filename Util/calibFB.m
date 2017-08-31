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

function [code,chunks,X] = calibFB(data,M,pFB_lr,K,bl_s,bl_p,geo,plotting)
    if(nargin < 9)
        plotting = 0;
    end
    if(strcmp(geo,'B11'))
        seq = [1 0 1 0 1 0 0 1 0 1 0 1 1 0 0 1 0 1 1 0 0 1];
        len = 11;
        numBlocks = 17;
    	nobauds = [1 1 1 1 1 2 1 1 1 1 2 2 1 1 2 2 1];

    elseif(strcmp(geo,'B13'))
        nobauds = [1 1 1 1 1 1 1 1 1 2 1 1 2 1 1 2 2 2 2];
        seq = [1 0 1 0 1 0 1 0 1 0 0 1 0 1 1 0 1 0 0 1 1 0 0 1 1 0];
        len = 13;
        numBlocks = 19;
    elseif(strcmp(geo,'B13i'))
        % this is flipped because I observe it in reverse
        nobauds = fliplr([1 1 1 1 1 1 1 1 2 1 1 2 1 1 2 2 2 2 1]);
        seq = ~([1 0 1 0 1 0 1 0 1 0 0 1 0 1 1 0 1 0 0 1 1 0 0 1 1 0]);
        len = 13;
        numBlocks = 19;
    elseif(strcmp(geo,'B7'))
        seq = [1 0 1 0 1 0 0 1 0 1 1 0 0 1];
        len = 7;
        numBlocks = 11;
        nobauds = [1 1 1 1 1 2 1 1 2 2 1];
    end
    bounds = [0.45 2.0];
    %% setup
    hdata = data;
    ldata = resample(hdata,1,M);
    gstart = 50;
    ldata = ldata(gstart:end-gstart);
    baseline = ASLS(ldata,bl_s,bl_p,10);
    bldata = ldata - baseline;
    %%
    Rtemp = applyFB(bldata,pFB_lr);
    [Rsort, Rind] = sort(Rtemp(:),1,'descend');
    k_ind = Rind(1:K);
    [mi,mj] = ind2sub(size(Rtemp),k_ind);
    if(plotting)
        figure
        subplot 211, imagesc(Rtemp);
        subplot 212, scatter(mj,mi);
    end
    
    % these are the triples
    ceng = Rsort(1:K);
    %%
    % now I have the top K detections, their indices (time,tau), and their
    % correlation energies. I must now selected them in time so they aren't
    % near each other.
    new_ind = []; % transit index, arrival index, correlation energy
    for ii = 1:length(ceng)
        temp_ind = find(abs(mj(ii)-mj) < 100);
        [ca,ti] = max(ceng(temp_ind));
        if(ca == ceng(ii)) % enforce its this iteration's component (otherwise the same term will be added many times)
            new_ind = [new_ind; [ca, mi(temp_ind(ti)), mj(temp_ind(ti))]];
        end
    end
%     new_ind
    
    kept_ind = [];
    for ii = 1:size(new_ind,1)
        if(all(abs(new_ind(ii,3)-new_ind(:,3)) > 500 | abs(new_ind(ii,3)-new_ind(:,3))==0))
            kept_ind = [kept_ind ii];
        end
    end
    
    if(plotting)
        hold on, scatter(new_ind(:,3),new_ind(:,2),'g*');
        hold on, scatter(new_ind(kept_ind,3),new_ind(kept_ind,2),'k*');
    end
    table = new_ind(kept_ind,:);
    %%
    % Here I need to threshold them in amplitude to get their intertiming
    % at high sampling rate??? (maybe I didn't do this)
    X = [];
    chunks = [];
    for ii = 1:size(table,1)
%         start = M*(table(ii,3) - 300 + gstart);
%         stop =  M*(table(ii,3) + 300 + gstart);
        start = table(ii,3) - 300;
        stop = table(ii,3) + 300;
        if(start < 1 ||stop > length(ldata)), continue, end;
        temp = ldata(start:stop);
        mthr = mean(temp);

        ttcurr = find(temp>mthr);
        if(isempty(ttcurr) || ttcurr(1)-1 <= 0 || ttcurr(end)+1 >= length(temp))
            continue;
        end
        
        ttpp = temp(ttcurr-1)<mthr | temp(ttcurr+1)<mthr;
        tt = ttcurr(ttpp);


        itt = diff(tt);
%         disp(num2str(length(itt)));
%         disp(numBlocks);
        if(length(itt) == numBlocks)
            
            ttemp = itt./sum(itt);
            ctemp = struct('data',temp,'cstr',itt./sum(itt));
            chunks = [chunks; ctemp];
            X = [X ttemp];
%             disp([num2str(itt') '\n']);
            if(plotting)
                figure, plot(temp); hold on, plot([1 length(temp)],mthr*ones(2,1),'r--');
                hold on, plot(tt,mthr*ones(size(tt)),'c*');
            end
            
        end
        
        
%         if(length(itt) >= numBlocks)
%             final = [];
%             for jj = 1:length(nobauds)
%                 if(itt(jj) > 5)
%                     final = [final itt(jj)];
%                     tempratio = final(end)/itt(1);
%                     if((bounds(1)*nobauds(jj) < tempratio) && (tempratio < bounds(2)*nobauds(jj)))
%                         break;
%                     end
%                 end
%                 if(length(final)==numBlocks)
%                     break;
%                 end
%             end
%             size(final)
%             
%         end
        
    end
    
    %%
    % Then I will tabulate and normalize the inter timings and solve the
    % robust regression (l1)
    [M,N] = size(X);
    cvx_begin quiet
        variable x(M)
        res = X-repmat(x,1,N);
        minimize norm(res(:),1)
        subject to
            sum(x) == 1
    cvx_end
    code = x;
%     disp(['Final calibrated code ratios: ' num2str(code')]);
%     disp(['Check sum: ' num2str(sum(code))]);
    
end