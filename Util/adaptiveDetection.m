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

function [detections,PSLR,med] = adaptiveDetection(R,dict,thr,nobauds)
    detections = [];
    for jj = 1:size(R,1)
        L = dict(jj,1);
        lwin = round(L/nobauds);
        for ii = ceil(1+L/2):floor(size(R,2)-L/2)
            win = R(jj,ii-ceil(L/2):ii+floor(L/2));
            
            % takes the main lobe from center value
            [peak] = win(round(L/2));
            
            % finds the maximum sidelobe (non-robust to multiple users)
%             [sidelobe,~] = max([win(1:round(L/2) - lwin-1) win(round(L/2) + lwin+1:end)]);

            % estimates the sidelobe energy (robust)
%             sortValues= sort(abs([win(1:round(L/2) - lwin-1) win(round(L/2) + lwin+1:end)]));
%             sidelobe = sortValues(end-round(length(sortValues)/4));
            
            
            % estimates the sidelobe energy (robust)
            sortValues= (abs([win(1:round(L/2) - lwin-1) win(round(L/2) + lwin+1:end)]));
            sidelobe = median(sortValues);
            med(jj,ii) = sidelobe;
            
            PSLR(jj,ii) = peak/sidelobe;
            if(PSLR(jj,ii) > thr)
                detections = [detections; [jj,ii]];
            end
        end
    end
end

