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

function [x,enorm1,enorm2,xi] = RWLS(A,y,epsilon,maxIter,rwElemLen)
%RWLS Summary of this function goes here
%   Detailed explanation goes here

    % init weights
%     W = eye(length(y));
    W = ones(length(y),1); % no need for matrix
    for ii = 1:maxIter
        disp(['IRLS Iteration: ' num2str(ii)]);
        
        % fit
%         x = (W.^(1/2)*A)\(W.^(1/2)*y);
        x = bsxfun(@times, sqrt(W), A) \ (sqrt(W).*y); % faster than line above
%         for jj = 1:rwElemLen
%             W(jj,jj) = 1/(abs(A(jj,:)*x-y(jj)) + epsilon);
%         end
        W(1:rwElemLen) = 1 ./ (abs(A(1:rwElemLen,:) .* x - y(1:rwElemLen)) + epsilon .* ones(rwElemLen, 1)); % no need for loop
        
        % error
        enorm1(ii) = norm(A*x - y,2);
        enorm2(ii) = norm(A*x - y,1);
        xi(ii) = x(1);
        
    end
end

