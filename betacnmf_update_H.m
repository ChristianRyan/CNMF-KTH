% Copyright (c) 2018 Pedro Villasana and Stanislaw Gorlow
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [W, H, U, cost] = betacnmf(V, W, H0, beta, nu, maxiter, showcost)
%BETACNMF Convolutional nonnegative matrix factorization under beta-divergence
%   Exact updates

%   Authors: Pedro Villasana, Stanislaw Gorlow
%   Date:    5/3/2018

M = numel(W);
%W = W0;
H = H0;
U = pconv(W, H);

cost = NaN(1, maxiter);

switch beta
  case 0
    % Itakura-Saito divergence
    div = @div_is;
  case 1
    % Generalized Kullback-Leibler divergence
    div = @div_kl;
  case 2
    % Euclidean distance
    div = @div_ls;
  otherwise
    % Beta-divergence
    div = @div_beta;
end

if showcost
  % display cost function
  f = figure('Name', 'Cost Function');
end

% initial cost
C = Inf;

for iter = 1 : maxiter
  % current cost
  D = sum(sum(div(V, U, beta)));
  D = abs(D);
  
  % relative deviation from previous cost
  delta = abs(D - C) / D;
  C = D;
  
  cost(iter) = C / numel(V);
  
  % break condition
  if delta < nu
    break
  end
  
  % multiplicative updates
  % W
  U1 = (U + eps).^(beta - 1);
  U2 = (U + eps).^(beta - 2);
  V2 = (V + eps) .* U2;
% 
%   for m = 1 : M
%     Hm = shr(H, m - 1);
%     W{m} = W{m} .* ((V2 * Hm' + eps) ./ (U1 * Hm' + eps));
%   end
%   
  U = pconv(W, H);
  
  % H
  U1 = (U + eps).^(beta - 1);
  U2 = (U + eps).^(beta - 2);
  V2 = (V + eps) .* U2;

  A = W{1}' * V2;
  B = W{1}' * U1;
  
  for m = 2 : M
    A = A + W{m}' * shl(V2, m - 1);
    B = B + W{m}' * shl(U1, m - 1);
  end
  
  H = H .* ((A + eps) ./ (B + eps));
  U = pconv(W, H);
  
  if showcost
    % display cost function
    figure(f)
    loglog(cost), xlabel('Iteration'), ylabel('Cost'), xlim([0, maxiter])
    drawnow limitrate nocallbacks
  end
end

end

function V = pconv(W, H)

M = numel(W);
V = W{1} * H;

for m = 2 : M
  V = V + W{m} * shr(H, m - 1);
end

end

function B = shl(A, k)

if k < 0
  B = shr(A, -k);
  return
end

B = A(:, 1 + k : end);
B = cat(2, B, zeros(size(A, 1), size(A, 2) - size(B, 2)));

end

function B = shr(A, k)

if k < 0
  B = shl(A, -k);
  return
end

B = A(:, 1 : size(A, 2) - k);
B = cat(2, zeros(size(A, 1), size(A, 2) - size(B, 2)), B);

end

function d = div_is(p, q, beta)

assert(beta == 0, 'Bad value for beta.')

p = p + eps;
q = q + eps;
d = p ./ q - log(p ./ q) - 1;

end

function d = div_kl(p, q, beta)

assert(beta == 1, 'Bad value for beta.')

p = p + eps;
q = q + eps;
d = p .* log(p ./ q) - p + q;

end

function d = div_ls(p, q, beta)

assert(beta == 2, 'Bad value for beta.')

d = 0.5 * (p - q).^2;

end

function d = div_beta(p, q, beta)

if beta < 0 || beta > 2
  error ('Bad value for beta.')
end

switch beta
  case 0
    d = div_is(p, q, 0);
  case 1
    d = div_kl(p, q, 1);
  case 2
    d = div_ls(p, q, 2);
  otherwise
    gamma = beta - 1;
    p = p + eps;
    q = q + eps;
    d = 1 / (beta * gamma) * (p.^beta + gamma * q.^beta - beta * p .* q.^gamma);
end

end

% [EOF]