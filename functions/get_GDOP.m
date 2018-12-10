function GDOP = get_GDOP(G)
%% DESCRIPTION
%
%  Given the geometry matrix of the form:
%  [los1', 1
%   los2', 1
%   ...
%   losN', 1]
%  Compute the GDOP. This is an  efficient implementation based on the 
%  analytical formula based on:
%  http://link.springer.com/article/10.1007/s10291-008-0111-2
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    October 15, 2015
%  Last Modified: October 22, 2015 - changed to GDOP (not squared)
%
%% INPUTS
% 
%  G = Geometry matrix.
%
%% OUTPUTS
%
%  GDOP = The Geometric Dilution of Precision. This is computationally more
%         efficient for this specific purpose.  
%
%% IMPLEMENTATION

% Compute the eigenvalues.
% lambda = eigs( G'*G );

% The trace of a matrix is the sum of its eigenvalues. We want the trace of
% (G'*G)^-1 as this will be our GDOP^2. The eigenvalues of G'*G are the 
% reciprocal of the eigenvalues of (G'*G)^-1. Thus the sum of the inverse 
% of the eigenvalues of G'*G is the GDOP^2. Computing the eigenvalues is 
% more efficient than computing the inverse of the matrix directly. 
% GDOP^2 = sum( 1./lambda ); That being said, we employ an analytical
% method here based on:
% http://link.springer.com/article/10.1007/s10291-008-0111-2

M = G'*G;
M2 = M*M;
M3 = M*M2;

h1 = trace(M);
h2 = trace(M2);
h3 = trace(M3);
h4 = det(M);

GDOP = sqrt( (0.5*h1^3 - 1.5*h1*h2 + h3) / 3 / h4 );

