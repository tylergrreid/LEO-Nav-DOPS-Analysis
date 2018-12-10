%% TEST OF get_PDOP FUNCTION
%
% DESCRIPTION
%
%  This is a test script for the get_GDOP function.
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    October 15, 2015
%  Last Modified: October 22, 2015 - changed to GDOP (not squared)
%
%% TEST

close all
clear
clc

% Test geometry matrix.
G = [
    0.408248 0.816497 0.408248 1;
    0.666667 0.333333 0.666667 1;
    -0.408248 0.816497 0.408248 1;
    0.816497 -0.408248 0.408248 1;
    0.408248 -0.816497 0.408248 1
    ];

% Define the measurement matrix.
M = G'*G;

% Compute GDOP directly by matrix inversion
GDOP_direct = sqrt(trace(inv(M)))

GDOP_function = get_GDOP(G)

% what about PDOP?
Q = inv(M);
Q = diag(Q);
PDOP_direct = sqrt(sum(Q(1:3)))
