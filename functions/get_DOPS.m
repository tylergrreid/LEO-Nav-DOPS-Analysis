function [GDOP, PDOP, TDOP, HDOP, VDOP] = get_DOPS(G, R_l)
%% DESCRIPTION
%
%  Given the geometry matrix of the form:
%  [los1', 1
%   los2', 1
%   ...
%   losN', 1]
%  Compute the DOPs. This is a direct computation by matrix inversion. For
%  a more efficient implementation that computes GDOP only see getGDOP.m
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    November 17, 2015
%  Last Modified: November 17, 2015
%
%% INPUTS
% 
%  G        = Geometry matrix.
%  R_l = Transformation matrix from ECEF to ENU coordinates.
%
%% OUTPUTS
%
%  GDOP = Geometric Dilution of Precision. 
%  PDOP = Position Dilution of Precision.
%  TDOP = Time Dilution of Precision.
%  HDOP = Horizontal Dilution of Precision.
%  VDOP = Vertical Dilution of Precision.
%
%% IMPLEMENTATION

% NOTE: see page 208 - 209 of Per Enge's GPS book.

% Create the augmented transformation matrix needed for the similarity
% transformation.
R_l_tilde = eye(4);
R_l_tilde(1:3, 1:3) = R_l;

% Take inverse. 
H = inv(R_l_tilde*(G'*G)*R_l_tilde');

% Get the diagonal elements.
DOPs = diag(H);

% Compute GDOP.
GDOP = sqrt( sum(DOPs) );

% Compute PDOP.
PDOP = sqrt( sum(DOPs(1:3)) );

% Compute TDOP.
TDOP = sqrt( DOPs(4) );

% Compute HDOP.
HDOP = sqrt( sum(DOPs(1:2)) );

% Compute VDOP.
VDOP = sqrt( DOPs(3) );
