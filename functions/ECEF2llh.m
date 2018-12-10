function [lat,long,h_ellp] = ECEF2llh(X)
%% DESCRIPTION
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%
%       Project Start Date:   March 28, 2011
%       Last updated:         April 14, 2011
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Given the spacecraft position the ECEF coordinates system, determine its
% geodetic latitude, longitude, and height above an ellipsoidal model of
% the Earth.
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
%
%       X = Position vector of the spacecraft expressed in      [length]*
%           the ECEF coordinate frame.
%
% -------------------------------------------------------------------------
% OUPUT
% -------------------------------------------------------------------------
%
%       ouput format:
%       [lat,long,h_ellp]
%
%       where:
%            lat = geodetic latitude                             [rad]
%           long = longitude                                     [rad]
%         h_ellp = height above an ellipsoidal model of Earth   *[length]
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%
% * this quantity can be expressed in either m or km or etc as long
%   as the global value of R_e (the Earth's Radius) is in consitant units.
%
%% DEFINE GLOBAL VARIABLES TO BE USED

global R_e Earth_E2

%% DETERMINE THE ALTITUDE ABOVE THE ELLIPSOIDAL EARTH (h_ellp)

% Define r_I, r_J, r_K -> for consistancy with Vallado (2007)
r_I = X(1);
r_J = X(2);
r_K = X(3);

% handle special case of undefined longitude
if r_I == 0 && r_J == 0 % you are at one of the poles
    
    long    = 0;
        
    % polar radius of the Earth:
    r_p = R_e*sqrt(1-Earth_E2);
        
    if r_K == 0 % at center 
        h_ellp = -R_e;
        lat = 0;
    elseif r_K > 0 % in Northern hemisphere
        h_ellp = r_K - r_p;
        lat = pi/2;
    else % r_K <0, in southern hemisphere
        h_ellp = -(r_K + r_p);
        lat = -pi/2;
    end
    
else
    
    r_delta = sqrt((r_I^2)+(r_J^2));
    
    % compute alpha
    s_alp = r_J/r_delta;
    c_alp = r_I/r_delta;
    
    alpha = atan2(s_alp,c_alp);
    
    % compute delta
    delta = atan(r_K/r_delta);
    
    % compute longitude
    long = alpha;
    
    % compute geodetic latitude:
    %
    phi     = delta;   % initial guess
    phi_old = 1e5;     % dummy value
    tol     = 1e-10;    % tolerence to be used in the iteration
    
    % iterate
    while abs(phi_old-phi)>tol
        
        temp = phi;
        C_E     = R_e/sqrt(1-Earth_E2*sin(temp)^2);
        phi     = atan((r_K+C_E*Earth_E2*sin(temp))/r_delta);
        phi_old = temp;
        
    end
    
    lat = phi;

    % Using the values of phi_gd and C_E determined above, compute h_ellp
    h_ellp = (r_delta/cos(phi)-C_E);

end
