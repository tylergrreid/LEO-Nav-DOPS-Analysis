function [rho,beta,el,R_ENU] = rho_az_el(lat,long,h_ellp,R_sat)
%% DESCRIPTION
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab

%       Last updated:         April 09, 2014
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
%    
%            lat = geodetic latitude of the site               [rad]
%           long = longitude of the site                       [rad]
%         h_ellp = height above an ellipsoidal model of Earth  [length]
%                  of the site
%          R_sat = position vector of the satellite            [length]
%                  in the ECEF frame.
%          V_sat = velocity vector of the satellite in         [length/time]
%                  the ECEF frame.
%         
% *note: vectors are assumed to be row vectors
%
% -------------------------------------------------------------------------
% OUPUT
% -------------------------------------------------------------------------
%
%
%% IMPLEMENTATION

% compute the position vector of the site in the ECEF frame
R_site = llh2ECEF(lat,long,h_ellp);

% compute the range and range rate vectors in the ECEF frame
rho_ECEF  = R_sat - R_site;
%drho_ECEF = V_sat;

% define the tranformation to the south east and up (SEZ) coordinate system
mat(1,1) =  sin(lat)*cos(long);
mat(1,2) =  sin(lat)*sin(long);
mat(1,3) = -cos(lat);

mat(2,1) = -sin(long);
mat(2,2) =  cos(long);
mat(2,3) =  0;

mat(3,1) =  cos(lat)*cos(long);
mat(3,2) =  cos(lat)*sin(long);
mat(3,3) =  sin(lat);

% transform to the SEZ coordinate system
rho_SEZ  = mat*rho_ECEF';
%drho_SEZ = mat*drho_ECEF';
R_SEZ = rho_SEZ;
R_ENU = [R_SEZ(2), -R_SEZ(1), R_SEZ(3)];

% compute the range and range rate
rho  = norm(rho_SEZ);
%drho = dot(rho_SEZ,drho_SEZ)/rho;

% compute the elevation angle
el  = asin(rho_SEZ(3)/rho);

temp  = sqrt( rho_SEZ(1)^2 + rho_SEZ(2)^2 );
%dtemp = sqrt( drho_SEZ(1)^2 + drho_SEZ(2)^2 );

%del = (drho_SEZ(3) - drho *sin(el)) / temp;

% compute the azimuth angle
if el == pi/2
    
    sinb =  drho_SEZ(2)/dtemp;
    cosb = -drho_SEZ(1)/dtemp;
    beta =  atan2(sinb,cosb); 

    
else
    
    sinb =  rho_SEZ(2)/temp;
    cosb = -rho_SEZ(1)/temp;
    beta =  atan2(sinb,cosb);
    
end

if beta < 0
    beta = beta + 2*pi;
end

%dbeta = (drho_SEZ(1)*rho_SEZ(2)-drho_SEZ(2)*rho_SEZ(1))/temp^2;


