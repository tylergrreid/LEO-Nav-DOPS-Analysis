function [num_SV_map, GDOP_map, PDOP_map, ...
                TDOP_map, HDOP_map, VDOP_map] = analyze_constellation(...
                alm_param, t_max, el_mask, R_site, ECEF2ENU)
%% DESCRIPTION
%
%  Given the satellites in the constellation and their orbital parameters,
%  determine the constellation performance in terms of geometry (GDOP) and
%  availability. 
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    November 09, 2015
%  Last Modified: November 17, 2015 - modified to output all DOPS not just
%                                     GDOP.
%
%% INPUTS
% 
%  alm_parameters = Almanac parameters following the form from MAAST. This
%      is a matrix of the form:
%      1    2    3    4     5    6       7       8          9     10  11 
%     PRN ECCEN TOA INCLIN RORA SQRT_A R_ACEN ARG_PERIG MEAN_ANOM AF0 AF1
%      -    -   sec  rad    r/s m^(1/2) rad     rad        rad     s  s/s
%
%  t_max = The max time to be used in the analysis. This should be the
%          repeat period of the satellite ground tracks. 
%
%  el_mask = Elevation mask to the satellite before it is considered out of
%            view. 
%
%  R_site  = The ECEF positions of equally spaced grid (over a spherical
%            Earth) that represent the sites on Earth for the geometry /
%            availability to be computed.
%
%  ECEF2ENU = Transformation matrices from ECEF to ENU (pre-computed).
%
%% OUTPUTS
%
%   num_SV_map = This is a map for each site specified on Earth and
%                at each time specified. This just records the number of 
%                satellites in view and will be used along with GDOP to 
%                compute availability. 
%
%   GDOP_map =   This is a map for each site specified on Earth and
%                at each time specified. This records the geometry in terms
%                of Geometric Dilution of Precision. NaN implied no GDOP
%                computable due to an insufficient number of satellites.
%
%   PDOP_map  = Similar to the GDOP map but for Position Dilution of
%               Precision.
%
%   TDOP_map  = Similar to the GDOP map but for Time Dilution of
%               Precision.
%
%   HDOP_map  = Similar to the GDOP map but for Horizontal Dilution of
%               Precision.
%
%   VDOP_map  = Similar to the GDOP map but for Vertical Dilution of
%               Precision.
%
%% GLOBAL VARIABLES

global R_e mu

%% IMPLEMENTATION

% Get the number of sites on Earth in this analysis.
N_site = length(R_site);

% Compute the mean orbital rate.
a = alm_param(end, 6)^2;
n_sv = sqrt(mu / a^3);

% Compute the time step. This is based on the constellation altitude. 
n_MEO = sqrt(mu / 26560000^3); % [rad/s]
t_step_MEO = 1200; % [sec] - From MAAST which is also based on DOPS
t_step = n_MEO / n_sv * t_step_MEO; % This adjusts this for LEO such that 
                                    % the same angle is swept in the sky.

% Create a vector of times for analysis
t_vec = 0:t_step:t_max + min(alm_param(:,3));

% Number of time points.
N_t = length(t_vec);

% Initialize the availability map.
num_SV_map = NaN(N_site, N_t);

% Initialize the DOP maps.
GDOP_map = NaN(N_site, N_t);
PDOP_map = NaN(N_site, N_t);
TDOP_map = NaN(N_site, N_t);
HDOP_map = NaN(N_site, N_t);
VDOP_map = NaN(N_site, N_t);

% Compute the sin of the elevation mask angle
sin_el_mask = sin(el_mask);

for t = 1:N_t % For each time step.
    for site = 1:length(R_site) % For every location specified on Earth.
    
        % Compute all satellite ECEF positions. Simplify this by assuming
        % circular orbits as all in this analysis effectively are. This 
        % will greatly speed things up as we don't have to deal with 
        % Kepler's equation. 
        [~,R_sat,~] = alm2satposvel(t_vec(t), alm_param);
        radius_sat_squared = sum( R_sat.^2, 2 );
        rho_mask = ...
            -R_e * sin_el_mask +...
            (R_e ^ 2 * (sin_el_mask ^ 2 - 1) + radius_sat_squared).^0.5;
        
        % Find satellites in view based on page 32 in your notes. 
        % This is based on visibility and the elevation mask angle. 
        test = R_sat*R_site(site,:)';
        metric = ( R_e^2 + radius_sat_squared - rho_mask.^2 )./2;
        
        % Determine the satellites in view.
        sv_in_view = test>=metric;
        num_sv = sum(sv_in_view);
        
        % Record the number of satellites. 
        num_SV_map(site, t) = num_sv;
        
        % If there are enough satellites in view, compute the GDOP.
        if num_sv >= 4
        
            % Compute the geometry matrix. 
            
            % Get the positions of the satellites in view. 
            R_sat_hat = R_sat(sv_in_view,:);
            
            % Compute the line of sights.
            LOS = R_sat_hat - repmat(R_site(site,:), num_sv, 1);
            
            % Compute the ranges to each satellite so we can normalize the
            % LOS matrix. 
            rho = sum( LOS.^2, 2).^0.5;
            RHO = repmat(rho, 1, 3); 
            
            % Form the geometry matrix G. 
            G = [LOS./RHO, ones(num_sv,1)];
            
            % Compute the GDOP. (OLD)
            % GDOP_map(site, t) = get_GDOP(G);
            
            % Compute all DOPS.
            [GDOP, PDOP, TDOP, HDOP, VDOP] = get_DOPS(G, ECEF2ENU{site});
            
            % Save the data.
            GDOP_map(site, t) = GDOP;
            PDOP_map(site, t) = PDOP;
            TDOP_map(site, t) = TDOP;
            HDOP_map(site, t) = HDOP;
            VDOP_map(site, t) = VDOP;
            
        end % end num_sv
    
    end % end site
end % end t 

