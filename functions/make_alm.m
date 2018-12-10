function alm_param = make_alm(n_planes, n_per_plane, altitude, TOA, ...
    inc, polar_rules, repeat_f, repeat_f_mode)
%% DESCRIPTION
%
%  Given the number of planes, the number of satellite per plane, and the
%  constellation altitude, compute the yuma orbital elements. This assumes
%  a walker star (polar) constellation.
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    November 10, 2015
%  Last Modified: November 25, 2015
%
%  Modifications:
%  November 25 - Modified cross plane relative phasing to be optimal for
%  polar orbits. Still need to fix it for non-polar orbits as it depends on
%  the cost function being modified.
%
%% INPUTS
%
%  n_planes    = Number of orbital planes.
%  n_per_plane = Number of satellites per plane. 
%  altitude    = Orbital altitude. 
%  TOA         = Time of almanac given in GPS seconds of week.
%  inc         = Inclination [rad].
%  polar_rules = Use parameters for a Walker delta or star configuration.
%  repeat_f    = The number used in #sats(T) / num_planes(P) / repeat (f). 
%                Can be an integer between 0 and P-1
%  repeat_f_mode= Indicated whether we want to set this manually 'f_manual' 
%                or use the optimal for polar orbits 'f_polar'.
%
%% OUTPUTS
% 
%  alm_parameters = Almanac parameters following the form from MAAST. This
%      is a matrix of the form:
%      1    2    3    4     5    6       7       8          9     10  11 
%     PRN ECCEN TOA INCLIN RORA SQRT_A R_ACEN ARG_PERIG MEAN_ANOM AF0 AF1
%      -    -   sec  rad    r/s m^(1/2) rad     rad        rad     s  s/s
%
%% GLOBAL VARIABLES

global R_e

%% IMPLEMENTATION

% Compute the total number of satellites.
N_total = n_planes * n_per_plane;

% Compute the in plane mean anomaly step.
delta_M0_in_plane = 2 * pi / n_per_plane;

% Compute the step in RAAN.
if polar_rules
    delta_RAAN = pi / n_planes; % For polar orbits.
else
    delta_RAAN = 2 * pi / n_planes; % For non-polar orbits.
end

% Compute the mean anomaly offset. The relative phasing.
switch repeat_f_mode
    case 'f_polar'
    % This is optimal for the polar case but we set it apart as some
    % specify it separately. See the original paper by walker (p. 17):
    % http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix=html&identifier=AD0722776
    delta_M0_cross_plane = delta_M0_in_plane / 2;
    
    case 'f_manual'
    % TODO - figure this out! The optimal value for this for the non-polar
    % case can vary due to many factors depending on the particular cost
    % function. For now, use 2. 
    delta_M0_cross_plane = repeat_f * 2 * pi / N_total;
end


% Satellite number counter. 
sat_number = 1;

% Initialize array. 
alm_param = zeros(N_total, 11);

% Populate the orbital elements
for i = 1:n_planes
    for j = 1:n_per_plane
    
        % Satellite number. 
        alm_param(sat_number, 1) = sat_number;
        
        % Eccentricity.
        alm_param(sat_number, 2) = 0; % Assume circular.
        
        % TOA in GPS seconds of week.
        alm_param(sat_number, 3) = TOA;
        
        % Inclination (typically assumed polar). 
        alm_param(sat_number, 4) = inc;
        
        % Rate of right ascention (assume 0 for simplicity here).
        % TODO: Compute this based on Gauss variation of parameters.
        alm_param(sat_number, 5) = 0;
        
        % Semi-major axis. 
        alm_param(sat_number, 6) = sqrt(R_e + altitude); % [m^1/2]
        
        % Right ascention angle.
        alm_param(sat_number, 7) = delta_RAAN*(i - 1); % [rad]
        
        % Argument of perigee.
        alm_param(sat_number, 8) = 0; % [rad] 
        
        % Mean anomaly. 
        alm_param(sat_number, 9) = mod( ...
            delta_M0_in_plane*(j - 1) +  ...
            delta_M0_cross_plane*(i - 1), 2*pi); % [rad]
        
        % Clock parameters are not used here so are kept as 0 here. 
        
        % Update the satellite number. 
        sat_number = sat_number + 1;
        
    end % end j
end % end i
