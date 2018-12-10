function f = cost_function(M)
%% DESCRIPTION
%
%  Given the constellation parameters and phase angles, plot the
%  constellation layout in the form given for GPS in Green et. Al 1989.
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    October 15, 2015
%  Last Modified: October 15, 2015
%
%% INPUTS
%
%  M = vector of mean anomalies in the form:
%      This assumes it was filled as:
%        for k = 1:constellation.num_planes
%            for j = 1:constellation.num_SV_per_plane
%
%% OUTPUTS
%
% f = cost function for the given M. This takes the form of the original
%     paper entitled, 'The GPS 21 Primary Satellite Constellation' by Green
%     et al. (1989).
%
%% PARAMETERS / SETUP

global mu

% Load the constellation configuration.
constellation_config;

% Time intervals / period to average over.
t_vec = 0:0.25:24; % [hours], this may change depending on the problem
t_vec = t_vec*3600; % [seconds]

% Threshold for PDOP. The orginal design for GPS was 10, the current 
% performance requires 6.
PDOP_max = 5;

% Threshold for number of satellites. This is taken as the minimum needed
% to compute a position solution.
SV_min = 4;

% Elevation mask angle.
el_min = 5 * pi/180; % [rad]

% Discretization for the Earth.
lats = [-90:5:90] * pi/180; % [rad]
lons = [-180:5:180] * pi/180; % [rad]

%% IMPLEMENTATION

% square the PDOP value (more computationally efficient)
PDOP_max_squared = PDOP_max^2;

% Number of time points.
N_t = length(t_vec);

% compute the mean orbital rate.
mean_orb_rate = sqrt(mu / constellation.a^3); % [rad/s]

% Initialize the array used to keep track of the availability.
map_availability_average = zeros( length(lats), length(lons) );

% For every failure mode.
for failure_index = 0:constellation.num_sats
    % The 0 index does a check when all satellites are functioning. 
    failure_index
    % Compute the map of 'constellation value' this is essentially the
    % availability given the constraints specified above.
%     map_availability = constellation_availability(constellation, ...
%         M*pi/180, t_vec, lats, lons, PDOP_max, SV_min, el_min,...
%         failure_index);
    
    map_availability = constellation_availability(constellation, ...
        M*pi/180, t_vec, lats, lons, PDOP_max_squared, SV_min, el_min, ...
        failure_index, N_t, mean_orb_rate);
    
    % Update the availability map.
    map_availability_average = ...
        map_availability_average + map_availability;
end  % end failures

% Normalize the availability map.
map_availability_average = ...
    map_availability_average / (constellation.num_sats + 1);

% Compute the cost function.
f = sum( map_availability_average(:) ) / length(lats) / length(lons);

