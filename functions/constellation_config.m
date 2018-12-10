%% CONSTELLATION CONFIGURATION FILE
%% DESCRIPTION 
%  This file constains the orbital parameters defining the number of planes
%  and Walker Constellation parameters used in initialization of the
%  optimization problem.
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    October 15, 2015
%  Last Modified: October 22, 2015

%% DEFINE WALKER DELTA CONSTELLATION PARAMETERS

global R_e

% GPS constellation
% using Walker notation this is a 24/6/2 configuration
constellation.GPS.a = 20.180e6+R_e; % [m], double check
constellation.GPS.e = 0.0; % [-]
constellation.GPS.inc = 55*pi/180; % [rad]
constellation.GPS.omega = 0; % [rad]
constellation.GPS.RAAN_ref = 0; % [rad]
constellation.GPS.num_planes = 6;
constellation.GPS.num_SV_per_plane = 4;
constellation.GPS.repeat = 2; % repeats pattern after 2 planes
constellation.GPS.num_sats = ...
    constellation.GPS.num_planes*...
    constellation.GPS.num_SV_per_plane;

% Galileo constellation
% using Walker notation this is a 27/3/1 configuration
constellation.Galileo.a = 23.222e6+R_e; % [m], double check
constellation.Galileo.e = 0.0; % [-]
constellation.Galileo.inc = 56*pi/180; % [rad]
constellation.Galileo.omega = 0; % [rad]
constellation.Galileo.RAAN_ref = 0; % [rad]
constellation.Galileo.num_planes = 3;
constellation.Galileo.num_SV_per_plane = 10;
constellation.Galileo.repeat = 1; % repeats pattern after 1 planes
constellation.Galileo.num_sats = ...
    constellation.Galileo.num_planes*...
    constellation.Galileo.num_SV_per_plane;

% Iridium constellation
% using Walker notation this is a 66/6/2 configuration
constellation.Iridium.a = 781+R_e; % [m], double check
constellation.Iridium.e = 0.0; % [-]
constellation.Iridium.inc = 86.4*pi/180; % [rad]
constellation.Iridium.omega = 0; % [rad]
constellation.Iridium.RAAN_ref = 0; % [rad]
constellation.Iridium.num_planes = 6;
constellation.Iridium.num_SV_per_plane = 11;
constellation.Iridium.repeat = 2; % repeats pattern after 2 planes
constellation.Iridium.num_sats = ...
    constellation.Iridium.num_planes*...
    constellation.Iridium.num_SV_per_plane;

% OneWeb constellation
% using Walker notation this is a 648/18/2 (most likely)
constellation.OneWeb.a = 1100e3+R_e; % [m]
constellation.OneWeb.e = 0.0; % [-]
constellation.OneWeb.inc = 87.9*pi/180; % [rad]
constellation.OneWeb.omega = 0; % [rad]
constellation.OneWeb.RAAN_ref = 0; % [rad]
constellation.OneWeb.num_planes = 18;
constellation.OneWeb.num_SV_per_plane = 36;
constellation.OneWeb.repeat = 2; % repeats pattern after 2 planes
constellation.OneWeb.num_sats = ...
    constellation.OneWeb.num_planes*...
    constellation.OneWeb.num_SV_per_plane;

% SpaceX constellation (nneds to be revised / verified)
% using Walker notation this is a 3888/54/2 (guess - no details yet)
constellation.SpaceX.a = 1100e3+R_e; % [m]
constellation.SpaceX.e = 0.0; % [-]
constellation.SpaceX.inc = 87.9*pi/180; % [rad]
constellation.SpaceX.omega = 0; % [rad]
constellation.SpaceX.RAAN_ref = 0; % [rad]
constellation.SpaceX.num_planes = 54;
constellation.SpaceX.num_SV_per_plane = 72;
constellation.SpaceX.repeat = 2; % repeats pattern after 2 planes
constellation.SpaceX.num_sats = ...
    constellation.SpaceX.num_planes*...
    constellation.SpaceX.num_SV_per_plane;

%% COMPUTE MAPPING FROM (K,J) INDEX TO SATELLITE NUMBER AND BACK

% Get the constellations.
fields = fieldnames(constellation);

for i = 1:length(fields)

    % Get the constellation field.
    field_i = fields{i};
    
    % Set the mapping.
    constellation.(field_i).map_kj_2_num = zeros(...
        constellation.(field_i).num_planes,...
        constellation.(field_i).num_SV_per_plane );
    
    constellation.(field_i).map_num_2_kj = ...
        zeros( constellation.(field_i).num_sats, 2 );

    count = 1;
    for k = 1:constellation.(field_i).num_planes
        for j = 1:constellation.(field_i).num_SV_per_plane
            constellation.(field_i).map_num_2_kj(count,:) = [k,j];
            constellation.(field_i).map_kj_2_num(k,j) = count;
            count = count + 1;
        end
    end

end

%% COMPUTE INITIAL CONDITIONS

% Create assuming a Walker constellation. 
for i = 1:length(fields)
    
    % Get the constellation field.
    field_i = fields{i};
    
    % Initialize array
    M0 = zeros( constellation.(field_i).num_planes*...
        constellation.(field_i).num_SV_per_plane);
    M0_plane1_ref = 0;

    % Based on wikipedia but also matches SMAD (page 195).
    delta_in_plane = 360/constellation.(field_i).num_SV_per_plane;
    delta_cross_plane = constellation.(field_i).repeat * 360 /...
        constellation.(field_i).num_sats;

    count = 1;
    for k = 1:constellation.(field_i).num_planes
        for j = 1:constellation.(field_i).num_SV_per_plane

            % Compute the mean anomaly.
            M0(count) = M0_plane1_ref + delta_in_plane*(j-1) + ...
                delta_cross_plane*(k-1);

            % Put this between 0 and 360.
            M0(count) = mod(M0(count),360);

            count = count + 1;
        end % end j
    end % end k
    
    constellation.(field_i).M0 = M0;

end % end i

