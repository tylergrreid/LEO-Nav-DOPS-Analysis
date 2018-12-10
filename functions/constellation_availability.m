function map_availability = constellation_availability(constellation, ...
    M, t_vec, lats, lons, PDOP_max_squared, SV_min, el_min, ...
    failure_index, N_t, mean_orb_rate)
%% DESCRIPTION
%
%  Compute the constellation availability over a certain period given the
%  constraints on PDOP, number of satellites, and evelation mask angle over
%  the specified period of time.
%
%  Author:        Tyler Reid (tyreid@stanford.edu)
%  Affiliation:   Stanford University GPS Lab
%  Start Date:    October 15, 2015
%  Last Modified: October 15, 2015
%
%% INPUTS
%
%  constellation = data strcuture containing the information on the Walker
%                  constellation.
%  M             = Vector of mean anomalies [rad].
%  t_vec         = Vector of time instances for analysis [seconds].
%  lats          = Vector of latitudes for analysis [rad].
%  lons          = Vector of longitudes for analysis [rad].
%  PDOP_max_squared = Maximum acceptable Position Dilution Of Precision.
%  SV_min        = Minimum number of satellites required in view.
%  el_min        = Minimum satellite elevation angle [rad].
%  failure_index = Index or indices of satellite under failure.
%  N_t           = Number of time intervals.
%  mean_orb_rate = Mean orbital rate n [rad/s].
%
%% OUTPUTS
%
%  map_availability = This is a 2D map (lat/long) of availability for the
%                     specified time interval and restrictions on PDOP,
%                     number of satellites, and elevation mask.
%
%% IMPLEMENTATION

global R_e omega_e

% Initialize the array to keep track of the average 'constellation value'.
% This is effectively the availability of the constellation over a given
% period of time under certain performance constraints.
map_availability = zeros( length(lats), length(lons) );

% For every point on Earth.
for x = 1:length(lons)
    for y = 1:length(lats)
        
        % Compute the ECEF position vector of the site
        % (assumes a spherical Earth for the moment).
        r_site_ECEF = R_e * [
            cos(lats(y)) * cos(lons(x));
            cos(lats(y)) * sin(lons(x));
            sin(lats(y))];
        
        for n = 1:N_t % For each time index.
            
            % Compute the greenwich sidereal time.
            GMST = omega_e * t_vec(n); % [rad]
            
            % Get all the mean anomalies at epoch.
            M_t = M + mean_orb_rate * t_vec(n); % [rad]
            
            % Run through all satellites to determine coverage.
            G = []; % geometry matrix
            for sv_num = setdiff(1:constellation.num_sats,failure_index)
                
                % get the satellite plane number
                k = constellation.map_num_2_kj(sv_num,1);
                
                % Compute the ECI position of the satellite.
                [r_sat_ECI] = COE2RV_circ( constellation.a,...
                    constellation.inc,...
                    constellation.RAAN_ref + ...
                    2*pi/constellation.num_planes*(k-1),...
                    M_t(sv_num) );
                
                % Compute the ECEF position of the satellite.
                r_sat_ECEF = ECI2ECEF(r_sat_ECI, GMST);
                
                % Compute the elevation angle to the satellite
                % (simplified for spherical Earth for efficiency).
                [range, el] = rho_az_el(lats(y), lons(x), ...
                    r_site_ECEF', r_sat_ECEF');
                
                % Determine if we are in view of the site or not based
                % on our elevation mask.
                if el >= el_min
                    % Compute the line of sight vector.
                    los_vec = ( r_sat_ECEF - r_site_ECEF )./range;
                    
                    % Add this to our geometry matrix.
                    G = [G; los_vec'];
                end
                
            end % end sv_num
            
            % Default value for the availability at this instant of time
            % and location. If we meet the criteria, we bump this to an
            % availability of 1.
            constellation_value_added = 0.0;
            
            % If we meet our minimum number of satellites in view,
            % compute the PDOPs.
            [num_sv_in_view, ~] = size(G);
            if num_sv_in_view >= SV_min
                % compute the PDOPs
                PDOP = get_PDOP(G); % TODO - write this function
                
                % If the PDOPs are within our threshold, we have
                % availability at this point on the Earth at this time.
                if PDOP <= PDOP_max_squared
                    % If we meet all these conditions, this instance adds
                    % value to the constellation and is available.
                    constellation_value_added = 1.0;
                end
            end
            
            % Add the case to this case if we meet the criteria.
            map_availability(y,x) = ...
                map_availability(y,x) + constellation_value_added;
            
        end % end n (time indices)
        
    end % end y
end % end x



% normalize the 2D map of constellations values.
map_availability = map_availability / N_t;
