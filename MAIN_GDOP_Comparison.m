%% MAIN - COMPARISON OF GDOP PERFORMANCE
%
%  This code performs a first order analysis of GDOP comparison of
%  constellations and their relative performance.
%
%  Written by:    Tyler Reid (tyreid@stanford.edu)
%  Start Date:    October 22, 2015
%  Last Modified: December 9, 2018
%  
%  Modifications:
%
%  November 25, 2015 - Fixed optimal relative phasing between planes for
%  polar orbits. Also fixed progress indicator and changed to parrellel by
%  using parfor.
%
%  December 9, 2018 - Modified for upload to github in a standalone format.
%
%  NOTES:
%  (1) This works in GPS time for simplicity in dealing with Yuma files.
%
%% WORKSPACE SETUP

% Clear command prompt.
clc

% Clear workspace variables.
clear

% Close all figures / windows.
close all

% Format screen output to long floating point.
format longG

% Determine if this is mac/linux/windows. 
if ispc
    delimit = '\';
else 
    delimit = '/'; % This applies for both unix and mac. 
end

% Add path to orbit plotting tools and other required functions and data.
p_functions = genpath([pwd, delimit, 'functions']);
addpath(p_functions);
p_data = genpath([pwd, delimit, 'data']);
addpath(p_data);

% Load physical constants such as the radius of the Earth.
physical_constants

% File name for saving variables. 
file_name_simulation_data = 'simulation_data_el_5';

% Define figure dimensions. 
fig_height = 9;
fig_width = 12;

% Start the clock.
tic;

% Compute a sidereal day. This is used in the calculation of the orbit
% ground track repear period.
day_sidereal = 2 * pi / omega_e; % [seconds]

%% RUN CONFIGURATION

% Set flags of which items to calculate. 
core_constellation_flag = true;
Optimal_All_flag = true;

%% SET THE EARTH GRID FOR GDOP CALCULATION

% Set up the grid over the Earth. This will be uniformly spaced over a
% spherical Earth for efficiency. Uniform grid spacing over lat/long is not
% uniform over the Earth's surface and makes it longer to compute.

% Create the grid of lat / long ECEF points evenly distributed on a
% spherical surface.
Number_sites = 1200
[latGridInDegrees, longGridInDegrees] = GridSphere(Number_sites);

% Plot the grid points on a map.
figure; hold all;
load coast
plot(long,lat,'k')
plot(longGridInDegrees, latGridInDegrees,'o')
grid on
xlim([-180,180])
ylim([-90,90])
axis equal

% Compute ECEF positions and plot them. This assumes a spherical Earth for
% the time being. 
R_site = R_e*[
    cos(latGridInDegrees*pi/180).*cos(longGridInDegrees*pi/180),...
    cos(latGridInDegrees*pi/180).*sin(longGridInDegrees*pi/180),...
    sin(latGridInDegrees*pi/180)];

% Plot ECEF site positions.
figure; hold all;
[x,y,z] = sphere(100);
surf(R_e*x,R_e*y,R_e*z,'EdgeColor','none')
plot3(R_site(:,1), R_site(:,2), R_site(:,3),'.','MarkerSize',30)
colormap bone
axis equal
axis off

% Plot grid of points.
dir_save = [pwd, delimit, 'results', delimit, 'plots', delimit];
file_name = [dir_save, 'grid_of_points.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,...
    'fontsize',22,'resolution',220);

% Compute the transformation matrices from ECEF to ENU. 
ECEF2ENU = cell( length(R_site) );
for k = 1:length(R_site)
    % Convert to radians.
    lat = latGridInDegrees(k)*pi/180;
    lon = longGridInDegrees(k)*pi/180;
    
    % Enter the matrix.
    ECEF2ENU{k} = [
        -sin(lon)         ,  cos(lon)         , 0;
        -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat);
         cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)];
end

disp(['Items precomputed, elapsed time is: ',num2str(toc/60),' [min]'])

%% SET SATELLITE VISIBILITY / AVAILABILITY PARAMETERS

% Elevation mask angle.
el_mask = 5*pi/180; % [rad]

%% GPS, Galileo, GLONASS, BeiDou

% Define the yuma files.
files_yuma = {'almgps_nov_17_2015.txt', ...
    'almgalileo.txt', 'almglonass_nov_17_2015.txt',...
    'almbeidouFullFuture.txt'};

% Compute or load the data for the MOPS 24 GPS constellation.
if core_constellation_flag
    
    % Initialize arrays.
    num_SV_map_core_constellation = cell( length(files_yuma), 1);
    GDOP_map_core_constellation = cell( length(files_yuma), 1);
    PDOP_map_core_constellation = cell( length(files_yuma), 1);
    TDOP_map_core_constellation = cell( length(files_yuma), 1);
    HDOP_map_core_constellation = cell( length(files_yuma), 1);
    VDOP_map_core_constellation = cell( length(files_yuma), 1);
    number_satellites_core_constellation = zeros(length(files_yuma), 1);
    
    for k = 1:length(files_yuma)

        % Get orbital parameters from yuma file. This makes use of the 
        % reader from MAAST. This has form:
        %  1    2    3    4     5    6       7       8          9     10  11 
        % PRN ECCEN TOA INCLIN RORA SQRT_A R_ACEN ARG_PERIG MEAN_ANOM AF0 AF1
        %  -    -   sec  rad    r/s m^(1/2) rad     rad        rad     s  s/s
        if k == 3 % GLONASS
            % Number of planes.
            n_planes = 3;
            
            % Number of satellites per plane.
            n_per_plane = 8;
            
            % Orbital altitude.
            height = 19100e3; % [m]
            
            % Time of Almanac (GPS seconds of week) - doesn't matter here.
            TOA = 0;
            
            % Inclination.
            inc = 64.8*pi/180; % [rad]
            
            % Do not apply polar rules.
            polar_rules = false;
            
            % Repeat of pattern.
            repeat_f = 1;
            repeat_f_mode = 'f_manual';
            
            % Make an idealized Walker constellation.
            alm_param = make_alm(...
                n_planes, n_per_plane, height, TOA, inc, ...
                polar_rules, repeat_f, repeat_f_mode);
        else % Otherwise read the yuma file.
            % Read yuma file. 
            alm_param = read_yuma(files_yuma{k});
        end
        
        % Compute the orbital period.
        a = alm_param(end, 6)^2;
        period = 2 * pi * sqrt(a ^ 3 / mu); % [sec]
        
        % Define the time interval for analysis. This is based on the 
        % repeat period of the orbit.
        [repeat,~] = rat(period / day_sidereal,1e-4);
        t_max = 24*3600*repeat; % [sec]

        % Sanity check plot of orbits.
        % [~,R_sat,~] = alm2satposvel(0, alm_param);
        % figure; hold all;
        % [x,y,z] = sphere(100);
        % surf(R_e*x,R_e*y,R_e*z,'EdgeColor','none')
        % plot3(R_sat(:,1), R_sat(:,2), R_sat(:,3), 'k*','MarkerSize',15)
        % axis equal
        % title(files_yuma{k})
        
        % Compute GDOP and availability maps.
        [num_SV_map, GDOP_map, PDOP_map, ...
                TDOP_map, HDOP_map, VDOP_map] = analyze_constellation(...
                alm_param, t_max, el_mask, R_site, ECEF2ENU);

        % Store in a larger array.
        num_SV_map_core_constellation{k} = num_SV_map;
        GDOP_map_core_constellation{k} = GDOP_map;
        PDOP_map_core_constellation{k} = PDOP_map;
        TDOP_map_core_constellation{k} = TDOP_map;
        HDOP_map_core_constellation{k} = HDOP_map;
        VDOP_map_core_constellation{k} = VDOP_map;
        number_satellites_core_constellation(k) = length(alm_param);

        % Output the computation time.
        disp([files_yuma{k}(1:end-4),' complete; Elapsed time is ',...
            num2str(toc/60),'  [min]']);
        
        % Output the GDOP.
        disp(['GDOP is ',num2str(prctile(GDOP_map(:),99))])

    end % end k
    
    % Save core constellation data in a .mat file.
    save([pwd,delimit,'results',delimit,'simulation_data',delimit,...
        file_name_simulation_data,'_core_constellations.mat'],...
        'num_SV_map_core_constellation', ...
        'GDOP_map_core_constellation', 'PDOP_map_core_constellation', ...
        'TDOP_map_core_constellation', ...
        'HDOP_map_core_constellation', 'VDOP_map_core_constellation', ...
        'number_satellites_core_constellation',...
        'latGridInDegrees', 'longGridInDegrees', 'R_site');         
end

%% Optimal - ALL

if Optimal_All_flag
    
    % Define the altitudes to consider.
    % 300 = ISS
    % 500 = Low LEO
    % 700 = Iridium
    % 900 = Cicada
    % 1100 = OneWeb / SpaceX / Transit
    % 1400 = Teledesic / Global Star
    % 26500 = MEO / GPS
    % 35786 = GSO
    altitudes = [500, 700, 1100, 1200, 1400, 20181.477, 35786]*1000; % [m]
    
    % Use walker delta or star pattern?
    polar_rules = true;
    
    % Define the range on the number of satellites.
    % number_satellites = [10, 30, 66, 120, 288, 648, 1400, 2500, 4000];
    number_satellites = [10, 20, 24, 30, 45, 66, 96, 120, 180, ...
        220, 288, 400, 500, 648, 800, 900, 1100, 1400, 2000, 2500, 4000];
    
        % Compute the total number of cases.
    num_cases = length(altitudes) * length(number_satellites);
    
    % Assumed polar inclination.
    inc = pi/2; % [rad]
    
    % Define repeat number. This is the number used in the convention of
    % number_sats / number_planes / repeat.
    % See: https://en.wikipedia.org/wiki/Satellite_constellation
    repeat_f = 2; % This is arbitrary at the moment, 
                  % for polar we use the optimal.
                  
    % Turn on the repeat_f mode which tells us to use the optimal for polar
    % orbits.
    repeat_f_mode = 'f_polar';
    
    % TOA, seconds of GPS week.
    TOA = 0; % [sec]
    
    % Start a count variable to keep track of how far we are.
    count_cases = 1;
    
    % Boeing request.
    % NOTE: You also need to specify the planes manually to 25.
    %altitudes = 1200*1000;
    %number_satellites = 600;
    %inc = 60*pi/180;
    %polar_rules = false;
    %repeat_f_mode = 'f_manual';
       
    % Initialize arrays. 
    num_SV_map_walker = cell( ...
        length(altitudes), length(number_satellites) );
    GDOP_map_walker = cell( ...
        length(altitudes), length(number_satellites) );
    PDOP_map_walker = cell( ...
        length(altitudes), length(number_satellites) );
    TDOP_map_walker = cell( ...
        length(altitudes), length(number_satellites) );
    HDOP_map_walker = cell( ...
        length(altitudes), length(number_satellites) );
    VDOP_map_walker = cell( ...
        length(altitudes), length(number_satellites) );
    actual_number_satellites = cell( ...
        length(altitudes), length(number_satellites) );
    plane_config = cell( ...
        length(altitudes), length(number_satellites) );
    
    % Run through all scenarios and compute the GDOP performance using the
    % optimized scheme outlined in your notes.
    for alt = 1:length(altitudes)
        for N = 1:length(number_satellites)
            
            % Output the case being started to screen as well as progress.
            disp(['Altitude ',num2str(altitudes(alt) / 1000), ...
                ' [km], Number Sats ',num2str(number_satellites(N))])
            disp(['Progress: ',...
                num2str(count_cases / num_cases * 100),' percent'])
            disp(['Elapsed time is ', num2str(toc / 60), '[min]']) 
            
            % Update the counter.
            count_cases = count_cases + 1;
            
            % TODO: Determine the bounds we should bother looking at based
            % on the spherical capp area. This is the lower bound. We also
            % shouldn't look at really high numbers either. This will just
            % help us reduce the number of cases we need to look at and the
            % overall computation time. 
            
            % We will assume a star pattern Walker constellation, i.e. a
            % Polar constellation with i = 90 degrees. To find the number
            % of planes and number of satellites per plane, refer to your
            % notes on page 57 - 58. This is somewhat simplified but good
            % enough to start with and matches well with Iridium, Teledesic
            % and OneWeb.
            %n_plane = 25; % From Boeing.
            n_planes = round( sqrt( number_satellites(N) / 2) );
            n_per_plane = round( number_satellites(N) / n_planes );
            
            % Record the number of satellites. Not all combinations are
            % truely divisible into a proper number of planes. This allows
            % us to break it down in a reasonable way. 
            actual_number_satellites{alt, N} = n_planes*n_per_plane;
            
            % Record the constellation configuration.
            plane_config{alt, N} = [n_planes, n_per_plane];
            
            % Create the almanac parameters based on this. 
            alm_param = make_alm(...
                n_planes, n_per_plane, altitudes(alt), TOA, inc, ...
                polar_rules, repeat_f, repeat_f_mode);

            % Sanity check plot of orbits.
            %[~,R_sat,~] = alm2satposvel(0, alm_param);
            %figure; hold all;
            %[x,y,z] = sphere(100);
            %surf(R_e*x,R_e*y,R_e*z,'EdgeColor','none')
            %plot3(R_sat(:,1), R_sat(:,2), R_sat(:,3), 'k*','MarkerSize',15)
            %axis equal
            
            % Compute the orbital period.
            a = alm_param(end, 6)^2;
            period = 2 * pi * sqrt(a ^ 3 / mu); % [sec]

            % Define the time interval for analysis. This is based on the 
            % repeat period of the orbit.
            % [repeat,~] = rat(period / day_sidereal,1e-4);
            repeat = 1.0; %min([14, repeat]);
            disp(['Number of simulation days: ',num2str(repeat)])
            t_max = 24*3600*repeat; % [sec]
            
            % Compute GDOP and number of satellites maps.
            [num_SV_map, GDOP_map, PDOP_map, ...
                TDOP_map, HDOP_map, VDOP_map] = analyze_constellation(...
                alm_param, t_max, el_mask, R_site, ECEF2ENU);
            
            % Record data.
            num_SV_map_walker{alt, N} = num_SV_map;
            GDOP_map_walker{alt, N} = GDOP_map;
            PDOP_map_walker{alt, N} = PDOP_map;
            TDOP_map_walker{alt, N} = TDOP_map;
            HDOP_map_walker{alt, N} = HDOP_map;
            VDOP_map_walker{alt, N} = VDOP_map;

            % Output results to screeen.
            % disp(['GDOP is ',num2str(prctile(real(GDOP_map (:)),99))])
            
        end % N       
    end % alt
    
    % Record data in a .mat file.
    save([pwd,delimit,'results',delimit,'simulation_data',delimit,...
        file_name_simulation_data,'_leo_walker.mat'],...
        'num_SV_map_walker', 'GDOP_map_walker', ...
        'PDOP_map_walker', 'TDOP_map_walker', ...
        'HDOP_map_walker', 'VDOP_map_walker', ...
        'actual_number_satellites', 'altitudes', 'number_satellites', ...
        'plane_config', 'latGridInDegrees', 'longGridInDegrees', 'R_site',...
        'repeat_f', 'repeat_f_mode')
end

