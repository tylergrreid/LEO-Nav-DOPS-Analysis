function plot_dops_latitude_Boeing(number_satellites, altitudes, ...
    latGridInDegrees, latBins,...
    data_map_walker, data_map_core_constellation, data_map_Boeing,...
    percentile, legend_pos)
%% PLOT CDF OF DOPS
%
%  Given the CDF of certain LEOs and MEOs, plot the CDF of the DOPS. This
%  is set up for this particular scenario but is meant to ease in the
%  process of plotting GDOP, TDOP, PDOP, HDOP, VDOP.
%
%  Author: Tyler Reid
%
%  Start Date: November 18, 2015
%  Last Modified: November 7, 2016
%
%  Modified specifically for the Boeing constellation analysis. 
%
%  This plotting scheme is based on:
%  http://stackoverflow.com/questions/19746970/custom-scaling-on-y-axis-in-matlab
%
%% INPUTS
%
%  number_satellites = vector of number of satellites used in the analysis.
%  altitudes = altitudes under consideration [m].
%  latGridInDegrees = List of latitudes corresponding to the points.
%  latBins = Bins to use for latitude.
%  data_map_walker = Data map of the walker constellation grid.
%  data_map_core_constellation = Data map of the core constellations.
%  percentile = desired percentile of the quantity being plotted.
%  legend_pos = legend positioning parameters.
%
%% OUTPUTS
%
%  The desired plot.
%
%% IMPLEMENTATION 

% Determine the latitude groups
lat_groups = ones(size(latGridInDegrees));
for i = 1:length(latGridInDegrees)  
    % Determine the closest bin.
    test = abs(latGridInDegrees(i) - latBins);
    [~,idx] = min(test);
    
    % Record it.
    lat_groups(i) = idx;
end % end i

% Interpolation points.
latGridInDegrees_line = [-90:2.5:90];
latGridInDegrees_markers = [-90:9:90];

% Open figure to start plotting.
figure; hold all;

% Core constellations to plot.
core_keep = [1,2,3,4];

% Colors for plotting.
colors_plot = lines(20);
colors_core = colors_plot(1:4,:); 
colors_leo = colors_plot(5:7,:); 
colors_boeing = [colors_plot(4,:); 0.5 * [1, 1, 1]];

colors_core(1,:) = [0, 0 , 0]; % GPS
colors_core(2,:) = colors_plot(2,:); % Galileo
colors_core(3,:) = colors_plot(1,:); % GLONASS
colors_core(4,:) = colors_plot(3,:); % BeiDou

% Define line styles for the core constellations.
linestyle_plot_core = {'-','-.','-','-.'};
linewidth_plot_core = [4,3,3,3]*2.5/3;
markersize_plot_core = [8,30,8,8];
marker_plot_core = {'','.','s',''};

linestyle_plot_leo = {'-','-',':','-.'};
linewidth_plot_leo = [3,3,4,3]*2.5/3;
markersize_plot_leo = [8,8,8,8];
marker_plot_leo = {'v','','','^'};

linestyle_plot_boeing = {'-.','-'};
linewidth_plot_boeing = [3,3]*2.5/3;
markersize_plot_boeing = [8,8];
marker_plot_boeing = {'^','o'};

% Make some dummy plots to spoof the legend.
for k = 1:4
    plot(100000,1000000, ...
        [linestyle_plot_core{k}, marker_plot_core{k}], ...
        'color', colors_core(k,:), ...
        'LineWidth', linewidth_plot_core(k), ...
        'MarkerFaceColor',colors_core(k,:), ...
        'MarkerSize', markersize_plot_core(k))
end % end k

for k = 1:3
    plot(100000,1000000, ...
        [linestyle_plot_leo{k}, marker_plot_leo{k}], ...
        'color', colors_leo(k,:), ...
        'LineWidth', linewidth_plot_leo(k), ...
        'MarkerFaceColor',colors_leo(k,:), ...
        'MarkerSize', markersize_plot_leo(k))
end % end k

for k = 1:2
    plot(100000,1000000, ...
        [linestyle_plot_boeing{k}, marker_plot_boeing{k}], ...
        'color', colors_boeing(k,:), ...
        'LineWidth', linewidth_plot_boeing(k), ...
        'MarkerFaceColor',colors_boeing(k,:), ...
        'MarkerSize', markersize_plot_boeing(k))
end % end k

% Compute the latitude dependence the core constellations.
for k = core_keep

    % Initialize array.
    data_lat = NaN(size(latBins));
    
    % Compute data as a function of latitude.
    for i = 1:length(latBins)
        data = data_map_core_constellation{k}(lat_groups == i,:);
        
        data_lat(i) = prctile( data(:), percentile );
    end

    % Interpolate for the data markers.
    data_lat_markers = interp1(...
        latBins, data_lat, latGridInDegrees_markers);
    
    % Interpolate for lines.
    data_lat_line = interp1(...
        latBins, data_lat, latGridInDegrees_line);
    
    % Plot.
    plot(latGridInDegrees_line, data_lat_line, linestyle_plot_core{k}, ...
        'color',colors_core(k,:),...
        'LineWidth', linewidth_plot_core(k));
    if ~isempty(marker_plot_core{k})
        plot(latGridInDegrees_markers, data_lat_markers, ...
            marker_plot_core{k},...
            'color', colors_core(k,:),...
            'MarkerFaceColor',colors_core(k,:), ...
            'MarkerSize', markersize_plot_core(k));
    end % end if
end % end k

% Compute statistics on certain LEO constellations.
% Let these be Iridium, Teledesic, OneWeb, SpaceX
% Produce the curve for each altitude and label it.
LEO_cases = [
%     1200, 600;  % Boeing
    700,  66;   % Iridium.
%     1400, 288;  % Teledesic.
    1100, 648;  % OneWeb.
    1100, 4000; % SpaceX.
    ];

for l = 1:length(LEO_cases)
    % Initialize array.
    data_lat = NaN(size(latBins));
    
    % Get the number of satellites index.
    n_idx = find(number_satellites == LEO_cases(l, 2));
    
    % Get the altitude index.
    alt_idx = find(altitudes/1000 == LEO_cases(l, 1));
    
    % Compute GDOP as a function of latitude.
    for i = 1:length(latBins)
        data = data_map_walker{alt_idx, n_idx}(lat_groups == i,:);
        
        data_lat(i) = prctile( data(:), percentile );
    end
    
    % Interpolate for the data markers.
    data_lat_markers = interp1(...
        latBins, data_lat, latGridInDegrees_markers);
    
    % Interpolate for lines.
    data_lat_line = interp1(...
        latBins, data_lat, latGridInDegrees_line);
    
    % Plot.
    plot(latGridInDegrees_line, data_lat_line, linestyle_plot_leo{l}, ...
        'color',colors_leo(l,:),...
        'LineWidth', linewidth_plot_leo(l));
    if ~isempty(marker_plot_leo{l})
        plot(latGridInDegrees_markers, data_lat_markers, ..., 
            marker_plot_leo{l},...
            'color', colors_leo(l,:),...
            'MarkerFaceColor', colors_leo(l,:), ...
            'MarkerSize', markersize_plot_leo(l));
    end % end if
end % end l

% Compute the latitude dependence of the Boeing constellations.
for k = 1:2

    % Initialize array.
    data_lat = NaN(size(latBins));
    
    % Compute data as a function of latitude.
    for i = 1:length(latBins)
        data = data_map_Boeing{k}(lat_groups == i,:);
        
        data_lat(i) = prctile( data(:), percentile );
    end

    % Interpolate for the data markers.
    data_lat_markers = interp1(...
        latBins, data_lat, latGridInDegrees_markers);
    
    % Interpolate for lines.
    data_lat_line = interp1(...
        latBins, data_lat, latGridInDegrees_line);
    
    % Plot.
    plot(latGridInDegrees_line, data_lat_line, linestyle_plot_boeing{k}, ...
        'color',colors_boeing(k,:),...
        'LineWidth', linewidth_plot_boeing(k));
    if ~isempty(marker_plot_boeing{k})
        plot(latGridInDegrees_markers, data_lat_markers, ...
            marker_plot_boeing{k},...
            'color', colors_boeing(k,:),...
            'MarkerFaceColor',colors_boeing(k,:), ...
            'MarkerSize', markersize_plot_boeing(k));
    end % end if
end % end k


% Turn the grid on.
grid on


% Tranform the desired y_tick using inverse function.
xlim([-90,90])
set(gca(), 'xtick', -90:15:90);  
xlabel('Latitude [deg]')

h = legend('GPS','Galileo','GLONASS','BeiDou',...
    'Iridium','OneWeb','SpaceX','Boeing Initial', 'Boeing Final');
legend boxoff
set(h, 'position', legend_pos)


