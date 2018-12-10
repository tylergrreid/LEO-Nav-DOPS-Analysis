function plot_dops_cdf_Boeing(number_satellites, altitudes, ...
    DOP_map_walker, DOP_map_core_constellation, DOP_map_boeing)
%% PLOT CDF OF DOPS
%
%  Compute the CDF of certain LEOs and MEOs, plot the CDF of the DOPS. This
%  is set up for this particular scenario but is meant to ease in the
%  process of plotting GDOP, TDOP, PDOP, HDOP, VDOP.
%
%  Author: Tyler Reid
%
%  Start Date: November 18, 2015
%  Last Modified: November 18, 2015
%
%  This plotting scheme is based on:
%  http://stackoverflow.com/questions/19746970/custom-scaling-on-y-axis-in-matlab
%
%% INPUTS
%
%  number_satellites = vector of number of satellites used in the analysis.
%  altitudes = altitudes under consideration [m].
%  DOP_map_walker = DOP map of the walker constellation grid.
%  DOP_map_core_constellation = DOP map of the core constellations.
%
%% OUTPUTS
%
%  The desired plot.
%
%% IMPLEMENTATION 

% Desired y_ticks.
y_axis = [0.9999 0.999 0.99 0.98 0.95 0.9 0.8 0.7 0.5...
    0.3 0.2 0.1 0.05 0.02 0.01 0.001 0.0001];

% Order in increasing way.
y_val = fliplr (y_axis);

% The Sigmoid function and its inverse
% py_x = @(x) 0.5*erf(x*sqrt(pi/8)) + 0.5; % inverse function - unused
% here.
px_y = @(x) sqrt (8/pi)*erfinv (2*x - 1);

% Tranfrom the data using inverse function
% P_b2 = px_y(P_b1);

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

% Compute the empirical CDF for the core constellations.
for k = core_keep
    % Get the DOP data.
    data = DOP_map_core_constellation{k};
    
    % Fill the gaps with really large numbers.
    %disp(['Data points with no availability ', ...
    %    num2str( sum(isnan(data(:))) + sum(isinf(data(:))) )])
    
    % Compute the CDF.
    [f,x] = ecdf(data(:));
    
    % Transform the data for plotting.
    f_transformed = px_y(f);
    
    % Resample for plotting markers.
    x_interp = interp1(...
        f_transformed(~isinf(f_transformed)),...
        x(~isinf(f_transformed)), px_y(y_axis));
    
    % Plot.
    plot(x,f_transformed*100, linestyle_plot_core{k}, ...
        'color',colors_core(k,:),...
        'LineWidth', linewidth_plot_core(k));
    if ~isempty(marker_plot_core{k})
        plot(x_interp,px_y(y_axis)*100, marker_plot_core{k},...
            'color', colors_core(k,:),...
            'MarkerFaceColor',colors_core(k,:), ...
            'MarkerSize', markersize_plot_core(k));
    end % end if
end % end k

% Compute statistics on certain LEO constellations.
% Let these be Iridium, Teledesic, OneWeb, SpaceX
% Produce the curve for each altitude and label it.
LEO_cases = [
%     1200, 600; % boeing
    700,  66;   % Iridium.
%     1400, 288;  % Teledesic.
    1100, 648;  % OneWeb.
    1100, 4000; % SpaceX.
    ];

for l = 1:length(LEO_cases)
    % Get the number of satellites index.
    n_idx = find(number_satellites == LEO_cases(l, 2));
    
    % Get the altitude index.
    alt_idx = find(altitudes/1000 == LEO_cases(l, 1));
    
    % Get the DOP data.
    data = DOP_map_walker{alt_idx, n_idx};
        
    % Fill the gaps with really large numbers.
    %disp(['Data points with no availability ', ...
    %   num2str(sum(isnan(data(:)))),' out of ', ...
    %   num2str(length(data(:)))])
    data(isnan(data)) = 1e10;
    
    % Compute the CDF.
    [f,x] = ecdf(data(:));
    x(x == 1e10) = NaN;
    
    % Transform the data for plotting.
    f_transformed = px_y(f);
    
    % Resample for plotting makers.
    x_interp = interp1(...
        f_transformed(~isinf(f_transformed)),...
        x(~isinf(f_transformed)), px_y(y_axis));
    
    % Plot.
    plot(x,f_transformed*100, linestyle_plot_leo{l}, ...
        'color',colors_leo(l,:),...
        'LineWidth', linewidth_plot_leo(l));
    if ~isempty(marker_plot_leo{l})
        plot(x_interp,px_y(y_axis)*100, marker_plot_leo{l},...
            'color', colors_leo(l,:),...
            'MarkerFaceColor', colors_leo(l,:), ...
            'MarkerSize', markersize_plot_leo(l));
    end % end if
end % end l


% Compute the empirical CDF for the Boeing constellations.
for k = 1:2
    % Get the DOP data.
    data = DOP_map_boeing{k};
    
    % Fill the gaps with really large numbers.
    %disp(['Data points with no availability ', ...
    %    num2str( sum(isnan(data(:))) + sum(isinf(data(:))) )])
    
    % Compute the CDF.
    [f,x] = ecdf(data(:));
    
    % Transform the data for plotting.
    f_transformed = px_y(f);
    
    % Resample for plotting markers.
    x_interp = interp1(...
        f_transformed(~isinf(f_transformed)),...
        x(~isinf(f_transformed)), px_y(y_axis));
    
    % Plot.
    plot(x,f_transformed*100, linestyle_plot_boeing{k}, ...
        'color',colors_boeing(k,:),...
        'LineWidth', linewidth_plot_boeing(k));
    if ~isempty(marker_plot_boeing{k})
        plot(x_interp,px_y(y_axis)*100, marker_plot_boeing{k},...
            'color', colors_boeing(k,:),...
            'MarkerFaceColor',colors_boeing(k,:), ...
            'MarkerSize', markersize_plot_boeing(k));
    end % end if
end % end k


% Turn the grid on.
grid on

% Set the y-limits based on the transformation.
ylim ([ px_y(0.0001) px_y(0.9999) ]*100);

% Tranform the desired y_tick using inverse function.
set(gca(), 'ytick', px_y (y_val)*100);  
set (gca (), 'yticklabel', num2cell(y_val*100));
ylabel('Probability [%]')

h = legend('GPS','Galileo','GLONASS','BeiDou',...
    'Iridium','OneWeb','SpaceX','Boeing Initial','Boeing Final');
legend boxoff
set(h, 'position', [0.75, 0.5, 0.1, 0.1])

