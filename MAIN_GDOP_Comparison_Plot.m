%% MAIN - COMPARISON OF GDOP PERFORMANCE - ANALYSIS OF RESULTS
%
%  This code plots the data output produced by MAIN_GDOP_Comparison.m
%
%  Written by:    Tyler Reid (tyreid@stanford.edu)
%  Start Date:    November 11, 2015
%  Last Modified: December  8, 2018
%
%% SET UP WORKSPACE

clc
clear
close all


% Determine if this is mac/linux/windows. 
if ispc
    delimit = '\';
else 
    delimit = '/'; % This applies for both unix and mac. 
end

% Add path to orbit plotting tools and other required functions and data.
p_functions = genpath([pwd, delimit, 'functions']);
addpath(p_functions);

% Add path to simulation data. 
p_functions = genpath([pwd, delimit, 'results']);
addpath(p_functions);

% Directory of where to save figures. 
dir_save = [pwd,delimit,'results',delimit,'plots',delimit];

% Define figure dimensions. 
fig_height = 9;
fig_width = 12;

%% PARAMETERS

% Level of percentile used to compute the GDOP.
percentile_level = 95;

% Required availability of position solution.
availability_thresh = 95;

% Make colors for plotting.
plot_colors = lines(50);
plot_colors_core = zeros(5, 3);
% plot_colors_core = [plot_colors(7,:);
%     plot_colors(2,:);
%     plot_colors(7,:);
%     plot_colors(5,:);
%     plot_colors(7,:)]; % BeiDou

%% WALKER RESULTS

% Load the Walker Constellation results.
load('simulation_data_el_5_leo_walker.mat') % Example only. 

% Get the number of scenarios.
[N_alt, N_sats] = size(GDOP_map_walker);

% Make figure.
figure; hold all;

% Produce the curve for each altitude and label it.
for i = 1:N_alt
    
    % Initialize arrays.
    num_sats = zeros(N_sats, 1);
    GDOP = NaN(N_sats, 1);
    
    for j = 1:N_sats
        
        % Get the actual number of satellites.
        num_sats(j) = actual_number_satellites{i, j};
        
        % Compute the GDOP.
        availability_test = sum(~isnan(GDOP_map_walker{i, j}(:))) /...
            numel(GDOP_map_walker{i, j}(:)) * 100;
        if availability_test >= availability_thresh
            GDOP(j) = ...
                prctile(real(GDOP_map_walker{i, j}(:)), percentile_level);
        end
        
%         figure;
%         hist(GDOP_map_walker{i, j}(:), 100)
        
    end % end i
    
    % Plot the curve for the particular altitude.
    plot(num_sats, GDOP, 'k--', 'color', ...
        plot_colors(i, :), 'LineWidth', 3)
    hold on
    
end % end j

%% CORE CONSTELLATIONS

% Load the core constellation data.
load('simulation_data_el_5_core_constellations.mat'); % Example only. 

% Give the names of each of the core constellations.
core_constellation_text = {'GPS', 'Galileo', 'GLONASS', 'BeiDou'};

% Get the number of core constellations.
N_core_constellations = length(GDOP_map_core_constellation);

% Define label types.
core_label = {'^', 'o', 'p', 'x'};
core_label_size = [10, 10, 12, 12];
core_line_width = [1, 2, 1, 2];

% Over plot these with their label.
for k = 1:N_core_constellations
    
    % Compute GDOP.
    GDOP = prctile(GDOP_map_core_constellation{k}(:), percentile_level);
    
    % Get the number of satellites in the constellation.
    num_sats = number_satellites_core_constellation(k);
    
    % Plot histogram. 
%     figure;
%     hist(GDOP_map_core_constellation{k}(:))
%     
    % Plot point and label.
    if k == 2
        plot( num_sats, GDOP, core_label{k}, 'color', ...
            plot_colors_core(k, :), ...
            'MarkerSize', core_label_size(k),...
            'LineWidth', core_line_width(k) );
    else
        plot( num_sats, GDOP, core_label{k}, 'color', ...
            plot_colors_core(k, :), ...
            'MarkerFaceColor', plot_colors_core(k, :), ...
            'MarkerSize', core_label_size(k),...
            'LineWidth', core_line_width(k) );
    end
    
end % end k

% legend(core_constellation_text{2:end})
% legend boxoff

%% ALTITUDE LABELS

text(288, 24.17 + 3, '500 km', 'color', plot_colors(1,:),...
    'HorizontalAlignment', 'center', 'LineWidth', 3)
text(220-15, 17.68 + 2, '700 km', 'color', plot_colors(2,:),...
    'HorizontalAlignment', 'center')
text(129 + 2, 9.546 + 1.0, '1100 km', 'color', plot_colors(3,:),...
    'HorizontalAlignment', 'center')
text(98 - 30, 8.594 + 0.0, '1400 km', 'color', plot_colors(4,:),...
    'HorizontalAlignment', 'center')
text(66 - 0, 1.35 + 0.2, 'MEO', 'color', plot_colors(5,:),...
    'HorizontalAlignment', 'left')
text(66 - 0, 1.258 - 0.2, 'GSO', 'color', plot_colors(6,:),...
    'HorizontalAlignment', 'right', 'LineWidth', 3)

%% SPECIFIC LEOs

% Teledesic
plot(288, 1.582, 'ks', 'MarkerFaceColor', [0, 0, 0],...
    'MarkerSize', 12, 'LineWidth', 1)
text(288-3, 1.582-0.2, 'Teledesic',...
    'HorizontalAlignment', 'right', 'LineWidth', 3)

% OneWeb
plot(648, 1.183,  'kv', 'MarkerFaceColor', [0, 0, 0],...
    'MarkerSize', 12, 'LineWidth', 1)
text(648-3, 1.183-0.2, 'OneWeb',...
    'HorizontalAlignment', 'right', 'LineWidth', 3)

% SpaceX
plot(4005, 0.4589,  'kh', 'MarkerFaceColor', [0, 0, 0],...
    'MarkerSize', 12, 'LineWidth', 1)
text(4005-3, 0.4589-0.1, 'SpaceX',...
    'HorizontalAlignment', 'right', 'LineWidth', 3)

%% SAVE GDOP COMPARISON FIGURE

% Set figure limits and scales.
set(gca,'XScale','log');
set(gca,'YScale','log');
grid on
xlim([2e1, 4010])
ylim([3e-1, 40])
ylabel('Geometric Dilution of Precision (GDOP)')
xlabel('Number of Satellites')
% legend(['500', '700', '1100', '1400', 'MEO', 'GSO', core_constellation_text])

% Use exportfig function.
file_name = [dir_save, 'GDOP_Comparison.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% SET LATITUDE DEPENDANT PARAMETERS TO APPLY 

% Set the percentile to take for the data.
percentile = 95; 

% Set the legend position.
legend_pos = [0.72, 0.73, 0.1, 0.1];

% Define the latitude bins to use.
latBins = intersect(latGridInDegrees,latGridInDegrees);
latBins = -90:6:90;

%% LATITUDE DEPENDANCE - GDOP

% Create the plot.
plot_dops_latitude(number_satellites, altitudes, ...
    latGridInDegrees, latBins,...
    GDOP_map_walker, GDOP_map_core_constellation, percentile, ...
    legend_pos)

% Set y axis parameters.
ylim([0,6])
ylabel('Geometric Dilution of Precision (GDOP)')

% Use exportfig function.
file_name = [dir_save, 'Latiude_GDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% LATITUDE DEPENDANCE - PDOP

% Create the plot.
plot_dops_latitude(number_satellites, altitudes, ...
    latGridInDegrees, latBins,...
    PDOP_map_walker, PDOP_map_core_constellation, percentile,...
    legend_pos)

% Set y axis parameters.
ylim([0,6])
ylabel('Position Dilution of Precision (PDOP)')

% Use exportfig function.
file_name = [dir_save, 'Latiude_PDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% LATITUDE DEPENDANCE - VDOP

% Create the plot.
plot_dops_latitude(number_satellites, altitudes, ...
    latGridInDegrees, latBins,...
    VDOP_map_walker, VDOP_map_core_constellation, percentile,...
    legend_pos)

% Set y axis parameters.
ylim([0,6])
ylabel('Vertical Dilution of Precision (VDOP)')

% Use exportfig function.
file_name = [dir_save, 'Latiude_VDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% LATITUDE DEPENDANCE - HDOP

% Create the plot.
plot_dops_latitude(number_satellites, altitudes, ...
    latGridInDegrees, latBins,...
    HDOP_map_walker, HDOP_map_core_constellation, percentile, ...
    legend_pos)

% Set y axis parameters.
ylim([0,3])
ylabel('Horizontal Dilution of Precision (HDOP)')

% Use exportfig function.
file_name = [dir_save, 'Latiude_HDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% LATITUDE DEPENDENCE - NUMBER OF SATELLITES IN VIEW

% Redefine the percentile.
percentile = 50;

% Define the latitude bins to use.
latBins = intersect(latGridInDegrees,latGridInDegrees);
% latBins = -90:6:90;

% Set the legend position.
legend_pos = [0.775, 0.19, 0.1, 0.1];

% Create the plot.
plot_dops_latitude(number_satellites, altitudes, ...
    latGridInDegrees, latBins,...
    num_SV_map_walker, num_SV_map_core_constellation, percentile,...
    legend_pos)

% Set y axis parameters.
ylim([0.1,1000])
ylabel('Number of Satellites in View')
set(gca,'YScale','log');
% set(gca(), 'ytick', [1, 100:1:1000]);  
set (gca (), 'yticklabel', [' ',num2cell([1, 10, 100, 1000])]);

% Use exportfig function.
file_name = [dir_save, 'Latiude_Num_Sats.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% CDF - NUMBER OF SATELLITES IN VIEW

% Plot
plot_dops_cdf(number_satellites, altitudes, ...
    num_SV_map_walker, num_SV_map_core_constellation)

% Set the x-axis parameters based on the DOPS under analysis.
% xlim([0, 100])
set(gca,'XScale','log');
xlim([0.9, 1000])
xlabel('Number of Satellites in View')

% Use exportfig function.
file_name = [dir_save, 'CDF_Num_Sats.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% CDF - PDOP

% Use GPS SPS Performance Standards:
% http://www.gps.gov/technical/ps/2008-SPS-performance-standard.pdf
% This is PDOP<=6 Global Coverage >= 98%

% Plot
plot_dops_cdf(number_satellites, altitudes, ...
    PDOP_map_walker, PDOP_map_core_constellation)

% Set the x-axis parameters based on the DOPS under analysis.
xlim([0, 5])
xlabel('Position Dilution of Precision (PDOP)')

% Use exportfig function.
file_name = [dir_save, 'CDF_PDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% CDF - GDOP

% Plot
plot_dops_cdf(number_satellites, altitudes, ...
    GDOP_map_walker, GDOP_map_core_constellation)

% Set the x-axis parameters based on the DOPS under analysis.
xlim([0, 5])
xlabel('Geometric Dilution of Precision (GDOP)')

% Use exportfig function.
file_name = [dir_save, 'CDF_GDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% CDF - TDOP

% Plot
plot_dops_cdf(number_satellites, altitudes, ...
    TDOP_map_walker, TDOP_map_core_constellation)

% Set the x-axis parameters based on the DOPS under analysis.
xlim([0, 5])
xlabel('Time Dilution of Precision (TDOP)')

% Use exportfig function.
file_name = [dir_save, 'CDF_TDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% CDF - HDOP

% Plot
plot_dops_cdf(number_satellites, altitudes, ...
    HDOP_map_walker, HDOP_map_core_constellation)

% Set the x-axis parameters based on the DOPS under analysis.
xlim([0, 5])
xlabel('Horizontal Dilution of Precision (HDOP)')

% Use exportfig function.
file_name = [dir_save, 'CDF_HDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);

%% CDF - VDOP

% Plot
plot_dops_cdf(number_satellites, altitudes, ...
    VDOP_map_walker, VDOP_map_core_constellation)

% Set the x-axis parameters based on the DOPS under analysis.
xlim([0, 5])
xlabel('Vertical Dilution of Precision (VDOP)')

% Use exportfig function.
file_name = [dir_save, 'CDF_VDOP.png'];
exportfig(gcf,file_name,'height',fig_height,'width',fig_width,'fontsize',22,...
    'resolution',220);
 