%% TEST - plot_constellation.m
%
%  This is a test script for the plot_constellation.m function which is
%  intented to be used to create 
%
%  Written by:    Tyler Reid (tyreid@stanford.edu)
%  Start Date:    October 12, 2015
%  Last Modified: October 22, 2015 (made into a test function)
%
%% WORKSPACE SETUP

% Clear command prompt.
clc

% Clear workspace variables.
clear

% Close all figures / windows.
close all

% Add path to orbit plotting tools and other astrodynamics functions.
addpath(...
    '/Users/TylerReid/Google Drive/Stanford Current/GPS Lab/TLE Reader')

% Load physical constants such as the radius of the Earth.
physical_constants

%% INITIALIZE THE WALKER CONSTELLATION

% Load the constellations configuration file.
constellation_config;

% Plot all constellations, parameters.
dir_save = [pwd,'/Results/Constellation Plots/'];
flag_export = true;
plot_marker_size = [40,40,30,30,10];
plot_color = [0,0,1];

% Get constellations.
fields = fieldnames(constellation);

% Plot all constellations.
for i = 1:length(fields)
    % Get the field.
    field_i = fields{i};
    
    % Set the file name.
    file_name = [dir_save, field_i, '.tiff'];

    % Plot.
    plot_constellation( constellation.(field_i).M0, constellation.(field_i),...
        flag_export, file_name, plot_color, plot_marker_size(i) )
end
