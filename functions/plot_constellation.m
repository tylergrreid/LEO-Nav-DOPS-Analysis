function plot_constellation(M0, constellation,...
    flag_export, file_name, plot_color, plot_marker_size)
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
%  M0               = Initial phasing of all satellites in a column vector.
%                     This assumes it was filled as:
%                     for k = 1:constellation.num_planes
%                       for j = 1:constellation.num_SV_per_plane
%  constellation    = This is a data structure of the form:
%                     constellation.a, semi-major axis [km]
%                     constellation.e, eccentricity [-]
%                     constellation.inc, inclination [rad]
%                     constellation.omega, argument of perigee [rad]
%                     constellation.num_planes, number of planes
%                     constellation.num_SV_per_plane, num sats per plane
%                     constellation.repeat, Walker repeat pattern
%                     constellation.num_sats, total number of satellites
%  flag_export      = boolean of whether or not to export a figure.
%  plot_color       = row vector of [red,green,blue].
%  plot_marker_size = plot marker size.
%
%% OUTPUTS
%
%  Figure containing the plot.
%
%% IMPLEMENTATION

% open figure and hold
figure; hold all;

% generate points for axes and planes
t = -180:10:180;
x = t.*cos(constellation.inc);
y = t.*sin(constellation.inc);

% plot y axis with
plot(x,y,'k','Linewidth',4)

% plot x axis
max_x = floor(360 + 360/constellation.num_planes);
plot(0:1:max_x,zeros(max_x+1),'k','Linewidth',4)

% plot each plane
for k = 1:constellation.num_planes
    plot(x+k*360/constellation.num_planes,y,'k','Linewidth',2)
end

% determine plot bounds
buffer = 10;
ymin = -abs( 200 * sin(constellation.inc) );
ymax = abs(ymin);
xmin = -abs(200 * cos(constellation.inc)) - buffer;
xmax = 360 + abs(xmin) + buffer;

% y axis delta
y_delta = 40; % [deg]

% y-tick length
y_tick_length = 7;

% yaxis labels and ticks
% this is broken up this way to match what is done by GPS
for i = 0:y_delta:160
    x_text = i * cos(constellation.inc);
    y_text = i * sin(constellation.inc);
    text(x_text - 10,y_text,[num2str(i),'^o'],'HorizontalAlignment','right');
    plot([x_text,x_text + y_tick_length],...
        [y_text,y_text],'k','LineWidth',2);
end

for i = 200:y_delta:340
    x_text = -(360 - i) * cos(constellation.inc);
    y_text = -(360 - i) * sin(constellation.inc);
    text(x_text-10,y_text,[num2str(i),'^o'],'HorizontalAlignment','right');
    plot([x_text,x_text + y_tick_length],...
        [y_text,y_text],'k','LineWidth',2);
end

% y-label
h = text(-75,0,'Argument of Latitude','HorizontalAlignment','center');
set(h, 'rotation', constellation.inc*180/pi);

% x-label
text(max_x + 10,0,{'\Omega','RAAN'},'HorizontalAlignment','left')

% plot the satellites within the planes
for sv_num = 1:constellation.num_sats
    
    % get the plane number
    k = constellation.map_num_2_kj(sv_num,1);
    
    x_plot = M0(sv_num) * cos(constellation.inc)+...
        k * 360 / constellation.num_planes;
    y_plot = M0(sv_num) * sin(constellation.inc);
    
    if M0(sv_num)>180
        x_plot = -(360 - M0(sv_num)) * cos(constellation.inc)+...
            k * 360/constellation.num_planes;
        y_plot = -(360 - M0(sv_num)) * sin(constellation.inc);
    end
    
    plot(x_plot,y_plot,'k.','MarkerSize',plot_marker_size + 10)
    plot(x_plot,y_plot,'.','MarkerSize',plot_marker_size, ...
        'color',plot_color)
    
end

% set the plot bounds
set(gca,'ytick',[])
set(gca,'xtick',[])
axis off
ylim([ymin,ymax]);
xlim([xmin,xmax]);

% set the axes equal, this makes sure all text is aligned with axes. It
% also makes axis tilt representative of actual inclination.
axis equal

% export figure
if flag_export
    exportfig(gcf,file_name,'height',12,'width',18,'fontsize',22,...
        'resolution',220);
end
