function [latsInDegrees, longsInDegrees] = Cartesian2Spherical(Xs, Ys, Zs)
% Converts Cartesian coordinates to latitudes and longitudes.
% Cartesian coordinates should lie on a unit sphere. The origin of the
% coordinate system is at the center of the sphere. The positive x-axis passes
% through 0 degrees longitude, the positive y-axis passes through +90 degrees
% longitude, and the positive z-axis passes through +90 degrees latitude.
% Latitude and longitude are given in degrees with latitudes of [-90, 90] and
% longitudes of (-180, 180]. The point (0, 0, 0) at the center of the sphere is
% defined to have a latitude and longitude of NaN (not a number). Refer to
% http://www.geomidpoint.com/calculation.html for more information.
%
%    usage: [latsInDegrees, longsInDegrees] = Cartesian2Spherical(Xs, Ys, Zs)

    latsInDegrees = asind(Zs ./ realsqrt((Xs .^ 2) + (Ys .^ 2) + (Zs .^ 2)));

    centers = AreBadValues(latsInDegrees);
    % Create a mask for points for which division by zero occurred. These
    % points correspond to the center of the sphere.

    latsInDegrees(centers) = NaN;
    % Overwrite whatever bad value indicator (i.e., NaN + NaNi) resided in these
    % points with the more appropriate indiciator, NaN (not a number).

    longsInDegrees = Radians2Degrees(atan2(Ys, Xs));

    longsInDegrees(centers) = NaN;
    % Overwrite the longitude values corresponding to bad longitude values with
    % NaN (not a number). This prevents the caller from thinking the longitudes
    % of points at the center of the sphere are defined.

end
