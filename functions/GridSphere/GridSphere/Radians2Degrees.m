function anglesInDegrees = Radians2Degrees(anglesInRadians)
% Converts the given angles in radians to equivalent angles in degrees.
% Expects anglesInRadians to be floating point values.
% 
%   usage: anglesInDegrees = Radians2Degrees(anglesInRadians)

    DEGREES_PER_RADIAN = 180 / pi();
    anglesInDegrees = anglesInRadians * DEGREES_PER_RADIAN;
    % There are 180 / pi radians in a degree.

end
