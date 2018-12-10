function areWithin = AreWithin(ones, twos, tolerance)
% Returns true if and only if the two values lie within the given tolerance.
% ones and twos may be arrays, in which case they are compared element-wise.
%
%    usage: areWithin = AreWithin(ones, twos, tolerance)

    areWithin = (abs(ones - twos) <= tolerance);
    % abs(ones - twos) is the separation between the values.

end
