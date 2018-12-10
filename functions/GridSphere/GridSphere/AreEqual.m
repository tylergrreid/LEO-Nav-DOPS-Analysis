function areEqual = AreEqual(ones, twos)
% Returns true if and only if the given doubles are equal.
% Compares the given arrays element-wise. Unless extended computation has
% produced large errors or extremely dense sets of numbers are involved, doubles
% separated by values close to the machine precision usually differ only because
% of rounding errors. These values can be considered equal. If you are comparing
% scalars and are concerned about efficiency, use Equals instead.
%
%    usage: areEqual = AreEqual(ones, twos)

    TOLERANCE_FACTOR = 32;
    % The largest multiple of the machine tolerance that should be considered a
    % rounding error rather than a substantive distinction. If this value is
    % too large, then an excessively broad range of values will be considered
    % equal. If this value is too small, then an excessively small range of
    % values will be considered equal.

    maxMagnitudes = max(abs(ones), abs(twos));
    % For each pair of corresponding elements, find the larger of their
    % magnitudes.

    tolerance = TOLERANCE_FACTOR * arrayfun(@eps, maxMagnitudes);
    % For each pair of elements in the arrays, find the tolerance within which
    % two doubles should be considered equal. The built-in function eps returns
    % the machine tolerance, but can only operate on scalars, not arrays, so it
    % must be called with arrayfun.

    areEqual = AreWithin(ones, twos, tolerance);
    % Determine whether the doubles are close enough to be considered equal.

end
