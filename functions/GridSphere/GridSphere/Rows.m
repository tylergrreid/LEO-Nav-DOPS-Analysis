function numRows = Rows(array)
% Returns the number of rows in the given array as type uint32.
% Assumes the given value is a one-dimensional (or higher-dimensional) array.
%
%    usage: numRows = Rows(array)

    ROW_DIMENSION = 1;
    % The number of the dimension corresponding to rows.

    numRows = uint32(size(array, ROW_DIMENSION));
    % By convention, rows are the first dimension of an array. The number of
    % rows cannot be negative or fractional, so it is converted to a uint32.

end
