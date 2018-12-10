function range = RangeToNext(values, index)
% Creates a range between the value at the given index and the next value.
% The vector starts at values(index), and each subsequent element is one greater
% than the one before it. The range continues until it reaches a value within 1
% of values(index + 1), which can never be values(index + 1) itself. Produces an
% index out of bounds error if index < 1 or index >= numel(values).
%
%    usage: range = RangeToNext(values, index)

    range = values(index):(values(index + 1) - 1);

end
