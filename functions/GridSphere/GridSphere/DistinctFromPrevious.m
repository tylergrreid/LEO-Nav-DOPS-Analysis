function distinctFromPrevious = DistinctFromPrevious(values)
% Returns the indices of the elements that do not equal those preceding them.
% Assumes values is a numeric column vector. distinctFromPrevious is a column
% vector of type uint32. If values is sorted, then this function is equivalent
% to finding the indices of the unique elements. If there are multiple instances
% of a value, then the lowest index at which that value is found is chosen.
%
%    usage: distinctFromPrevious = DistinctFromPrevious(values)

    shiftedVals = circshift(values, 1);
    % Rotate the column vector one element forwards.

    if (isinteger(values))
        equalsPrevious = (values == shiftedVals);
    else
        equalsPrevious = AreEqual(values, shiftedVals);
    end
    % Determine whether each element is equal to the previous element in the
    % list.

    equalsPrevious(1) = false;
    % Since the circular shift compared the first element to the last element,
    % the first element could mistakenly have been marked equal to the
    % non-existent "previous" element.

    distinctFromPrevious = uint32(find(~equalsPrevious));
    % Get the indices of the values distinct from the element preceding them.

end
