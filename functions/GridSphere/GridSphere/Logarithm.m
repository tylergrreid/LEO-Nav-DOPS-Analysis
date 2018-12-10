function logarithms = Logarithm(values, bases)
% Computes the element-wise logarithms of the given values.
% bases are the bases for which the logarithms are taken.
%
%    usage: logarithms = Logarithm(values, bases)

    logarithms = log(values) ./ log(bases);
    % This is a logarithmic identity.

end
