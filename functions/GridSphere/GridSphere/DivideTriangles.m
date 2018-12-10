function subTriangles = DivideTriangles(triangles, numDimensions, ...
                            numVerticesInATriangle, numNewTrianglesPerTriangle)
% Splits each given triangle into numNewTrianglesPerTriangle equal triangles.
% Assumes triangles is a 3-D array. Each row represents a triangle. Columns
% represent x, y, and z coordinates. Depth represents the A, B, and C triangle
% vertices. Ensures that every vertex is an A, B, and C vertex in at least one
% triangle. Work in Cartesian coordinates, projecting onto a unit sphere at each
% iteration. The resulting points are evenly spaced in great circle angle, not
% latitude and/or longitude, which ensures correct handling of the poles.
% numDimensions is the number of dimensions of space. Assumes triangles have
% numVerticesInATriangle vertices each. Does not handle arbitrary input for the
% last three parameters. Only handles numDimensions = 3,
% numVerticesInATriangle = 3, and numNewTrianglesPerTriangle = 4.
%
% -----------------------------------------------------------------------
% N. Teanby   13-01-04    Original IDL (Interactive Data Language) code
% (available at http://www.atm.ox.ac.uk/user/teanby/software.html#icos)
%
% Kurt von Laven    18-09-13    Ported simplified version to GNU Octave.
% -----------------------------------------------------------------------
%
%    usage: subTriangles = DivideTriangles(triangles, numDimensions, ...
%                            numVerticesInATriangle, numNewTrianglesPerTriangle)

    oldAs = triangles(:, :, 1);
    oldBs = triangles(:, :, 2);
    oldCs = triangles(:, :, 3);
    % Get the original triangle vertices.

    Ps = ElementWiseMean(oldAs, oldBs);
    Qs = ElementWiseMean(oldBs, oldCs);
    Rs = ElementWiseMean(oldCs, oldAs);
    % Find the midpoints of each side.

    scaling = 1 / norm(Ps(1, :));
    unitPs = scaling * Ps;
    unitQs = scaling * Qs;
    unitRs = scaling * Rs;
    % Normalize midpoints onto the surface of a unit sphere.

    newAs = [oldAs; unitPs; unitRs; unitQs];
    newBs = [unitPs; oldBs; unitQs; unitRs];
    newCs = [unitRs; unitQs; oldCs; unitPs];
    % Find the sub-triangles' vertices. Ensure that every point gets used as an
    % A, B, and C in at least one triangle.

    subTriangles = reshape([newAs, newBs, newCs], ...
        numNewTrianglesPerTriangle * Rows(triangles), numDimensions, ...
        numVerticesInATriangle);
    % Put the sub-triangles' vertices in an array of the same form as the
    % original triangles.

end
