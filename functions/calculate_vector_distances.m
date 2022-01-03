function D = calculate_vector_distances(X)
    %ABOUT: This function takes rank vectors and calculates distances,
    %treating each vector as a description of a point in n-dimensional
    %space.
    %
    %INPUTS:
    %   1. X = a matrix that is n x m in size, where n is the number of
    %           dimensions, and m is the number of points to compare.
    %OUTPUTS:
    %   1. D = a matrix of distances calculated between points, of size
    %           [m x m], with a diagonal of 0s.
    
    [~,m] = size(X);
    D = zeros(m,m);
    for i = 1:m %Index of point 1
        for j = i+1:m %Index of point 2
            vec_1 = X(:,i);
            vec_2 = X(:,j);
            dist = sqrt(sum((vec_1 - vec_2).^2));
            D(i,j) = dist;
            D(j,i) = dist;
        end    
    end
end