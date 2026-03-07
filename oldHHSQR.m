function [Q,R] = oldHHSQR(A)
    
    [n,m] = size(A);
    % Stores all reflection matrices, sequentially left-multiplied.
    H_combined = eye(n);
    
    % The starting matrix is simply A. In each iteration, we are going to
    % left-multiply the submatrix by the computed reflection matrix, and
    % then crop it by 1 row and 1 column. Therefore, we will store the
    % starting matrix A in a separate variable "submatrix" and operate on
    % that only.
    submatrix = H_combined * A;
    
    % Need to compute a reflection matrix for every diagonal. Since there
    % are more columns than rows, we'll iterate over rows.
    for col_idx = 1:1:m
        
        H_curr = eye(size(H_combined)); % full reflection matrix shape
        H_small = eye(n+1-col_idx); % small reflection sub-matrix

        % Get first column of sub-matrix
        first_a_col = submatrix(:,1);
        % Compute V = +/- norm(a)*[1 0 0...]
        v = eye(size(first_a_col)) * norm(first_a_col);
        
        if first_a_col(1) >=0
            v = v + first_a_col;
        else
            v = v - first_a_col;
        end
        
        % Compute the sub-matrix reflection
        H_small = H_small - 2 * ((v * transpose(v)) / (transpose(v) * v));
        
        % Fit sub-reflection matrix into original sized reflection matrix
        H_curr(col_idx:end, col_idx:end) = H_small;
        % Left multiply with previous reflection matrices, stores progress
        H_combined = H_curr * H_combined;
        
        % Compute the result of reflection * submatrix. In next iterations, 
        % need to get the sub-submatrix by removing first row & column.
        submatrix = H_small * submatrix;
        submatrix = submatrix(2:end, 2:end);
    end

    % Final R is all the reflections * original matrix A.
    R = H_combined * A;
    % Final Q is transpose == inverse of reflections 
    Q = transpose(H_combined);
end