function [Q,R] = myGSQR(A)
    % Get sizes and initialize matrices to zero
    [n,m] = size(A);
    Q = zeros(n,m);
    R = zeros(m,m);
    
    % Go column by column, left to right.
    for column = 1:1:m
        % Retrieve the A column
        a_m = A(:,column);
        
        % For every row above the diagonal at the index of the column,
        % going from top to bottom, calculate the R entries. For each
        % entry, multiply it by the corresponding Q column of the same
        % index. Add this product to a "sub" vector. Subtract the sub
        % vector from the A column to get the V vector.
        sub = zeros(size(a_m));
        for row = 1:1:column-1
            % Get associated Q column and compute the r entry
            prev_q = Q(:,row);
            r_row = transpose(prev_q) * a_m;

            % Save computed r_entry
            R(row, column) = r_row;

            % Add to total subtraction
            sub = sub + r_row*prev_q;
        end

        % Compute orthogonal vector V. The diagonal R entry being its norm,
        % and Q column being its unit vector.
        v = a_m - sub;
        r_m = norm(v);
        qm = v/r_m;
        
        % Save the results
        R(column, column) = r_m;
        Q(:, column) = qm;
    end
end