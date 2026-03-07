function [Q,R] = myHHSQR(A)
    [m,n] = size(A);
    % Stores all reflection matrices, sequentially left-multiplied.
    
    % The starting matrix is simply A. In each iteration, we are going to
    % left-multiply the submatrix by the computed reflection matrix. Therefore,
    % we will store the starting matrix A in a separate variable
    % "R" and operate on that only.
    R = A;
    Q = eye(m);


    % Need to compute a reflection matrix for every diagonal
    for col_idx = 1:n
        
        % Take the specific sub-column of current column, starting at the
        % row of the loop index
        x = R(col_idx:end, col_idx);

        % Pick sign that avoids catastrophic cancellation, therefore
        % minimizing potential error when subtracting
        v = x;
        v(1) = x(1) + sign(x(1)) * norm(x);

        % Take the scale 2/(v_T * v) from lectures
        scale = 2 / (v' * v);
        
        % Apply the computed reflection directly to the submatrix of the
        % overall R. 
        R(col_idx:end,:) = R(col_idx:end,:) - (scale * v) * (v' * R(col_idx:end,:));
        
        % Similarly, accumulate the reflector into Q.
        Q(:, col_idx:end) = Q(:,col_idx:end) - (Q(:, col_idx:end) * v) * (scale*v)';

    end

    % Trim Q and R to economy form
    %Q = Q(:, 1:n);
    %R = R(1:n, :);
end
