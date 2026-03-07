function [Q,R] = myMGSQR(A)
    [n,m] = size(A);
    
    % defaults
    Q = zeros(n,m);
    R = zeros(m,m);

    % Moving column by column
    for curr_col = 1:m
        % Take current column, its norm as the r_diag
        A_col = A(:,curr_col);
        r_diag = norm(A_col);
        % Q_col will just be the column normalized
        Q_col = A_col / r_diag;

        % For all future A columns, remove projection of this Q_col vector
        for next_col = curr_col+1:1:m
            r_entry = transpose(Q_col) * A(:, next_col);
            R(curr_col, next_col) = r_entry;
            A(:, next_col) = A(:, next_col) - r_entry * Q_col;
        end

        % add Q_col to Q and r_diag to R
        Q(:, curr_col) = Q_col;
        R(curr_col, curr_col) = r_diag;
    end
end