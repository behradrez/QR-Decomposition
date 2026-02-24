function testmyHHSQR()
    n = 10;
    m = 5;
    A = randi(1000, [n,m]);
    
    [Q,R] = myHHSQR(A);
    
    recomposed = Q*R;
    tol = 1e-10;
    is_proper_decompose = abs(recomposed - A) < tol;
    if ~is_proper_decompose
        for i = 1:1:n
            for j = 1:1:m
                entry_equal = abs(recomposed(i,j) - A(i,j)) < tol;
                if ~entry_equal
                    fprintf("Entry at %d %d not equal: Q*R: %f A: %f\n", i,j,recomposed(i,j), A(i,j));
                end
            end
        end
        disp("Q*R:");
        disp(Q*R);
        disp("A");
        disp(A);
        disp("Q*R not equal A");
        assert(is_proper_decompose);
    end

    is_orthonormal = abs(transpose(Q)*Q - eye(n)) < tol;
    if ~is_orthonormal
        disp("Q not orthonormal: ");
        disp(transpose(Q)*Q);
        assert(is_orthonormal);
    end

    disp("myMGSQR Tests Passed!");
    

end