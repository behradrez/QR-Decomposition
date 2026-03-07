%% Make Krylov
disp("Making krylov subspace...");
function [krylov] = makeKrylov(n,q)

    % To make our well-conditioned A matrix, we'll use an orthogonal matrix
    % and a diagonal matrix of small values serving as the eigenvalues.
    % Multiplying the orthogonal's tranpose by the diagonals by the
    % orthogonal will give a well-conditioned A matrix.
    [Q, ~] = qr(randn(n));                  
    D = diag(1 + 0.1*rand(n,1));          
    A = Q * D * Q';                        
    
    % Random starting vector b
    b = randn(n,1);
    
    % Construct Krylov matrix K = [b, Ab, A^2b, ..., A^(q-1)b]
    krylov = zeros(n, q);
    krylov(:,1) = b;
    
    % Each column is A * previous column, so the next multiple of A
    for j = 2:q
        krylov(:,j) = A * krylov(:,j-1);
    end

    disp(['Estimated cond(A): ', num2str(cond(A))])
    disp(['Size of K: ', mat2str(size(krylov))])
    disp(['Estimated cond(krylov): ', num2str(cond(krylov))])
end
n = 5000;
q = 400;
krylov = makeKrylov(n,q);

%% Effect of Columns on Condition Number
disp("Checking Column effect on condition number...");
function [vals] = checkColumnEffect(krylovMat)
    
    [~,c] = size(krylovMat);
    vals = zeros(1,c+1);

    for numColsRemoved = 0:c
        subKry = krylovMat(:,1:c-numColsRemoved);
        condition = cond(subKry);
        vals(numColsRemoved+1) = condition;
    end
    vals = fliplr(vals);

end

vals = checkColumnEffect(krylov);

%% Plot values

[r,c] = size(vals);
figure;
plot(0:c-1, vals, 'b-', LineWidth=1.5);
set(gca, 'YScale', 'log');
xlabel('Number of Columns Present');
ylabel('Condition Number (logarithmic)');
title('Krylov Submatrix Condition Number vs Number of Columns')
grid on;

%% Check CPU cost and Accuracy of QR methods

% Evaluation function
function [time, acc, ortho_acc] = evaluateQRFunc(func, krylov, krylovNorm)
    
    
    %time = timeit(@() func(krylov));
    tic;
    [Q, R] = func(krylov);
    time = toc;

    acc = norm(krylov - Q*R, 'fro') / krylovNorm;
    ortho_acc = norm(eye(size(Q,2)) - (transpose(Q) * Q),'fro');

    fprintf("Time taken to complete: %.6e, recontruction acc: %.6e, Q orthogonality: %.6e\n",time, acc, ortho_acc);
end

krylovNorm = norm(krylov,'fro');

% GSQR
disp("Evaluating myGSQR");
evaluateQRFunc(@myGSQR, krylov, krylovNorm);

% MGSQR
disp("Evaluating myMGSQR");
evaluateQRFunc(@myMGSQR, krylov, krylovNorm);

% HHS
disp("Evaluating myHHSQR");
evaluateQRFunc(@myHHSQR, krylov, krylovNorm);

% built-in qr
disp("Evaluating built-in qr()");
evaluateQRFunc(@qr, krylov, krylovNorm);
