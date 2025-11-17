% Driver File - run file once and results for every scenario will load 

% Tests the LUdense routine under:
%   (1) No Pivoting
%   (2) Partial Pivoting
%   (3) Complete Pivoting

%clear previous console
clear; clc; close all;

% Testing Parameters
nvals = 5:5:50;          % Matrix sizes
ntests = 5;              % Number of repeated tests per size
routinenums = [1, 2, 3]; 
Mrdim = 50; Mcdim = 50;  
rng(1);                 

fprintf('LUdense Structured Tests (Empirical Task Set 1)\n');

% Storage Initialization
n_sizes = length(nvals);
n_routines = length(routinenums);
n_tasks = 8;

task1_mean_err   = zeros(n_tasks, n_routines, n_sizes);
task1_max_err    = zeros(n_tasks, n_routines, n_sizes);
task1_mean_gamma = zeros(n_tasks, n_routines, n_sizes);
task1_max_gamma  = zeros(n_tasks, n_routines, n_sizes);
task1_mean_resid = zeros(n_tasks, n_routines, n_sizes);
task1_max_resid  = zeros(n_tasks, n_routines, n_sizes);
task1_condnum    = zeros(n_tasks, 1, n_sizes); % store condition numbers

for task = 1:8
    switch task
        case 1, label = 'Task 1.1: Diagonal';
        case 2, label = 'Task 1.2: Antidiagonal';
        case 3, label = 'Task 1.3: Diagonal + Antidiagonal';
        case 4, label = 'Task 1.4: Unit Lower Triangular';
        case 5, label = 'Task 1.5: Lower Triangular';
        case 6, label = 'Task 1.6: Tridiagonal (Diagonally Dominant)';
        case 7, label = 'Task 1.7: Growth Factor Matrix';
        case 8, label = 'Task 1.8: Symmetric Positive Definite Matrix (A = L~L~^T)';
    end

    % Table to capture resulting metrics
    fprintf('\n%s\n', label);
    fprintf(' n  | cond(A)    | piv | mean_err_fac | max_err_fac | mean_gamma | max_gamma | mean_resid | max_resid\n');
    fprintf('----|-------------|-----|--------------|-------------|------------|-----------|------------|----------\n');

    for idx_n = 1:n_sizes
        n = nvals(idx_n);

        % Construct Atrue and b
        switch task
            case 1 % Diagonal
                Atrue = diag(1:n);
                b = (1:n)';

            case 2 % Antidiagonal
                Atrue = zeros(n);
                for i=1:n
                    Atrue(i,n+1-i)=i;
                end
                b = (1:n)';

            case 3 % Diagonal + Antidiagonal
                D = diag(1:n);
                T = zeros(n);
                for i=1:n
                    T(i,n+1-i)=i;
                end
                Atrue = D + T;
                b = (1:n)';

            case 4 % Unit Lower Triangular
                Ltrue = eye(n);
                for i=2:n 
                    for j=1:i-1 
                        Ltrue(i,j)=2*rand-1;
                    end 
                end
               Atrue = Ltrue; 
               b = Atrue * ones(n,1);

            case 5 % Lower Triangular diag>0 off>1
                Atrue = zeros(n);
                for i=1:n
                    Atrue(i,i)=2;
                end

                for i=2:n
                    for j=1:i-1
                        Atrue(i,j)=i-j+2; 
                    end
                end
                b = (1:n)';

            case 6 % Tridiagonal
                Atrue = zeros(n);
                for i=1:n
                    Atrue(i,i)=4;
                    if i<n
                        Atrue(i,i+1)=1;
                    end
                    if i>1
                        Atrue(i,i-1)=1;
                    end
                end
                b = (1:n)';

            case 7 % Growth Factor Test
                Atrue = zeros(n);
                for i=1:n
                    Atrue(i,i)=1; 
                    Atrue(i,n)=1;
                    for j=1:i-1
                        Atrue(i,j)=-1;
                    end
                end
                b = (1:n)';

            case 8 % Special Matrix 
                Ltilde = eye(n);
                for i=1:n
                    Ltilde(i,i)=i;
                    for j=1:i-1
                        Ltilde(i,j)=0.5*(i-j);
                    end
                end
                Atrue = Ltilde * Ltilde';
                b = (1:n)';
        end

        % Compute condition number once per Atrue
        condA = cond(Atrue);
        task1_condnum(task, 1, idx_n) = condA;

        % Run all pivoting methods
        for idx_r = 1:n_routines
            routinenum = routinenums(idx_r);
            [mean_err, max_err, mean_gamma, max_gamma, mean_resid, max_resid] = run_LU_test(Atrue, b, routinenum, n, ntests, Mrdim, Mcdim);

            % Store results
            task1_mean_err(task, idx_r, idx_n)   = mean_err;
            task1_max_err(task,  idx_r, idx_n)   = max_err;
            task1_mean_gamma(task, idx_r, idx_n) = mean_gamma;
            task1_max_gamma(task,  idx_r, idx_n) = max_gamma;
            task1_mean_resid(task, idx_r, idx_n) = mean_resid;
            task1_max_resid(task,  idx_r, idx_n) = max_resid;

            fprintf('%4d | %1.2e |  %d  |  %10.2e |  %10.2e |  %8.2f |  %8.2f |  %10.2e |  %10.2e\n',n, condA, routinenum, mean_err, max_err, mean_gamma, max_gamma, mean_resid, max_resid);
        end
    end
end

function [mean_err, max_err, mean_gamma, max_gamma, mean_resid, max_resid] = run_LU_test(Atrue, b, routinenum, n, ntests, Mrdim, Mcdim)

    err_fac_vals = zeros(ntests,1);
    gamma_vals   = zeros(ntests,1);
    resid_vals   = zeros(ntests,1);

    for itest = 1:ntests
        M = zeros(Mrdim,Mcdim);
        M(1:n,1:n) = Atrue;
        Prow = (1:n)'; Pcol = (1:n)';

        [Mout, Prow_out, Pcol_out] = LUdense(M, n, routinenum, Prow, Pcol);

        % Extract L and U 
        L = tril(Mout(1:n,1:n), -1) + eye(n);
        U = triu(Mout(1:n,1:n));

        % relative factorization error
        Aperm = Atrue(Prow_out, Pcol_out);
        err_fac_vals(itest) = norm(Aperm - L*U, 2) / norm(Atrue, 2);

        % growth factor
        gamma_vals(itest) = norm(abs(L)*abs(U), 2) / norm(Atrue, 2);

        % Solve Ax=b
        f = b(Prow_out);
        y = Lvsolve_row(Mout, f, zeros(n,1), n);
        z = Uvsolve_row(Mout, y, zeros(n,1), n);
        xcomp = zeros(n,1);
        xcomp(Pcol_out) = z;

        % relative residual
        resid_vals(itest) = norm(b - Atrue*xcomp, 2) / norm(b, 2);
    end

    % Mean and max
    mean_err = mean(err_fac_vals);  max_err = max(err_fac_vals);
    mean_gamma = mean(gamma_vals);  max_gamma = max(gamma_vals);
    mean_resid = mean(resid_vals);  max_resid = max(resid_vals);
end


