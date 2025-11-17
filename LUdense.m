function [M, Prow, Pcol] = LUdense(M, n, routinenum, Prow, Pcol)
    
    % No Pivot Method (A = LU)
    if routinenum == 1
    % Perform LU factorization without pivoting
        for k = 1:n-1
            % Error Guard 
            if abs(M(k,k)) < 1e-12
                disp('Matrix is singular or nearly singular')
            end

            for i = k+1:n
            % Calculate the multiplier for each row
            M(i,k) = M(i,k) / M(k,k);
            % Update the matrix rows based on the multiplier
            M(i, k+1:n) = M(i, k+1:n) - M(i,k) * M(k, k+1:n);
            end
        end 

    % Note that L and U are stored in A for computational efficiency
    end 
    
    % Partial Pivoting Method (PA = LU)
    if routinenum == 2

        % initialize permutation vector 
        Prow = (1:n)'; % P = [1,2,...,n] where n is the assciated column for the e vector 

        for k = 1:n-1 
            % find max element in column k and the associated pivot row
            max_value = 0;
            associated_row = k; 
            for i = k:n 
                if abs(M(i,k)) > max_value
                    max_value = abs(M(i,k));
                    associated_row = i;
                end  
            end 

            % Error Guard 
            if max_value < 1e-12
                disp('Matrix is singular or nearly singular')
            end
            
            % Check for row swap 
            if associated_row ~= k
                % create copy of row getting replaced 
                row_copy = M(k,1:n);
                % replace with better row
                M(k,1:n) = M(associated_row,1:n);
                % replace old position of better row with the copy made 
                M(associated_row,1:n) = row_copy;
                
                % replicate process for permutation vector
                row_copy = Prow(k);
                Prow(k) = Prow(associated_row);
                Prow(associated_row) = row_copy; 
            end 

            for i = k+1:n
            % Calculate the multiplier for each row
            M(i,k) = M(i,k) / M(k,k);
            % Update the matrix rows based on the multiplier
            M(i, k+1:n) = M(i, k+1:n) - M(i,k) * M(k, k+1:n);
            end
        end 
        
    end 

    % Complete Pivoting (PAQ = LU)
    if routinenum == 3

        % initalize Permutations Vectors
            Prow = (1:n)';
            Pcol = (1:n)'; 

        for k = 1:n-1 
            max_value = 0;  
            associated_row = k; 
            associated_column = k; 
            
            % find max A(i,j) in active Matrix A(k:n,k:n) 
            for i = k:n
                for j = k:n 
                    if abs(M(i,j)) > max_value 
                        max_value = abs(M(i,j)); 
                        associated_row = i; 
                        associated_column = j;
                    end
                end
            end  

            % Error Guard 
            if max_value < 1e-12
                disp('Matrix is singular or nearly singular')
            end

            if associated_row ~= k
                % create copy of row getting replaced 
                row_copy = M(k,1:n);
                % replace with better row
                M(k,1:n) = M(associated_row,1:n);
                % replace old position of better row with the copy made 
                M(associated_row,1:n) = row_copy;
                
                % replicate process for permutation vector
                row_copy = Prow(k);
                Prow(k) = Prow(associated_row);
                Prow(associated_row) = row_copy; 
            end 

            % Similar to before swapped indexing
            if associated_column ~= k
                % Swap columns in A
                col_copy = M(1:n,k);
                M(1:n,k) = M(1:n,associated_column);
                M(1:n,associated_column) = col_copy;

                % Update permutation vector for columns
                col_copy = Pcol(k);
                Pcol(k) = Pcol(associated_column);
                Pcol(associated_column) = col_copy;
            end

            for i = k+1:n
            % Calculate the multiplier for each row
            M(i,k) = M(i,k) / M(k,k);
            % Update the matrix rows based on the multiplier
            M(i, k+1:n) = M(i, k+1:n) - M(i,k) * M(k, k+1:n);
            end
        end 
    end 
end





