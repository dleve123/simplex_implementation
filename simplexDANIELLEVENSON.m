function [ xsol, optimalobjective, basisfinal ] = simplexDANIELLEVENSON(A, b, c, BAS)
%simplexDANIELLEVENSON Solves LPs using the Simplex Method
%   Custom simplex method implementation for HW 5


    % utility methods
    num_of_rows = @(matrix) size(matrix, 1);
    num_of_cols = @(matrix) size(matrix, 2);

    % create simplex pretableau

    create_inner_pretableau = @(c,A) [-c'; A];
    create_pretableau = @(inner_tableau,b) ...
        [[1 zeros(1, num_of_rows(b))]' inner_tableau [0; b]];

    inner_pretableau = create_inner_pretableau(-c, A);
    pretableau = create_pretableau(inner_pretableau, b);

    current_tableau = pretableau;

    for i = 1:numel(BAS) % for every basis column

        basis_column_index = BAS(i);

        basis_column = current_tableau(2:end, basis_column_index+1); % add 1 because we want to index the inner tableau, not the 0 vector at the rhs of it.
        current_b = current_tableau(2:end, end);

        min_ratio = inf;
        pivot_index = -1;
        unbounded = true;

        % for every element in basis column
         for j = 1:num_of_rows(basis_column)
            lhs_el = basis_column(j);
            rhs_el = current_b(j);

            if (lhs_el == 0)
                fprintf('Degeneracy detected!\n');
                pivot_index = j;
                % Break out of loop???
            elseif (lhs_el < 0)
                fprintf('coefficient negative - skipping!\n');
                % TODO: Skip loop iteration!!
            else
                unbounded = false;

                ratio = lhs_el / rhs_el;
                if ratio < min_ratio
                    pivot_index = j;
                end
            end 

         end
         
         
         pivot_index = pivot_index + 1; % convert inner tableau index to full tableau index.
         unbounded
         
         % make pivoting element equal to 1;
         
         pivot_to_1 = eye(num_of_rows(current_tableau));
         pivot_to_1(pivot_index, pivot_index) = 1./current_tableau(pivot_index,basis_column_index+1);
         
         current_tableau = pivot_to_1 * current_tableau
         
    end
end


