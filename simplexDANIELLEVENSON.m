function [ xsol, optimalobjective, basisfinal ] = simplexDANIELLEVENSON(A, b, c, BAS)
%simplexDANIELLEVENSON Solves LPs using the Simplex Method
%   Custom simplex method implementation for HW 5

    % utility methods

    num_of_rows = @(matrix) size(matrix, 1);
    num_of_cols = @(matrix) size(matrix, 2);

    % create simplex pretableau

    inner_pretableau = [-c'; A];
    pretableau = [[1 zeros(1, num_of_rows(b))]' inner_pretableau [0; b]]

    current_tableau = pretableau;
    is_first_pass = true;

    % Round necessary to handle weird MATLAB numeric edge cases.
    while (round(max(current_tableau(1,2:end-1))) > 0) % While the largest cost vector is still possitive
        % TODO: Implement finding next basis!
        current_basis = BAS; % STUB CODE!!!

        for i = 1:numel(current_basis) % for every basis column
            basis_column_index = BAS(i) + 1; % add 1 because we want to index the inner tableau, not the 0 vector at the rhs of it.

            basis_column = current_tableau(2:end, basis_column_index);
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
                    % TODO (Tweak) : Skip loop iteration!!
                else
                    unbounded = false;
                    if (is_first_pass == true)
                        pivot_index = i + 1
                    else
                        ratio = rhs_el / lhs_el;
                        if ratio < min_ratio
                            pivot_index = j + 1;  % convert inner tableau index to full tableau index.
                            min_ratio = ratio;
                        end
                    end
                end

             end

             % TODO: throw exception/abort program because unbounded!
             unbounded

             % make pivoting element equal to 1;

             pivot_to_1 = eye(num_of_rows(current_tableau));
             pivot_to_1(pivot_index, pivot_index) = 1./current_tableau(pivot_index,basis_column_index);

             current_tableau = pivot_to_1 * current_tableau;

             for k=1:num_of_rows(current_tableau)
                 if (k == pivot_index)
                     continue
                 end

                 zeroify = eye(num_of_rows(current_tableau));
                 zeroify(k, pivot_index) = -current_tableau(k, basis_column_index);

                 current_tableau = zeroify * current_tableau;
             end

             current_tableau
        end
        is_first_pass = false;
    end

    fprintf('Optimal BFS and objective function value found!\n');

    xsol = zeros(num_of_cols(inner_pretableau),1);
    current_b = current_tableau(2:end,end);

    xsol(current_basis) = current_b;
    optimalobjective = current_tableau(1,end);
    basisfinal = current_basis;
end
