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
    next_basis_ids = BAS;
    max_cost_coeff = inf;

    while (max_cost_coeff > 0) % While the largest cost vector is still possitive

        for i = 1:numel(next_basis_ids) % for every basis column
            basis_column_index = next_basis_ids(i) + 1; % add 1 because we want to index the inner tableau, not the full one.

            basis_column = current_tableau(2:end, basis_column_index);
            current_b = current_tableau(2:end, end);

            min_ratio = inf;
            pivot_index = -1;
            unbounded = true;

             % for every element in basis column
             if (is_first_pass == true)
                pivot_index = i + 1;
                if (current_b(i) == 0)
                    fprintf('Degeneracy detected!\n');
                end
                if (max(basis_column) > 0)
                    unbounded = false;
                end
             else
                 for j = 1:num_of_rows(basis_column)
                    lhs_el = basis_column(j);
                    rhs_el = current_b(j);

                    if (rhs_el == 0)
                        fprintf('Degeneracy detected!\n');
                        pivot_index = j;
                        % Break out of loop???
                    elseif (lhs_el < 0)
                        continue;
                    else
                        unbounded = false;
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
        
        % TODO: Implement finding next basis!
        current_cost_vector = current_tableau(1,2:end-1);
        max_cost_coeff = max(current_cost_vector);
        next_basis_id = find(current_cost_vector == max_cost_coeff);
        
        next_basis_ids = [next_basis_id];

        is_first_pass = false;
    end

    fprintf('Optimal BFS and objective function value found!\n');

    xsol = zeros(num_of_cols(inner_pretableau),1);
    current_b = current_tableau(2:end,end);
    
    basisfinal = find(0 == current_tableau(1,2:end-1)); % Perhaps this isn't rigorous enough of a check??
    
    xsol(basisfinal) = current_b;
    optimalobjective = current_tableau(1,end);
end
