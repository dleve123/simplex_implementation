function [xsol, optimalobjective, basisfinal] = simplexdlevens1(A, b, c, BAS)
%simplexdlevens1 Solves LPs using the Simplex Method
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
    
    
    while (roundn(max_cost_coeff, -4) > 0) % While the largest cost vector is still possitive        
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
                if (roundn(current_b(i), -4) == 0)
                    fprintf('Degeneracy detected in below tableau!\n');
                end
                unbounded = false;
             else
                 for j = 1:num_of_rows(basis_column)
                    lhs_el = basis_column(j);
                    rhs_el = current_b(j);

                    
                    if (lhs_el > 0)
                        if (roundn(rhs_el, -4) == 0)
                            fprintf('Degeneracy detected in below tableau!\n');
                        end

                        unbounded = false;
                        ratio = rhs_el / lhs_el;
                        if ratio < min_ratio
                            pivot_index = j + 1;  % convert inner tableau index to full tableau index.
                            min_ratio = ratio;
                        end
                    end
                 end
             end

             if unbounded
                 error('Problem unbounded! - Exiting!');
             end

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
        
        current_cost_vector = current_tableau(1,2:end-1);
        max_cost_coeff = max(current_cost_vector);
        next_basis_id = find(current_cost_vector == max_cost_coeff);
        
        next_basis_ids = [next_basis_id];

        is_first_pass = false;
    end

    fprintf('Optimal BFS and objective function value found!\n');
    
    current_b = current_tableau(2:end,end);
    
    basisfinal = zeros(1,num_of_rows(b));
    
    for r = 1:num_of_rows(b)
        basisfinal(r) = find(roundn(current_tableau(r+1,2:end-1), -4) == 1);
    end
    
    xsol = zeros(num_of_rows(c),1);
    xsol(basisfinal) = current_b;
    
    optimalobjective = current_tableau(1,end);
end
