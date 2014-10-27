function [xsol, optimalobjective, basisfinal] = simplextwodlevens1(A, b, c)
%   Solves an LP using the Simplex Method starting from a Two Phase
%   initialization.

    % utility methods

    num_of_rows = @(matrix) size(matrix, 1);
    num_of_cols = @(matrix) size(matrix, 2);
    
    num_rows_of_A = num_of_rows(A);
    identity_for_slack = eye(num_rows_of_A);
    zero_buffer = zeros(1,num_of_cols(A));
    one_buffer = ones(1,num_of_cols(identity_for_slack));
    
    cost_vector = [zero_buffer one_buffer]';
    padded_A = [A identity_for_slack];
    
    basis_size = num_of_cols(identity_for_slack);
    last_basis_col_index = num_of_cols(padded_A);
    first_basis_col_index = last_basis_col_index - basis_size + 1;
    
    basis_indices = [first_basis_col_index:last_basis_col_index];
    
    [foo, bar, two_phase_bases] = simplexdlevens1(padded_A, b, cost_vector, basis_indices);
    fprintf('First Phase Complete! Starting at the following basis:\n');
    two_phase_bases
    
    [xsol, optimalobjective, basisfinal] = simplexdlevens1(A, b, c, two_phase_bases);

end

