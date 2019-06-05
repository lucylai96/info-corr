function t_matrix = trunc(matrix, num)
% truncate the data matrix
% input: matrix is the matrix you want to truncate, num is how many trials you want to keep
% output: t_matrix output matrix

t_matrix  = matrix(1:num,:);

end
