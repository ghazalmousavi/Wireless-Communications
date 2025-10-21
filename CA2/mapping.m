function [x_symbols, gray_matrix]  = mapping(x_bit, M)
    k = log2(M);
    gray_matrix = generate_gray_matrix(k);
    [~, idx] = ismember(x_bit, gray_matrix, 'rows');
    x_symbols = idx - 1;
end