function gray_matrix = generate_gray_matrix(k)

    if k == 1
        gray_matrix = [0; 1];
    else
        pre_gray = generate_gray_matrix(k - 1);
        n = size(pre_gray, 1);
        
        gray_matrix = zeros(2 * n, k);
        
        gray_matrix(1:n, :) = [zeros(n,1), pre_gray];
        
        for i = 1:n
            gray_matrix(n + i, :) = [1, pre_gray(n - i + 1, :)];
        end
    end

end