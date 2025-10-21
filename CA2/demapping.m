function x_hat = demapping(xr, constellation_point)      
    num_symbols = length(xr);
    x_hat = zeros(1, num_symbols);
    for k = 1:num_symbols
        [~, idx] = min(abs(xr(k) - constellation_point));
        x_hat(:, k) = idx - 1; 
        
    end
end