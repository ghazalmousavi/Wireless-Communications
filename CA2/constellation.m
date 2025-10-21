function cons = constellation(M, cons_name)

    %if strcmp(cons_name, "bpsk")
    %    m = 1:M;
    %    cons = exp(1j * (2 * pi * m / M )).';
    if strcmp(cons_name, "qpsk")
        m = 1 : M;
        cons = exp(1j * (2 * pi * m / M)).';


        real_part = real(cons);
        imag_part = imag(cons);
        
        % Quadrant classification
        q1 = cons(real_part > 0 & imag_part > 0);   % Quadrant I
        q2 = cons(real_part < 0 & imag_part > 0);   % Quadrant II
        q3 = cons(real_part < 0 & imag_part < 0);   % Quadrant III
        q4 = cons(real_part > 0 & imag_part < 0);       

        cons = [q1; q2; q3; q4];
    end
    %elseif strcmp(cons_name, 'odd-even')
    %    cons = [0, 1; 1, 0];
    %end
end