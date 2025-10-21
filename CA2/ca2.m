clc; clear; close all

%% Parameters:

snr_dB = 0: 2: 15;
ifft_length =  1024;
nc = 400; 
k = 1e7;
M = 4;
sfc = ceil(2 ^ 13 / nc); 
n_frame = ceil((k / 2) / (nc * sfc));

qpsk_constellation_point = [1; +1j; -1; -1j];

%% Part 1: AWGN Channel: snr = 20 

snr_dB1 = 20;
ber = ofdm(k, sfc, nc, n_frame, ifft_length, M, snr_dB1, qpsk_constellation_point);
fprintf("Bit error rate is = %.2f for SNR = %d dB\n", ber, snr_dB1);

%% Part 2: BER vs SND(dB)

ber_list = zeros(size(snr_dB));

for j = 1 : length(snr_dB)
    ber_list(j) = ofdm(k, sfc, nc, n_frame, ifft_length, M, snr_dB(j), qpsk_constellation_point);
end

figure;
semilogy(snr_dB, ber_list, 'bs-', 'LineWidth', 1.5)
title('Bit Error Rate', 'Interpreter', 'latex')
ylabel('BER', 'Interpreter', 'latex')
xlabel('SNR(dB)', 'Interpreter', 'latex')

%% Functions:

function z = ifft_bin_allocation(data, nc, ifft_length, sfc, n_frame)

    z = zeros(sfc + 1, ifft_length, n_frame);

    n = ceil((ifft_length  - 2*nc)/3);

    z(:, n + 1 : n + nc, :) = data;

    z(:, 2*n + nc + 1: 2*n + 2*nc, :) = conj(data);

end

function  y = cp_addition(data, ifft_length, sfc, n_frame)

        n = 0.25 * ifft_length;

        y = zeros(sfc + 1, n + ifft_length, n_frame);

        y(:, 1 : n, :) = data(:, end - n + 1 : end, :);

        y(:, n + 1 : end, :) = data;

end

function w = parallel_to_serial(data)
    [s, n, f] = size(data);
    w = zeros(f, s * n);
    
    for i = 1:f
        frame = data(:, :, i);
        w(i, :) = reshape(frame, [], 1)'; 
    end
end

function u = cascade_frame(data, sfc, ifft_length)    
    n_header = 8 * ifft_length * (sfc + 1);
    
    data_serial = reshape(data.', [], 1);
    u = zeros(1, length(data_serial) + 2 * n_header);

    u(n_header + 1 : n_header + length(data_serial)) = data_serial;

end

%%-----------------------------------------------------------------------%%
function ur = frames_detection(u, sfc, ifft_length, n_frame)

    n_header = 8 * ifft_length * (sfc + 1);
    ur = u(1, n_header + 1 : end - n_header);
    ur = reshape(ur, [], n_frame).';

end



function wr = serial_to_parallel(ur, sfc, ifft_length, n_frame)
    wr = reshape(ur', sfc + 1, 1.25 * ifft_length, n_frame);
    
end


function yr = cp_removal(wr, ifft_length)
    n_cp = 0.25 * ifft_length;
    yr = wr(:, n_cp + 1 : end, :);
end

function zr = extract_carriers_from_fft_bins(fft_data, nc, ifft_length)
     n = ceil( (ifft_length  - 2*nc) / 3);

     zr = fft_data(:, n + 1 : n + nc, :);

end

%%-----------------------------------------------------------------------%%

function ber = ofdm(k, sfc, nc, n_frame, ifft_length, M, snr_dB, qpsk_constellation_point)
    
    %%------------ Transmitter------------%%: 
    % Block 1: Generate Random Data:
    x_bits = reshape(randi([0, 1], k, 1), [], 2);
    
    % Block 2: QPSK Modulation:
    [x_symbols, gray_matrix]  = mapping(x_bits, M);
    
    % Block 3: Frame Divider
    data = zeros(sfc + 1, nc, n_frame);
    x_symbols_padded = [x_symbols(:); zeros(sfc*nc*n_frame - k/2, 1)];
    x_symbols_reshaped = reshape(x_symbols_padded, sfc, nc, n_frame);
    data(2 : end, :, :) = x_symbols_reshaped;
    
    % Block 4: Serial to Parrallel & Add Reference Row:
    r = ceil((k / 2 - nc*sfc*(n_frame -1))/ sfc) + nc*(n_frame -1);
    ref_row = randi([0, 3], 1, nc, n_frame);
    data(1, :, :) = ref_row;      
    data(1, r:end, end) = zeros(1, nc - r);
    
    % Block 5: DPSK Modulation:
    for frame = 1 : n_frame
        for row = 2 : sfc + 1 
            data(row, : , frame) = mod((data(row, : , frame) + data(row - 1, : , frame)), M);
        end
    end
    data = qpsk_constellation_point(data + 1);
    
    % Block 6: IFFT Bin Alloation:
    z = ifft_bin_allocation(data, nc, ifft_length, sfc, n_frame);
    
    % Block 7: IFFT:
    x = (ifft(z, ifft_length, 2));
    
    % Block 8: CP Addition:
    y = cp_addition(x, ifft_length, sfc, n_frame);
    
    % Block 9: Parallel to Serial:
    w = parallel_to_serial(y);
    
    % Block 10: Cascade Frames:
    u = cascade_frame(w, sfc, ifft_length);
    
    %%------------ Channel------------%%: 
    Es = mean(abs(u .^ 2));
    N0 = Es * 10 ^ (-snr_dB / 10);
    noise = (randn(size(u)) + 1j*randn(size(u))) * sqrt(N0/2);  
    p = u + noise;
    
    %%------------ Receiver------------%%: 
    % Block 1: Frames Detection:
    ur = frames_detection(p, sfc, ifft_length, n_frame);

    % Block 2: Serial to Parrallel
    wr = serial_to_parallel(ur, sfc, ifft_length, n_frame);
    
    % Block 3: CP Removal:
    yr = cp_removal(wr, ifft_length);
    
    % Block 4: FFT:
    xr = fft(yr, [], 2);
    
    % Block 5: Extract Carriers from FFT Bins:
    zr = extract_carriers_from_fft_bins(xr, nc, ifft_length);
    
    % Block 6: DPSK Demodulation:
    x_hat = zeros(size(zr));
    
    for frame = 1 : n_frame
        x_hat(1, :, frame) = zr(1, :, frame);
        for row = 1 : sfc + 1
            x_hat(row, :, frame) = demapping(zr(row, :, frame), qpsk_constellation_point);
        end
    end
    g = x_hat;
    for frame = 1 : n_frame
        for row = 2 : sfc + 1
            x_hat(row, : , frame) = mod((g(row, : , frame) - g(row - 1, : , frame)), M);
        end
    end    
    
    % Block 7: Parallel to Serial
    x_hat  = x_hat(2 : end, :, :);
    
    % Block 8: Cascade Frames:
    x_hat_symbols = reshape(x_hat, 1, []);    
    x_hat_symbols = x_hat_symbols(1, 1: k /2);
    
    % Block 9: QPSK Demodulation:
    x_bits_hat = zeros(k / 2, 2);
    x_bits_hat(:) = gray_matrix(x_hat_symbols + 1, :);

    ber = sum(x_bits(:) ~= x_bits_hat(:))/ (k);

end