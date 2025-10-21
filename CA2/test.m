clc; clear; close all
rng(1)

%% Parameters:

snr_dB = 2;
ifft_length =  1024;
nc = 400; 
k = 1e7;
M = 4;
sfc = ceil(2 ^ 13 / nc); 
n_frame = ceil((k / 2) / (nc * sfc));

qpsk_constellation_point = [1; +1j; -1; -1j];

%% Transmitter

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

v9 = data(:, 11, 1);

% Block 5: DPSK Modulation:
for frame = 1 : n_frame
    for row = 2 : sfc + 1 
        data(row, : , frame) = mod((data(row, : , frame) + data(row - 1, : , frame)), M);
    end
end

% 
% for frame = 1 : n_frame
%     for row = 2 : sfc + 1
%         data(row, : , frame) = mod((data(row, : , frame) - data(row - 1, : , frame)), M);
%     end
% end

v8 = data(:, 11, 1);

data = qpsk_constellation_point(data + 1);
v1 = data(:, 1, 1);

% Block 6: IFFT Bin Alloation:
z = ifft_bin_allocation(data, nc, ifft_length, sfc, n_frame);
v2 = z(:, 77, 1);

% Block 7: IFFT:
x = (ifft(z, ifft_length, 2));
v3 = x(:, 77, 1);

% Block 8: CP Addition:
y = cp_addition(x, ifft_length, sfc, n_frame);
v4 = y(:, 256 + 76, 1);

% Block 9: Parrallel to Serial:
w = parrallel_to_serial(y);
v5 = w(1, 1: 10);

% Block 10: Cascade Frames:
u = cascade_frame(w, sfc, ifft_length);
v6 = u(1, 180255: 180270);

%% AWGN Channel:
% 
h =  1;
Es = mean(abs(qpsk_constellation_point .^ 2));
N0 = Es * 10 ^ (-snr_dB / 10);
noise = (randn(size(u)) + 1j*randn(size(u))) * sqrt(N0/2);

p = u + noise;

%% Receiver:

% Block 1: Frames Detection:
[ur, v7] = frames_detection(p, sfc, ifft_length, n_frame);
v7;
% Block 2: Serial to Parrallel
wr = serial_to_parrallel(ur, sfc, ifft_length, n_frame);
v = wr(:, 256 + 76, 1);

% Block 3: CP Removal:
yr = cp_removal(wr, ifft_length);
v = yr(:, 77, 1);

% Block 4: FFT:
xr = fft(yr, [], 2);
v = xr(:, 77, 1);

% Block 5: Extract Carriers from FFT Bins:
zr = extract_carriers_from_fft_bins(xr, nc, ifft_length);
v = zr(:, 1, 1);

% Block 6: DPSK Demodulation:
x_hat = zeros(size(zr));

for frame = 1 : n_frame
    x_hat(1, :, frame) = zr(1, :, frame);
    for row = 1 : sfc + 1
        x_hat(row, :, frame) = demapping(zr(row, :, frame), qpsk_constellation_point);
    end
end
v = x_hat(:, 11 , 1);
g = x_hat;
for frame = 1 : n_frame
    for row = 2 : sfc + 1
        x_hat(row, : , frame) = mod((g(row, : , frame) - g(row - 1, : , frame)), M);
    end
end

v = x_hat(:, 11 , 1);


% Block 7: Parrallel to Serial
x_hat  = x_hat(2 : end, :, :);

% Block 8: Cascade Frames:
x_hat_symbols = reshape(x_hat, 1, []);

v = x_hat_symbols(1, 5);

x_hat_symbols = x_hat_symbols(1, 1: k /2);

% Block 9: QPSK Demodulation:
x_hat_bit = zeros(k / 2, 2);
x_hat_bit(:) = gray_matrix(x_hat_symbols + 1, :);

ber = sum(x_bits(:) ~= x_hat_bit(:))/ (k);

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

function w = parrallel_to_serial(data)
    [s, n, f] = size(data);
    w = zeros(f, s * n);
    
    for i = 1:f
        frame = data(:, :, i);
        w(i, :) = reshape(frame, [], 1)'; 
    end
end



function u = cascade_frame(data, sfc, ifft_length)
%     [f, n] = size(data);
    
    n_header = 8 * ifft_length * (sfc + 1);
    
    data_serial = reshape(data.', [], 1);
    u = zeros(1, length(data_serial) + 2 * n_header);

    u(n_header + 1 : n_header + length(data_serial)) = data_serial;

end

%%-----------------------------------------------------------------------%%
function [ur, v7] = frames_detection(u, sfc, ifft_length, n_frame)

    n_header = 8 * ifft_length * (sfc + 1);

    ur = u(1, n_header + 1 : end - n_header);
    ur = reshape(ur, [], n_frame).';
    v7 = ur(1, 1:10);

end



function wr = serial_to_parrallel(ur, sfc, ifft_length, n_frame)
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