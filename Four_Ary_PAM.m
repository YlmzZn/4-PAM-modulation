function Four_Ary_PAM

    clear all; close all; clc;

    % Parameters-----------------------------------------------------------
    N = 100; % Number of bits processed per packet
    N_0 = 1; % Noise variance
    no_of_packets = 1000; % Number of packets to be processed
    Max_Eb = 8; % Maximum Eb/N0 ratio in dB

    % Initialize the bit error rate array----------------------------------
    bit_error_rate = zeros(1, Max_Eb); 
    theoretical_ber = zeros(1, Max_Eb); 
    
    % Convert Eb/N0 values from dB to linear scale-------------------------
    for i = 1:Max_Eb
        Eb(i) = 10^(i/10);
    end

    % Loop over each Eb/N0 value-------------------------------------------
    for i = 1:Max_Eb
        total_error_count = 0;

        % Process each packet----------------------------------------------
        for j = 1:no_of_packets

            % Perform modulation and demodulation for 4-PAM----------------
            error_count = Mod_Demod_4PAM(N, Eb(i), N_0);
            total_error_count = total_error_count + error_count;
        end

        % Compute simulated bit error rate---------------------------------
        bit_error_rate(i) = (1/8)*(3*erfc(sqrt(Eb(i)/(N_0))) + ...
                                  2*erfc(3*sqrt(Eb(i)/(N_0))) - ...
                                  erfc(5*sqrt(Eb(i)/(N_0))));

        % Compute the theoretical bit error rate---------------------------
        theoretical_ber(i) = (3/8) * erfc(sqrt((Eb(i)) / (N_0)));
    end

    % Plot simulation results----------------------------------------------
    clf;

    % Plot simulated BER---------------------------------------------------
    semilogy(10*log10(Eb), bit_error_rate, "r*");
    hold on;

    % Plot theoretical BER-------------------------------------------------
    semilogy(10*log10(Eb), theoretical_ber, "b-");

    grid on;
    ylabel("Bit Error Probability");
    xlabel("E_b/N_o in dB");
    legend("Simulation", "Theory");
    title("4-Ary PAM");
end

function [error_count] = Mod_Demod_4PAM(N, E_b, N_0)
    % Initialize error count-----------------------------------------------
    error_count = 0;

    % 4-PAM symbol mapping-------------------------------------------------
    symbols = [-3*sqrt(E_b), -sqrt(E_b), sqrt(E_b), 3*sqrt(E_b)];
    
    % Generate random bit sequences (2-bit groups)-------------------------
    x = randi([0, 3], N, 1); % Generate N random values from {0,1,2,3}

    % Modulate symbols-----------------------------------------------------
    s = symbols(x + 1).';

    % Generate Gaussian noise----------------------------------------------
    n = sqrt(N_0 / 2) * randn(N, 1);
    r = s + n;

    % Demodulation process-------------------------------------------------
    x_hat = zeros(N, 1); % Initialize detected symbols array
    
    for i = 1:N
        % Find the nearest symbol to the received signal-------------------
        [~, idx] = min(abs(r(i) - symbols)); % Compute minimum distance
        x_hat(i) = idx - 1; % Convert index to original symbol (0-3)
    end
    
    % Count bit errors-----------------------------------------------------
    for i = 1:N
        if x(i) ~= x_hat(i) % If detected symbol is different from transmitted
            error_count = error_count + log2(4); % Each symbol represents 2 bits
        end
    end
end
