import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc


def mod_demod_4pam(N, E_b, N_0):
    # Define 4-PAM symbols
    symbols = np.array([-3*np.sqrt(E_b), -np.sqrt(E_b),
                       np.sqrt(E_b), 3*np.sqrt(E_b)])

    # Generate random bit sequences (2-bit groups mapped to {0,1,2,3})
    x = np.random.randint(0, 4, N)

    # Modulation: Map bits to symbols
    s = symbols[x]

    # Add Gaussian noise
    n = np.sqrt(N_0 / 2) * np.random.randn(N)
    r = s + n

    # Demodulation: Find the closest symbol
    x_hat = np.array([symbols[np.argmin(abs(r_i - symbols))] for r_i in r])

    # Convert received symbols back to indices
    x_hat = np.array([np.where(symbols == s_i)[0][0] for s_i in x_hat])

    # Count bit errors
    # Each symbol represents 2 bits
    error_count = np.sum(x != x_hat) * np.log2(4)

    return error_count


def four_ary_pam():
    # Parameters
    N = 100  # Number of bits per packet
    N_0 = 1  # Noise variance
    no_of_packets = 1000  # Number of packets
    Max_Eb = 8  # Maximum Eb/N0 ratio in dB

    # Convert Eb/N0 from dB to linear scale
    Eb = 10**(np.arange(1, Max_Eb+1) / 10)

    bit_error_rate = np.zeros(Max_Eb)
    theoretical_ber = np.zeros(Max_Eb)

    for i in range(Max_Eb):
        total_error_count = 0

        # Process each packet
        for _ in range(no_of_packets):
            error_count = mod_demod_4pam(N, Eb[i], N_0)
            total_error_count += error_count

        # Compute simulated BER
        bit_error_rate[i] = (1/8) * (3 * erfc(np.sqrt(Eb[i] / N_0)) +
                                     2 * erfc(3 * np.sqrt(Eb[i] / N_0)) -
                                     erfc(5 * np.sqrt(Eb[i] / N_0)))

        # Compute theoretical BER
        theoretical_ber[i] = (3/8) * erfc(np.sqrt(Eb[i] / N_0))

    # Plot results
    plt.figure()
    plt.semilogy(10 * np.log10(Eb), bit_error_rate, 'r*', label="Simulation")
    plt.semilogy(10 * np.log10(Eb), theoretical_ber, 'b-', label="Theory")
    plt.grid()
    plt.xlabel("E_b/N_o in dB")
    plt.ylabel("Bit Error Probability")
    plt.legend()
    plt.title("4-Ary PAM BER Simulation")
    plt.show()


four_ary_pam()
