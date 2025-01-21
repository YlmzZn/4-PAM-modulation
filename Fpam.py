import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

N = 10**5  # Number of symbols
alpha4pam = np.array([-3, -1, 1, 3])  # 4-PAM symbol set
Es_N0_dB = np.arange(0, 9)  # SNR values in dB

ipHat = np.zeros(N)  # To store demodulated symbols
nErr = np.zeros(len(Es_N0_dB))  # Error count array

for ii in range(len(Es_N0_dB)):
    # Random symbol generation
    ip = np.random.choice(alpha4pam, size=N)

    # Energy normalization
    s = (1 / np.sqrt(5)) * ip

    # Additive white Gaussian noise
    noise = (1 / np.sqrt(2)) * (np.random.randn(N) + 1j * np.random.randn(N))

    # Received signal
    r = s + 10**(-Es_N0_dB[ii] / 20) * noise

    # Demodulation (decision boundaries)
    r = np.real(r)  # Take only the real part

    ipHat[r < -2 / np.sqrt(5)] = -3
    ipHat[(r < 0) & (r >= -2 / np.sqrt(5))] = -1
    ipHat[(r >= 0) & (r < 2 / np.sqrt(5))] = 1
    ipHat[r >= 2 / np.sqrt(5)] = 3

    # Counting errors
    nErr[ii] = np.sum(ip != ipHat)

# Simulated BER
simBer = nErr / N

# Theoretical BER
theoryBer = 0.75 * erfc(np.sqrt(0.2 * 10**(Es_N0_dB / 10)))

plt.figure()
plt.semilogy(Es_N0_dB, theoryBer, 'b-', label='Theory')
plt.semilogy(Es_N0_dB, simBer, 'r*', label='Simulation')
plt.axis([1, 8, 1e-4, 1])
plt.grid(True, which='both')
plt.legend()
plt.xlabel('Es/No (dB)')
plt.ylabel('Symbol Error Rate')
plt.title('Symbol Error Probability Curve for 4-PAM Modulation')
plt.show()
