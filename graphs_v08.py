import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

x = np.linspace(-2.5, 1.5, 500)
ax1.plot(x, x**3 - 0.1*x**2 + 0.4*x + 2, label='y = x³ − 0.1x² + 0.4x + 2')
ax1.set_ylim(-8, 8)
ax1.set_title('Задание 2.1: x³ − 0.1x² + 0.4x + 2 = 0')

x = np.linspace(-1.5, 3, 500)
ax2.plot(x, (x - 1)**2, label='y = (x − 1)²')
ax2.plot(x, 0.5*np.exp(x), label='y = 0.5eˣ')
ax2.set_ylim(-0.5, 6)
ax2.set_title('Задание 2.2: (x − 1)² = 0.5eˣ')

for ax in (ax1, ax2):
    ax.axhline(0, color='black', lw=0.8)
    ax.axvline(0, color='black', lw=0.8)
    ax.grid(alpha=0.3)
    ax.legend()

plt.tight_layout()
plt.savefig('graphs_v08.png', dpi=150)
plt.show()
