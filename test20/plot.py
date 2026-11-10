import numpy as np
import matplotlib.pyplot as plt

# 1. 設定參數
npts = 150
mgapSq = 1.0  # 假設物理閾值質量平方為 1.0

# 2. 計算採樣數據
k_values = np.arange(npts)
phiSamples = np.pi/2 + np.pi/2 * np.cos((k_values + 0.5) * np.pi / npts)
massSamples = np.sin(phiSamples / 2)**2
massList = mgapSq / (1 - massSamples)

# 3. 開始繪圖 (包含線性軸與對數軸對比)
fig, (ax2) = plt.subplots(1, figsize=(10,8), dpi=100)

# 子圖 2：對數座標軸 (Log Scale) -> 最佳視覺化
ax2.scatter(k_values, massList, color='#d62728', s=25, alpha=0.8, edgecolors='none')
ax2.set_yscale('log')  # 將 Y 軸設為對數軸
ax2.set_title("Logarithmic Scale: $\log(m^2)$ vs $k$", fontsize=12, fontweight='bold')
ax2.set_xlabel("Index $k$", fontsize=11)
ax2.set_ylabel("Mass Variable ($m^2$) [Log Scale]", fontsize=11)
ax2.grid(True, which="both", linestyle='--', alpha=0.5) # 包含細部網格

# plt.suptitle(f"Chebyshev Mapping Scatter Plot (npts = {npts})", fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.show()
