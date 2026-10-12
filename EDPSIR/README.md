<div align="center">

# 🦠 Inverse PINNs with Mollifier Layers
### Inferring spatio-temporal reaction rates in a fourth-order SIR-type model with Physics-Informed Neural Networks

<br>

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![PyTorch](https://img.shields.io/badge/PyTorch-%23EE4C2C.svg?style=for-the-badge&logo=PyTorch&logoColor=white)
![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white)
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-%230C55A5.svg?style=for-the-badge&logo=scipy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-%23ffffff.svg?style=for-the-badge&logo=Matplotlib&logoColor=black)

<br>

</div>

---

## 📌 About This Notebook

This notebook develops, end to end, a **Physics-Informed Neural Network (PINN)** framework for an **inverse problem**: recovering an unknown, spatially and temporally varying transmission rate $\lambda(t,x)$ in a **fourth-order (bi-harmonic) reaction–diffusion epidemic model** from sparse, noisy observations. Its key numerical ingredient is the **mollifier layer**, a differentiable convolution that replaces nested automatic differentiation for high-order derivatives. It is organized as follows:

| Step | Section | Description |
|:-:|:---|:---|
| 🔹 | **Inverse problems for PDEs** | Forward vs. inverse problems, **Hadamard ill-posedness**, regularized variational formulation and the PINN loss. |
| 🔹 | **The high-order derivative problem** | Why computing $\partial_x^4 u$ by four nested `autograd.grad` calls is expensive and numerically fragile. |
| 🔹 | **Mollifier layers** | Rigorous justification via **Friedrichs mollification** and the commutation of differentiation and convolution. |
| 🔹 | **Fourth-order SIR model** | Biological motivation for a bi-harmonic term (long-range dispersal) and identifiability analysis of $\lambda(t,x)$. |
| 🔹 | **PyTorch implementation** | Spectral forward solver, Gaussian/bump kernels, `conv1d` mollifier layer, two-network PINN, Adam + L-BFGS training. |
| 🔹 | **Mollifier vs. autodiff** | Microbenchmark and matched-budget training comparison between both derivative strategies. |

---

## 🧮 Mathematical Highlights

**🦠 The governing PDE.** The infected fraction $u(t,x)\in[0,1]$ obeys an extended Fisher–KPP / SIR-type equation

$$u_t + D_4\,\Delta^2 u - D_2\,\Delta u = \lambda(t,x)\,u\,(1-u),$$

where $D_4\Delta^2 u$ models long-range dispersal and damps high-frequency modes, $D_2\Delta u$ is classical Fickian diffusion, and $\lambda(t,x)$ is the **unknown effective transmission rate**. Wherever $u\approx 0$ or $u\approx 1$ the reaction term vanishes regardless of $\lambda$, making the problem **locally unidentifiable** and motivating explicit regularization.

**📐 Derivatives through convolution.** Since differentiation commutes with convolution, $\partial^\alpha(u*\varphi_\epsilon) = u*\partial^\alpha\varphi_\epsilon$, so for small $\epsilon$

$$\frac{\partial^4 u}{\partial x^4}(x) \;\approx\; \int u(x-y)\,\partial_x^4\varphi_\epsilon(y)\,dy .$$

The right-hand side needs **no differentiation of the network**: only evaluations on a stencil, contracted against the analytically known kernel derivative (for a Gaussian, $\varphi_\epsilon^{(k)}(x) = (-1)^k\epsilon^{-k}He_k(x/\epsilon)\,\varphi_\epsilon(x)$). The mollifier width $\epsilon$ acts as a **regularization parameter** governing the bias–variance trade-off.

**🧠 PINN loss.** Two networks $\hat u_\theta$ and $\hat\lambda_\phi$ (with a periodic embedding that enforces the boundary condition by construction) are trained jointly on

$$\mathcal L = w_d\,\text{MSE}_{\text{data}} + w_r\,\text{MSE}_{\text{PDE}} + w_{ic}\,\text{MSE}_{\text{IC}} + w_\lambda \,\big\|\partial_x \hat\lambda_\phi\big\|^2,$$

with $u_t$, $u_{xx}$ from autodiff and $u_{xxxx}$ from the mollifier layer.

---

## 🔬 Experiments & Results

### 🌊 Synthetic ground truth
A **Fourier spectral solver** with an implicit stiff integrator (`BDF`) generates the reference solution on a periodic domain with $D_4 = 5\cdot10^{-4}$, $D_2 = 10^{-2}$ and $\lambda(t,x) = 1 + 0.5\sin(2\pi x)\cos(\pi t)$. Only **160 sparse observations** corrupted by Gaussian noise ($\sigma = 0.02$) are given to the PINN.

### 🎯 Validating the mollifier kernel
A sweep over $\epsilon$ on an analytic test function confirms the predicted bias behavior: the relative error of the 4th-derivative estimate drops from **24%** ($\epsilon = 0.06$) to **3%** ($\epsilon = 0.02$). It also shows that the kernel must be truncated at $6\epsilon$ (rather than the usual $4\epsilon$) because of the polynomially reweighted tail of $He_4$.

### 🧩 Inverse problem reconstruction
After 2,000 Adam iterations plus L-BFGS fine-tuning:
- **State** $u(t,x)$: relative $L^2$ error ≈ **3.3%**.
- **Reaction rate** $\lambda(t,x)$: relative $L^2$ error ≈ **25%**, with the largest errors in weakly identifiable regions.
- **Mollifier-estimated** $u_{xxxx}$ vs. finite-difference reference: correlation ≈ **0.95**.

### ⚔️ Mollifier vs. pure autodiff
Computing $u_{xxxx}$ (forward + backward) on 500 collocation points takes **~0.040 s** with the mollifier layer vs. **~0.067 s** with four nested `autograd.grad` calls. In a matched 400-iteration training run, both approaches reach the **same accuracy**, with the mollifier being about **15% faster per iteration**.

---

## ⚖️ Strengths & Limitations

| ✅ Strengths | ⚠️ Limitations |
|:---|:---|
| **Avoids nested autodiff**: one forward pass on a stencil plus a fixed linear contraction | **Kernel truncation & quadrature bias** depend on $\epsilon$ and the stencil spacing $h$ |
| **Fully differentiable** w.r.t. network parameters, batchable as a `conv1d` | **Choice of $\epsilon$** is a regularization hyperparameter, fixed here by a sanity check |
| **Mesh-free inversion**: recovers $u$ and $\lambda$ jointly without a forward solver in the loop | **Stencil cost** grows with dimension ($2K+1$ evaluations per point in 1D) |
| **Interpretable regularization** via $\epsilon$ and the PDE residual | **Ill-posedness**: $\lambda$ remains poorly identified where $u \approx 0$ or $u \approx 1$ |

---

## 🛠️ Tech Stack

<div align="center">

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![PyTorch](https://img.shields.io/badge/PyTorch-%23EE4C2C.svg?style=for-the-badge&logo=PyTorch&logoColor=white)
![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white)
![LaTeX](https://img.shields.io/badge/latex-%23008080.svg?style=for-the-badge&logo=latex&logoColor=white)

</div>

**Main libraries:** `torch` · `numpy` · `scipy.integrate` (`solve_ivp`) · `scipy.special` · `scipy.interpolate` · `matplotlib`
