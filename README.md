# GL-LRSS
**G**raph **L**earning based on **L**ow **R**ank and **S**patiotemporal **S**moothness

This algorithm estimates a graph Laplacian $\widehat{\mathbf{L}}$ and denoised matrix samples $\widehat{\mathbf{X}}$ given graph signal samples $\mathbf{X}\in\mathbb{R}^{n\times T}$.

The signal is considered to be low-rank and be smooth with the measure of the ground truth graph Laplacian $\mathbf{L}\in\mathbb{R}^{n\times n}$.

For more details please refer to the paper mentioned below.
## Keywords
- MATLAB
- graph learning
- spatiotemporal signal
- graph signal
- low rank
- spatiotemporal smoothness

# Dependencies
Please ensure to add the dependencies to the path of MATLAB before running the demo.

1. [MATLAB Odds-and-ends Supplimentaries](https://github.com/Mizera-Mondo/matlab-one-supp)
2. [nearestSPD - MATLAB File Exchange](https://ww2.mathworks.cn/matlabcentral/fileexchange/42885-nearestspd?s_tid=srchtitle)


# TODO
~~1. Refactor GL_LRSS.m with standalone low-rank decomposition.~~

# Notes and Known Problems
1. The choosing of $\alpha$ and $\beta$ should not be static, cuz the term $\mathrm{tr}\left\{\mathcal{D}(\mathbf{X})'\mathbf{L}\mathcal{D}(\mathbf{X})\right\}$ will inflate as the signal length grows.
2. The time consumed by Low-rank decomposition explodes as the length of signal grows.
 
# Reference
> Liu, Yueliang, et al. "Graph learning for spatiotemporal signals with long-and short-term characterization." IEEE Transactions on Signal and Information Processing over Networks 6 (2020): 699-713.


