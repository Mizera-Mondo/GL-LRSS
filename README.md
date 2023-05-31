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

# TODO
1. Refactor GL_LRSS.m with standalone low-rank decomposition.

# Reference
> Liu, Yueliang, et al. "Graph learning for spatiotemporal signals with long-and short-term characterization." IEEE Transactions on Signal and Information Processing over Networks 6 (2020): 699-713.

# Code Reference
This implementation of GL-LRSS includes a copy of nearestSPD.m written by John D'Errico. 
[nearestSPD - MATLAB File Exchange](https://ww2.mathworks.cn/matlabcentral/fileexchange/42885-nearestspd?s_tid=srchtitle)
