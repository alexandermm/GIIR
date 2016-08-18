# GIIR
MATLAB implementation of the Gradient-based Iterative Image Reconstruction (GIIR) method for optical tomography
using adjoint differentiation. 

This imaging algorithm is used to get an image inside tissue using infrared light. It does this by simulating the motion of photons from a set of light sources as a diffusion process, and then minimizing the least squares error between the simulated and measured light intensity distribution (spatially and over time) measured by a set of receivers. The minimization is done by changing the spatial optical properties of the medium being simulated, which can be related to different tissues. This can be done quickly using adjoint differentiation, in order to find out how the least squares error is related to the optical properties. The file giirReport.pdf describes all the techniques used for GIIR in detail. 

As shown in the presentation below, adding a regularizer to the function being minimized can give better results, and its also essential for dealing with noisy data.  

References:
* Prof. Hanson's [presentation] (http://kmh-lanl.hansonhub.com/talks/duke03-00.pdf) on optical tomography reconstruction
* The original [paper] (http://kmh-lanl.hansonhub.com/publications/tmi99a.pdf) on GIIR that implemented it using the BIE software from Los Alamos
