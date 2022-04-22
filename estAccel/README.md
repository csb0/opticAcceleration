## estAccel (ffs) #################################################
A method to estimate optic acceleration from videos.

## Authors ########################################################
Mengjian Hua &lt;<mh5113@nyu.edu>&gt; and Charlie S. Burlingham &lt;<cs.burlingham@gmail.com>&gt;

## Usage ##########################################################
You must first compile MatlabPyrTools by running compilePyrTools.m in the subdirectory of MEX files called MatlabPyrTools. 

The main function is estAccel.m. To use it, input the video name/path and number of frames for which you want to compute optic acceleration beyond frame 1. Use at least 30 frames or so, as the flow estimator takes 11 consecutive frames by default, and the derivative filters require 5 consecutive flow fields. The code will return an estimate of the optic acceleration field.

## Citing #########################################################
If you use estAccel (ffs), please reference the following:
*	Burlingham, C., Hua, M., Xu, O., Bonnen, K., & Heeger, D. (2022). Heading perception and the structure of the optic acceleration field. *arXiv*.

## References #####################################################
* Farnebäck, G. (2000). Fast and Accurate Motion Estimation Using Orientation Tensors and Parametric Motion Models. ICPR.
* Farid, H., & Simoncelli, E.P. (2004). Differentiation of Discrete Multi-Dimensional Signals.
* Raudies, F., & Neumann, H. (2012). A review and evaluation of methods estimating ego-motion. Comput. Vis. Image Underst., 116, 606-633.
* matlabPyrTools. Eero Simoncelli, Laboratory for Computational Vision, HHMI / NYU (2022). https://github.com/LabForComputationalVision/matlabPyrTools

## License ########################################################
This program is free software, licensed under the GNU Affero General Public License. See attached LICENSE.txt for full terms.
