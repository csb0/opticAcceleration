## Authors ########################################################
Oliver Xu &lt;<yx1797@nyu.edu>&gt;, Mengjian Hua &lt;<mh5113@nyu.edu>&gt;, and Charlie S. Burlingham &lt;<cs.burlingham@gmail.com>&gt;

## Usage ##########################################################
The main file used is function_test.m, which calls all of the requisite functions. Desired heading, fixation depth, and output directory can be adjusted through their respective variables. Fixation depth is given in terms of world coordinates, not displacement from observer. The observer is at -5 m in this coordinate system, so a 12.5 m fixation depth should be set to 7.5 m in the code, for example. function_test.m generates figures corresponding to Figures 3-6 in the paper.

Runs slowly due to very large plane size used, however, it is no longer storage intensive.

## Citing #########################################################
If you use the functions contained in this folder, please reference the following:
*	Burlingham, C., Hua, M., Xu, O., Bonnen, K., & Heeger, D. (2022). Heading perception and the structure of the optic acceleration field. *arXiv*.
* Matthis JS, Muller KS, Bonnen KL, Hayhoe MM. Retinal optic flow during natural locomotion. PLoS Comput Biol. 2022 Feb 22;18(2):e1009575. doi: 10.1371/journal.pcbi.1009575. PMID: 35192614; PMCID: PMC8896712.

## References #####################################################
* Matthis JS, Muller KS, Bonnen KL, Hayhoe MM. Retinal optic flow during natural locomotion. PLoS Comput Biol. 2022 Feb 22;18(2):e1009575. doi: 10.1371/journal.pcbi.1009575. PMID: 35192614; PMCID: PMC8896712.

## License ########################################################
This program is free software, licensed under the GNU Affero General Public License. See attached LICENSE.txt for full terms.
