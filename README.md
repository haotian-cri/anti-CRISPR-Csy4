# anti-CRISPR-Csy4
scripts and data for "Anti-CRISPR RNAs: designing universal riboregulators with deep learning of Csy4-mediated RNA processing"

This github repository includes matlab codes we used to analyze the data of pre-crRNA libraries of random or palindromic flanking sequences and codes we used to design anti-CRISPR RNAs, as well as the dataset for analysis.

Python codes of Seq2DFunc and 1D CNN deep learning are not included in this repository, but here: https://github.com/maincover/sequenceDLlearning

System requirements: all codes were tested on Macbook Pro (Mid-2015) with MacOS. They were not tested on Windows or Linux systems.

Install guideline:

Matlab R2017a or newer version, and the following Matlab toolboxes are required: Statistics and Machine Learning Toolbox, Bioinformatics Toolbox. The typical install time of Matlab and toolboxes can be estimated according to the offical website of Mathworks: https://www.mathworks.com/.

nupack 3.1.0 or newer version: http://www.nupack.org/.

viennarna-2.4.8 or newer version, which can be downloaded via anaconda2: https://www.anaconda.com/. It normally takes < 1 hour to download the codes and dataset.

Instruction for use:

Please first read through the codes and notations, and changes the paths of desired folders.

Most of the analysis take < 1 hour. Exceptions: Analysis of random libraries takes multiple hours to weeks, dependent of the size of dataset. It is normal, because sequence and structural analyses are time-comsuming.

Reproduction instruction: Random seeds were pre-defined in the codes for reproduction.
