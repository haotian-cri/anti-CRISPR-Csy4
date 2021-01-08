# anti-CRISPR-Csy4 source codes
scripts and data for "Anti-CRISPR RNAs: designing universal riboregulators with deep learning of Csy4-mediated RNA processing"

This github repository includes matlab scripts we used to analyze the data of pre-crRNA libraries of random or palindromic flanking sequences and codes we used to design anti-CRISPR RNA triggers and pre-cRNA sensors, as well as the dataset for analysis.

Python codes of Seq2DFunc (2D convolutional neural network) and 1D CNN (1D convolutional neural network) deep learning are not included in this repository, but here: https://github.com/maincover/sequenceDLlearning

System requirements: all codes were tested on Macbook Pro (Mid-2015) with MacOS. They were not tested on Windows or Linux systems.

Install guideline:

Matlab R2017a or newer version, and the following Matlab toolboxes are required: Statistics and Machine Learning Toolbox, Bioinformatics Toolbox. The typical install time of Matlab and toolboxes can be estimated according to the offical website of Mathworks: https://www.mathworks.com/.

nupack 3.1.0 or newer version: http://www.nupack.org/.

viennarna-2.4.8 or newer version, which can be downloaded via anaconda2: https://www.anaconda.com/. It normally takes < 1 hour to download the codes and dataset.

Instruction for use:

Please first read through the codes and notations, and changes the paths of desired folders.

Most of the analysis take < 1 hour. Exceptions: Feature analysis of random libraries takes multiple hours to weeks, dependent of the size of dataset, which is normal due to the time-consuming sequence and structural analysis. The dataset used in this study normally takes < 1 week to process.

Reproduction instruction: Random seeds were pre-defined in the codes for reproduction.
