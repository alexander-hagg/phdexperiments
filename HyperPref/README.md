# HyperPref

HyperPref is a cyclic co-creative process between human and computer, which triggers the user’s intuition about what is possible and is updated according to what the user wants based on their decisions. It combines quality diversity algorithms, a divergent optimization method that can produce many, diverse solutions, with variational autoencoders to both model that diversity as well as the user’s preferences, discovering the preference hypervolume within large search spaces.

I have removed everything but the bare necessities for HyperPref to run, so you will not be able to reproduce all results from the paper. If you use this code in any publication, please do cite:

Hagg, A., Asteroth, A., Bäck, T.: A Deep Dive Into Exploring the Preference Hypervolume. International Conference on Computational Creativity, ICCC 2020.

Author: Alexander Hagg
Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
email: info@haggdesign.de
Nov 2019; Last revision: 10-Sep-2020

# Demo

Run the demo.m script from the root directory of the repository. You can step through a (somewhat user unfriendly) version of HyperPref. Remember that in Matlab you can step through code blocks that are delimted by '%%' using CTRL-ENTER

1. Train a VAE on a random shape set and run QD based on the VAE's 2D latent space
2. Select shapes (sorry, you will have to manually change the selection IDs, line 30: selectionIDs = [20, 25];)
3. Retrain the VAE on a set with perturbed variations of the selected shapes and run QD again to get the final result.
