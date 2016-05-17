# PyMOL plugin - color_bfactor

This plugin is under LGPL 3.0. If you do not know what it means, you shouldn't be modifying it.

This plugin allows you to colorize multiple objects based on residue conservation using clustalO alignment (you must install clustalO by yourself, but usually it is available in your repository. 
You also have to download the zero_residues PyMOL plugin and put it on your pymol startup folder. 


USAGE --

As easy as it gets, select the objects in PyMOL. Active selection will be used to align and colorize. 

If you select only two objects, the plugin will also calculate RMSD of the sidechains (based on sidechain centeroid position) and add it to bfactor, so you can analyze it visually. 


It may contain errors, but not the ones which screw your molecules. The errors will probably be related to bfactor modification, but I'll get to it eventually.

You can put the script in your pymol startup folder (usually, $HOME/.pymol/startup/)



Gustavo Machado Alvares de Lima - 2016


COMPATIBILITY --

I tried it using PyMOL 1.7.2.1. It should work on other versions, but not all as the atom selection sintaxe is important. Again, it is related to the bfactor code. If you have and old PyMOL, choose the colorize only script.

It works under Linux, probably under MAC OS X (with a few modifications maybe, I'll test it someday) and should not work in Windows, but with this Windows bash stuff going on, who knows =D.



Feel free to copy, edit and redistribute. Keep the credits.
