A MATLAB implementation of an agent-based model of mechanically dominated vascular development.

-----------
Basic usage
-----------

Running main.m will run a sample simulation, save the output, and produce a plot of a grown network. All parameters for domain growth in an annulus can be set in the params structure within main.m. Helper functions for plotting can be found in the relevant subdirectory and are easily modified.

For full descriptions of the parameters and the context of the model, please refer to the included comments and the original publication (in preparation).

-----------
shuffle
-----------

For speed, this implementation makes use of 'shuffle', available at https://www.mathworks.com/matlabcentral/fileexchange/27076-shuffle. Precompiled mex files are included in this implementation for Windows (64 bit) and MacOS, which should work out-of-the-box (tested for MacOS). 

If you receive an error from shuffle about a missing mex file, follow the compilation guidance in helpers/shuffle/shuffle.c, which provides instructions for basic compilation. For most users, navigating to that directory and running 'mex -O shuffle.c' from within MATLAB will do the trick. Linux users may need to run 'mex -O CFLAGS="\$CFLAGS -std=c99" Shuffle.c' instead.

Benjamin Walker & Adriana Dawes
