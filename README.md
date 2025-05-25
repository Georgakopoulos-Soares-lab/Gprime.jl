# GPrime v0.1
Copyright 2025 Congzhou M Sha

License is in LICENSE.txt.

1. Install Julia (version >= 1.10).
2. Install required Julia packages. In a terminal, execute `julia`. Next, press the `]` key, and execute `activate .`, followed by `add ArgParse JLD2 ProgressBars DataFrames CSV GZip` and accept the installation. 
3. Press the backspace key to exit out of package mode and exit using `exit()`.
4. Run `julia Gprime.jl -h` for usage instructions. To run in multithreaded mode for classification, run `julia -t {NUMBER OF THREADS} Gprime.jl` instead. Note that if you run multi-threaded mode, you will not get updates via the progress bar. Additionally, {NUMBER_OF_THREADS} should be in the format `{N},1`, to prevent blocking of the master thread which coordinates all other threads.
