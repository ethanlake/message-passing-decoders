message-passing-decoders
=====
code for simulating various types of message-passing decoders for the 2D toric code and 1D repetition code. the main papers are:
- offline decoding: [arxiv:???](arxiv link) 
- real-time decoding: (not yet posted)

contents
=====
- `1d_local_simulation.py`, `2d_local_simulation.py`: error correction using local message passing fields 
- `1d_power_simulation.py`: error correction using instant powerlaw interactions  
- `simulation_plotter.py`: plots data output (as `.jld2` files, stored in `sim_data/`) from simulation codes
- `2d_history_visualizer.py`: creates animations of spin+field histories for 2D simulations
- `pair_distance_statistics.jl`: computes statistics of anyon pair separations for iid noise 
- `figure_commands.md`: commands used to reproduce figures

usage notes
=====
- Julia files are designed to be run in the REPL; changes to parameter values should be made by making direct edits to the relevant lines.

  

