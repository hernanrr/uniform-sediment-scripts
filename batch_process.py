
#!/usr/bin/env python
import profile_reader as pr
import profile_plotter as pp
import sidewall_correction as sw
import global_param_plotter as gp

# Read the profile data from the raw files
pr.main()
# Plot the profiles
pp.main()
# Compute the sidewall correction
#sw.main()
# Plot the global parameters
#gp.main()



