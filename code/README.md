
The code is developed with Julia v1.10.9. You should start Julia with multiple
threads to speed up the computations. You can do this by starting Julia with the
following command in the terminal:

```bash
$ julia --threads auto
```

Then, you can execute the following commands in the Julia REPL to reproduce the
numerical results.

The code is split into two parts
- `code/code.jl`: This file contains the main code for the numerical simulations.
- `visualization/visualization.jl`: This file contains the code for visualizing the results.

The splitting makes it easily possible to run the simulations on a headless server
and visualize them on a local machine.


## Example 1: Solutions including a shock and a stationary part

### 2D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_shock_and_stationary_2d(max_level = 9);
[...]
  3.690371 seconds (241.32 k allocations: 1.631 GiB, 1.12% gc time)

julia> @time save_shock_and_stationary_2d(max_level = 10);
[...]
  7.572707 seconds (43.48 M allocations: 8.465 GiB, 5.26% gc time, 0.00% compilation time)

julia> @time save_shock_and_stationary_2d(max_level = 11);
[...]
 29.188513 seconds (171.88 M allocations: 33.140 GiB, 4.50% gc time)

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_shock_and_stationary_2d(joinpath("data", "shock_and_stationary_P4estMesh2D_3_5_11.h5"))

```


### 3D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_shock_and_stationary_3d(max_level = 7);
[...]
129.610870 seconds (276.21 M allocations: 103.577 GiB, 5.51% gc time, 14.83% compilation time)

julia> @time save_shock_and_stationary_3d(max_level = 10); # roughly 18 hours with 4 threads

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_shock_and_stationary_3d(joinpath("data", "shock_and_stationary_P4estMesh3D_3_5_10.h5"))

```


### Check entropy conservation

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_check_ec_2d()
[...]
  0.979400 seconds (2.08 M allocations: 193.298 MiB, 6.47% gc time, 59.18% compilation time)

julia> @time check_ec_3d()
[...]
 11.517281 seconds (2.54 M allocations: 1.762 GiB, 1.19% gc time, 5.15% compilation time)

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_check_ec(joinpath("data", "entropy_conservation_test_P4estMesh2D.h5"))

julia> plot_check_ec(joinpath("data", "entropy_conservation_test_P4estMesh3D.h5"))

```



## Example 2: Self-Similar Expansion

### 2D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_self_similar_expansion_2d(base_level = 7, med_level = 11, max_level = 14);
[...]
 38.226044 seconds (281.79 M allocations: 53.698 GiB, 4.52% gc time, 0.11% compilation time)

julia> @time save_self_similar_expansion_2d(base_level = 8, med_level = 11, max_level = 14);
[...]
114.076531 seconds (747.95 M allocations: 142.494 GiB, 4.42% gc time)

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_self_similar_expansion_2d(joinpath("data", "self_similar_expansion_P4estMesh2D_8_11_14.h5"))

```


### 3D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_self_similar_expansion_3d(max_level = 7);
[...]
202.918851 seconds (522.96 M allocations: 234.680 GiB, 4.86% gc time, 0.02% compilation time)

julia> @time save_self_similar_expansion_3d(max_level = 9);
[...]
217.240046 seconds (567.75 M allocations: 255.106 GiB, 4.88% gc time)

julia> @time save_self_similar_expansion_3d(max_level = 11);
[...]
257.407259 seconds (629.30 M allocations: 282.834 GiB, 4.62% gc time)

julia> @time save_self_similar_expansion_3d(max_level = 13);
[...]
278.646031 seconds (692.17 M allocations: 311.355 GiB, 4.91% gc time, 0.00% compilation time)

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_self_similar_expansion_3d(joinpath("data", "self_similar_expansion_P4estMesh3D_3_5_13.h5"))

```



## Example 3: Expansion of a Spherical Bubble

### 2D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_expansion_spherical_bubble_2d(max_level = 9);
[...]
 19.438735 seconds (98.59 M allocations: 16.375 GiB, 3.30% gc time, 31.98% compilation time)

julia> @time save_expansion_spherical_bubble_2d(max_level = 10);
[...]
 50.953607 seconds (314.10 M allocations: 60.425 GiB, 4.63% gc time)

julia> @time save_expansion_spherical_bubble_2d(max_level = 11);
[...]
208.223712 seconds (1.25 G allocations: 239.832 GiB, 4.72% gc time)

julia> @time save_expansion_spherical_bubble_2d(max_level = 13); # roughly 7 hours with 4 threads

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_expansion_spherical_bubble_2d(joinpath("data", "expansion_spherical_bubble_P4estMesh2D_3_5_13.h5"))

```


### 3D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_expansion_spherical_bubble_3d(max_level = 7);
[...]
772.037618 seconds (1.48 G allocations: 669.293 GiB, 4.81% gc time, 0.20% compilation time)

julia> @time save_expansion_spherical_bubble_3d(max_level = 8); # roughly 6 hours with 4 threads

julia> @time save_expansion_spherical_bubble_3d(max_level = 9); # roughly 48 hours with 4 threads

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_expansion_spherical_bubble_3d(joinpath("data", "expansion_spherical_bubble_P4estMesh3D_3_5_9.h5"))

```


## Example 4: Collapse of a Spherical Bubble

### 2D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_collapse_spherical_bubble_2d(max_level = 9);
[...]
  6.547556 seconds (33.97 M allocations: 6.234 GiB, 3.82% gc time, 15.71% compilation time)

julia> @time save_collapse_spherical_bubble_2d(max_level = 10);
[...]
 19.895359 seconds (118.85 M allocations: 22.986 GiB, 4.92% gc time)

julia> @time save_collapse_spherical_bubble_2d(max_level = 11);
[...]
 76.389375 seconds (464.95 M allocations: 89.556 GiB, 4.61% gc time)

julia> @time save_collapse_spherical_bubble_2d(max_level = 12);

julia> @time save_collapse_spherical_bubble_2d(max_level = 13); # roughly 24 hours with 4 threads

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_collapse_spherical_bubble_2d(joinpath("data", "collapse_spherical_bubble_P4estMesh2D_3_5_13.h5"))

```


### 3D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_collapse_spherical_bubble_3d(max_level = 7);

julia> @time save_collapse_spherical_bubble_3d(max_level = 10); # roughly 50 hours with 4 threads

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_collapse_spherical_bubble_3d(joinpath("data", "collapse_spherical_bubble_P4estMesh3D_3_5_10.h5"))

```



## Example 5: Initially Sine–Shaped Radial Velocity

### 2D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_sine_velocity_2d(max_level = 9);
[...]
 31.117678 seconds (162.98 M allocations: 26.791 GiB, 3.76% gc time, 32.61% compilation time)

julia> @time save_sine_velocity_2d(max_level = 10);
[...]
 83.245538 seconds (522.54 M allocations: 100.495 GiB, 4.45% gc time, 0.00% compilation time)

julia> @time save_sine_velocity_2d(max_level = 11);
[...]
362.589236 seconds (2.15 G allocations: 408.878 GiB, 4.77% gc time, 2.67% compilation time)

julia> @time save_sine_velocity_2d(max_level = 13); # roughly 10 hours with 4 threads

julia> @time save_sine_velocity_2d(max_level = 14); # roughly 45 hours with conflicting load

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_sine_velocity_2d(joinpath("data", "sine_velocity_P4estMesh2D_3_5_14.h5"))

```


### 3D version

First, you can run the numerical simulations as follows:

```julia
julia> include("code/code.jl")

julia> @time save_sine_velocity_3d(max_level = 7);
[...]
373.252817 seconds (711.87 M allocations: 319.969 GiB, 3.52% gc time, 0.32% compilation time)

julia> @time save_sine_velocity_3d(max_level = 8); # roughly 31 hours with conflicting load

julia> @time save_sine_velocity_3d(max_level = 9); # roughly a few days with conflicting load

```

Afterward, you can visualize the results using the following commands:

```julia
julia> include("visualization/visualization.jl")

julia> plot_sine_velocity_3d(joinpath("data", "sine_velocity_P4estMesh3D_3_5_9.h5"))

```
