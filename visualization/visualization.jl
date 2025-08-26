# Install packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using HDF5

using LaTeXStrings
using GLMakie

set_theme!(theme_latexfonts();
           fontsize = 20,
           linewidth = 3,
           markersize = 16,
           Lines = (cycle = Cycle([:color, :linestyle], covary = true),),
           Scatter = (cycle = Cycle([:color, :marker], covary = true),))


const DATA_DIR = joinpath(dirname(@__DIR__), "data")
if !isdir(DATA_DIR)
    mkdir(DATA_DIR)
end



###########################################################
# Example 1 of Kunik, Kolb, Müller, and Thein (2024)
function plot_shock_and_stationary_2d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v')
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end

# Same setup as Example 1 of Kunik, Kolb, Müller, and Thein (2024),
# but in 3D instead of 2D
function plot_shock_and_stationary_3d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v')
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end


# Similar setup as Example 1 of Kunik, Kolb, Müller, and Thein (2024),
# but with smaller absolute value of the velocity to check entropy
# conservation
function plot_check_ec(filename)
    t, entropy_rate = h5open(filename, "r") do file
        read(file["time"]), read(file["entropy_rate"])
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = "Entropy Rate")
    lines!(ax, t, entropy_rate)
    @show extrema(entropy_rate)

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end



###########################################################
# Example 2 of Kunik, Kolb, Müller, and Thein (2024)
function plot_self_similar_expansion_2d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v'; colorrange = (0, 2))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end

function plot_self_similar_expansion_3d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v'; colorrange = (0, 1))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end



###########################################################
# Example 3 of Kunik, Kolb, Müller, and Thein (2024)
function plot_expansion_spherical_bubble_2d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v'; colorrange = (-1, 1))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p'; colorrange = (0, 3))
    Colorbar(fig[2, 2], hm)
    xlims!(ax_p, (5.0, 5.2))
    ylims!(ax_p, (0, 0.1))

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end

function plot_expansion_spherical_bubble_3d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v')
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end


###########################################################
# Example 4 of Kunik, Kolb, Müller, and Thein (2024)
function plot_collapse_spherical_bubble_2d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v'; colorrange = (-1, 0.5))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')#; colorrange = (0, 3))
    Colorbar(fig[2, 2], hm)
    xlims!(ax_p, (1.2, 1.5))
    ylims!(ax_p, (0, 0.25))

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end

function plot_collapse_spherical_bubble_3d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v')#; colorrange = (-1, 0.5))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')#; colorrange = (0, 3))
    Colorbar(fig[2, 2], hm)
    # xlims!(ax_p, (1.2, 1.5))
    # ylims!(ax_p, (0, 0.25))

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end



###########################################################
# Example 5 of Kunik, Kolb, Müller, and Thein (2024)
function plot_sine_velocity_2d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v'; colorrange = (-1, 1))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)
    xlims!(ax_p, (0.75, 1))
    ylims!(ax_p, (0, 0.25))

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end

function plot_sine_velocity_3d(filename)
    t, r, v, p = h5open(filename, "r") do file
        read(file["time"]), read(file["radius"]), read(file["velocity"]), read(file["pressure"])
    end

    fig = Figure(size = (800, 800))

    ax_v = Axis(fig[1, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Radial Velocity")
    hm = heatmap!(ax_v, t, r, v')#; colorrange = (-1, 1))
    Colorbar(fig[1, 2], hm)

    ax_p = Axis(fig[2, 1]; xlabel = L"Time $t$", ylabel = L"Radius $r$", title = "Pressure")
    hm = heatmap!(ax_p, t, r, p')
    Colorbar(fig[2, 2], hm)
    # xlims!(ax_p, (0.75, 1))
    # ylims!(ax_p, (0, 0.25))

    ax_v = Axis(fig[1, 3]; xlabel = L"Radius $r$", title = "Radial Velocity")
    lines!(ax_v, r, v[:, end]; label = "t = $(t[end])")

    ax_p = Axis(fig[2, 3]; xlabel = L"Radius $r$", title = "Pressure")
    lines!(ax_p, r, p[:, end]; label = "t = $(t[end])")

    figname = joinpath(DATA_DIR, (filename |> basename |> splitext |> first) * ".png")
    save(figname, fig, px_per_unit = 3)
    @info "Figure saved" figname

    return fig
end
