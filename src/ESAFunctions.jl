module ESAFunctions

# Import dependencies =========================================================
import DataFrames, CategoricalArrays
import SQLite, DBInterface, HDF5
import PlotlyJS, Plots
import Random, Distributions, Statistics, StatsBase

# Declare export ==============================================================
# Visualization functions
export plottimeseries, plotbar, plothistogram, plotcontour, plotsurface, plotheatmap, plotvolume, plotcluster, plotseries_percentile, plotbox, plotmap
# Data interface functions
export save_dftoh5, load_h5todf, loadall_h5todf, save_dftodb, load_dbtodf, list_dbtable, appendtxt
# Profile modification functions
export generatedailypattern, generatepoissonseries, synthesizeprofile, normalize
# Miscellaneous functions
export clippy, createdummydata, mergeposneg, averageprofile, equipmentloading, df_filter_collapse_plot, getcolor
# Pending retirement

# Visualization functions =====================================================
"""
    plottimeseries(df::DataFrame; xlab, ylab, title, col_time, col_variable, col_value, bstack, selectcolor, legendorientation)

Plot a line chart from `df`, a dataframe with columns `:time`, `:variable`, `:value`

# Keyword Arguments
- `col_time`, `col_variable`, `col_value` as `Symbol` : overwrite the default column names
- `bstack` as `Bool` : stack the components
- `selectcolor` : a function that returns a color given a variable name
- `legendorientation` : "h" or "l"
"""
function plottimeseries(df::DataFrames.DataFrame;
    xlab::String="Time", ylab::String="Power (kW)", title::Union{String, Missing}=missing,
    col_time=:time, col_variable=:variable, col_value=:value,
    bstack::Bool=false, selectcolor=missing,
    legendorientation="h",
    )
    # Handle when col_variable is missing.
    col_variable ∉ propertynames(df) && (df[!, col_variable] .= "")
    # Color palette
    ismissing(selectcolor) && (selectcolor = (x -> missing))
    # Plot settings
    stackgroup = (bstack ? "one" : missing)
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    for gd ∈ DataFrames.groupby(df, col_variable)
        push!(pTraces, PlotlyJS.scatter(x=gd[:, col_time], y=gd[:, col_value], name=gd[1, col_variable], mode="lines", stackgroup=stackgroup, line=PlotlyJS.PlotlyBase.attr(color=selectcolor(gd[1, col_variable]))))
    end
    showlegend = length(pTraces) > 1
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.0)", # Transparent plot BG
        paper_bgcolor="rgba(255,255,255,1.0)", # White paper BG
        title=title,
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=showlegend, legend=PlotlyJS.attr(orientation=legendorientation),
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end
# Dispatches of `plottimeseries`
"""
    plottimeseries(dt::Dict; kwargs...)

Format `dt::Dict{String, Vector{Real}}` into a `df` with a proper structure, then pass it to `plottimeseries`
"""
function plottimeseries(dt::Dict; kwargs...)
    df = DataFrames.DataFrame()
    [df = vcat(df, DataFrames.DataFrame(:time => 1:length(dt[k]), :variable => k, :value => dt[k])) for k ∈ keys(dt)]
    return plottimeseries(df; kwargs...)
end

"""
    plottimeseries(vt::AbstractVector; kwargs...)

Format a single timeseries vector into a `df` with a proper structure, then pass it to `plottimeseries`
"""
plottimeseries(vt::AbstractVector; kwargs...) = plottimeseries(DataFrames.DataFrame(:time => 1:length(vt), :variable => "", :value => vt); kwargs...)

"""
    plotbar(df::DataFrame; xlab, ylab, title, col_axis, col_variable, col_value, bstack, selectcolor, legendorientation)

Plot a bar chart from `df` which contains columns `:axis`, `:variable`, `:value`

# Keyword Arguments
- `col_axis`, `col_variable`, `col_value` as `Symbol` : overwrite the default column names
- `bstack` as `Bool` : stack the components
- `selectcolor` : a function that returns a color given a variable name
- `legendorientation` : "h" or "l"
"""
function plotbar(df::DataFrames.DataFrame;
    xlab::String="Scenario", ylab::String="", title::Union{String, Missing}=missing,
    col_axis=:axis, col_variable=:variable, col_value=:value,
    bstack::Bool=false, selectcolor=missing,
    legendorientation="h",
    )
    subfunction_xlabel(df, col_axis) = (isa(col_axis, Vector) ? [df[:, col] for col ∈ col_axis] : df[:, col_axis])
    # Color palette
    ismissing(selectcolor) && (selectcolor = (x -> missing))
    # Plot settings
    barmode = (bstack ? "stack" : missing)
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    for gd ∈ DataFrames.groupby(df, col_variable)
        push!(pTraces, PlotlyJS.bar(x=subfunction_xlabel(gd, col_axis), y=gd[:, col_value], name=gd[1, col_variable], marker=PlotlyJS.PlotlyBase.attr(color=selectcolor(gd[1, col_variable]))))
    end
    showlegend = length(pTraces) > 1
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.0)", # Transparent plot BG
        paper_bgcolor="rgba(255,255,255,1.0)", # White paper BG
        title=title,
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=showlegend, legend=PlotlyJS.attr(orientation=legendorientation),
        barmode=barmode,
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end

"""
    plothistogram(dt::Dict; xlab, ylab, title)

Plot histogram from `dt`, a dictionary of number vectors
"""
function plothistogram(dt::Dict; xlab::String="Value", ylab::String="", title::Union{String, Missing}=missing)
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    lvlOpacity = (length(dt) == 1 ? 1.00 : 0.60)
    for key ∈ keys(dt)
        push!(pTraces, PlotlyJS.histogram(x=dt[key], name=string(key), opacity=lvlOpacity))
    end
    showlegend = length(pTraces) > 1
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.0)", # Transparent plot BG
        paper_bgcolor="rgba(255,255,255,1.0)", # White paper BG
        title=title,
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=showlegend, legend=PlotlyJS.attr(orientation="h"),
        barmode="overlay",
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end

"""
    plothistogram(vt::Vector; kwargs...)

Format `vt` to a proper structure, and pass it to `plothistogram`
"""
plothistogram(vt::Vector; kwargs...) = plothistogram(Dict("" => vt); kwargs...)

"""
    plotcontour(x, y, z; xlab, ylab, title, zmin, zmax)

Create a contour plot given the range of `x` and `y`, and the matrix `z` (n_x x n_y)

# Keyword Arguments
- `zmin` and `zmax` : limits of `z` for the color gradient
"""
function plotcontour(x, y, z; xlab::String="x", ylab::String="y", title::Union{String, Missing}=missing, zmin=nothing, zmax=nothing)
    trace = PlotlyJS.contour(
        x=x, y=y, z=z, zmin=zmin, zmax=zmax,
        contours=PlotlyJS.attr(
            showlabels=true, 
            labelfont = PlotlyJS.attr(color="darkgray")
        )
    )
    layout = PlotlyJS.Layout(
        plot_bgcolor="rgba(255,255,255,0.0)", # Transparent plot BG
        paper_bgcolor="rgba(255,255,255,1.0)", # White paper BG
        title=title,
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showscale=true,
    )
    p = PlotlyJS.plot(trace, layout)
    return p
end

"""
    plotcontour(x, y, f::Function; kwargs...)

Plot a contour of a bivariante function `f` providing that its `x` and `y` ranges are given
"""
plotcontour(x, y, f::Function; kwargs...) = plotcontour(x, y, f.(x, y'); kwargs...)

"""
    plotsurface(x, y, z; xlab, ylab, zlab, title, zmin, zmax)

Create a contour plot given the range of `x` and `y`, and the matrix `z` (n_x x n_y)

# Keyword Arguments
- `zmin` and `zmax` : limits of `z` for the color gradient
"""
function plotsurface(x, y, z; xlab::String="x", ylab::String="y", zlab::String="z", title::Union{String, Missing}=missing, zmin=nothing, zmax=nothing)
    trace = PlotlyJS.surface(x=x, y=y, z=z, zmin=zmin, zmax=zmax)
    layout = PlotlyJS.Layout(
        plot_bgcolor="rgba(255,255,255,0.0)", # Transparent plot BG
        paper_bgcolor="rgba(255,255,255,1.0)", # White paper BG
        title=title,
        showscale=true,
        scene=PlotlyJS.attr(
            xaxis_title=xlab,
            xaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
                zerolinecolor="rgba(0,0,0,1.00)",
            ),
            yaxis_title=ylab,
            yaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
                zerolinecolor="rgba(0,0,0,1.00)",
            ),
            zaxis_title=zlab,
            zaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
                zerolinecolor="rgba(0,0,0,1.00)",
            ),
        ),
    )
    p = PlotlyJS.plot(trace, layout)
    return p
end

"""
    plotsurface(x, y, f::Function; kwargs...)

Plot a surface of a bivariante function `f` providing that its `x` and `y` ranges are given
"""
plotsurface(x, y, f::Function; kwargs...) = plotsurface(x, y, f.(x, y'); kwargs...)

"""
(Pending update) A custom heatmap plot
"""
function plotheatmap(X, Y, Z; title=nothing, xlab=nothing, ylab=nothing, zmin=nothing, zmax=nothing)
    trace = PlotlyJS.heatmap(
        x=X, y=Y, z=Z, zmin=zmin, zmax=zmax
    )
    layout = PlotlyJS.Layout(
        title=title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        showscale=false,
    )
    p = PlotlyJS.plot(trace, layout)
    return p
end

"""
    plotvolume(X, Y, Z, V; xlab, ylab, zlab, title, isomin, isomax, surface_count=10)

Plot a volume from 3D index matrixes `X`, `Y`, and `Z`, and the 3D value matrix `V`

# Keyword Arguments
- `isomin` and `isomax` : limits of `v` for the color gradient
- `surface_count` : number of separation surfaces
"""
function plotvolume(X, Y, Z, V; xlab::String="x", ylab::String="y", zlab::String="z", title::Union{String, Missing}=missing, isomin=nothing, isomax=nothing, surface_count=10)
    trace = PlotlyJS.volume(
        x=X[:], y=Y[:], z=Z[:], value=V[:],
        isomin=isomin, isomax=isomax,
        opacity=0.20, surface_count=surface_count,
    )
    layout = PlotlyJS.Layout(
        title=title,
        scene=PlotlyJS.attr(
            xaxis_title=xlab,
            xaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
                zerolinecolor="rgba(0,0,0,1.00)",
            ),
            yaxis_title=ylab,
            yaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
                zerolinecolor="rgba(0,0,0,1.00)",
            ),
            zaxis_title=zlab,
            zaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
                zerolinecolor="rgba(0,0,0,1.00)",
            ),
        ),
    )
    p = PlotlyJS.plot(trace, layout)
    return p
end

"""
    plotvolume(xrange, yrange, zrange, f::Function; kwargs...)

Plot a volume from range vectors and a function `f`
"""
function plotvolume(xrange, yrange, zrange, f::Function; kwargs...)
    mesh = Iterators.product(xrange, yrange, zrange)
    X = map(m -> m[1], mesh)
    Y = map(m -> m[2], mesh)
    Z = map(m -> m[3], mesh)
    V = f.(X, Y, Z)
    return plotvolume(X, Y, Z, V; kwargs...)
end

"""
    plotcluster(data::Matrix, label::Vector)

Plot 2D or 3D projected `data` with respective `label`
"""
function plotcluster(data::Matrix, label::Vector)
    p = Plots.plot()
    is3D = size(data)[1] == 3
    if is3D
        for l ∈ sort(unique(label))
            idx = findall(==(l), label)
            Plots.plot!(p, data[1, idx], data[2, idx], data[3, idx], seriestype=:scatter3d, label=string(l), opacity=0.75)
        end
    else
        for l ∈ sort(unique(label))
            idx = findall(==(l), label)
            Plots.plot!(p, data[1, idx], data[2, idx], seriestype=:scatter, label=string(l), opacity=0.75)
        end
    end
    return p
end

"""
    plotseries_percentile(mtx::Matrix; xlab, ylab::String, title, bsort::Bool=false, linealpha=0.10, ylims)

Plot multiple timeseries with the same length from `mtx`, a matrix of (time x observation)

# Keyword Arguments
- `bsort` : sort the series
- `linealpha` : transparency of the original series
"""
function plotseries_percentile(mtx::Matrix; xlab::String="Hours in a year", ylab::String="Power(kW)", title::String="", bsort::Bool=false, linealpha=0.10, ylims=:auto)
    pf_pct50 = Statistics.mean(mtx, dims=2)[:]
    pf_pct25 = zero(pf_pct50)
    pf_pct75 = zero(pf_pct50)
    for idx ∈ 1:length(pf_pct50)
        pf_pct25[idx] = StatsBase.quantile(mtx[idx, :], 0.25)
        pf_pct75[idx] = StatsBase.quantile(mtx[idx, :], 0.75)
    end
    if bsort
        sort!(pf_pct50); sort!(pf_pct25); sort!(pf_pct75)
        mtx = sort(mtx, dims=1)
    end
    p = Plots.plot(title=title, xlabel=xlab, ylabel=ylab, ylims=ylims)
    [Plots.plot!(p, mtx[:, iobj], label=nothing, linecolor=:black, linealpha=linealpha) for iobj ∈ 1:(size(mtx)[2])]
    Plots.plot!(p, pf_pct50, label="pct50", linecolor=Plots.palette(:default)[1])
    Plots.plot!(p, pf_pct75, label="pct75", linecolor=Plots.palette(:default)[2])
    Plots.plot!(p, pf_pct25, label="pct25", linecolor=Plots.palette(:default)[3])
    return p
end

"""
    plotbox(df::DataFrame; xlab, ylab, title, col_variable, col_value)

Plot a box charge from `df`, a dataframe with columns `:variable` and `:value`

# Keyword Arguments
- `col_variable` and `col_value` as `Symbol` : overwrite the default column names
"""
function plotbox(df::DataFrames.DataFrame; xlab::String="Scenario", ylab::String="", title::Union{String, Missing}=missing, col_variable=:variable, col_value=:value)
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    for gd ∈ DataFrames.groupby(df, col_variable)
        push!(pTraces, PlotlyJS.box(y=gd[:, col_value], name=gd[1, col_variable]))
    end
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.0)", # Transparent plot BG
        paper_bgcolor="rgba(255,255,255,1.0)", # White paper BG
        title=title,
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=false, legend=PlotlyJS.attr(orientation="h"),
        barmode="overlay",
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end

"""
    plotbox(vt::Vector; kwargs...)

Format `vt` to a proper structure, and pass it to `plotbox`
"""
plotbox(vt::Vector; kwargs...) = plotbox(DataFrames.DataFrame(:variable => "", :value => vt); kwargs...)

"""
    plotmap(df_point, df_line; minscattersize, maxscattersize, opacitymarker, zoom)

Under development

df_point: id (optional, needed for connection plot), cluster(optional), longitude, latitude, value(optional)

df_line: id_i, id_j, value(optional)

mapstyle "carto-positron", "open-street-map"
"""
function plotmap(df_point, df_line=missing; minscattersize=3., maxscattersize=8., opacitymarker=0.8, zoom=13, mapstyle="carto-positron")
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    df_point = deepcopy(df_point)
    if !ismissing(df_line)
        df_line = deepcopy(df_line)
        line_lat = []; line_lon = [];
        for r ∈ DataFrames.eachrow(df_line)
            idx_i = findfirst(==(r.id_i), df_point.id)
            idx_j = findfirst(==(r.id_j), df_point.id)
            push!(line_lat, df_point[idx_i, :latitude] , df_point[idx_j, :latitude] , missing)
            push!(line_lon, df_point[idx_i, :longitude], df_point[idx_j, :longitude], missing)
        end
        # If there is no value column, add one.
        push!(pTraces, PlotlyJS.scattermapbox(lat=line_lat, lon=line_lon, name="lines", mode="lines"))
    end
    (:id ∉ propertynames(df_point)) && (df_point.id .= missing)
    (:cluster ∉ propertynames(df_point)) && (df_point.cluster .= "")
    (:value ∉ propertynames(df_point)) && (df_point.value .= 0.)
    minvalue = minimum(df_point.value); maxvalue = maximum(df_point.value)
    df_point.scattersize = (df_point.value .- minvalue) / (==(minvalue, maxvalue) ? 1.0 : maxvalue-minvalue)
    df_point.scattersize = (maxscattersize-minscattersize) * df_point.scattersize .+ minscattersize
    for gd ∈ DataFrames.groupby(df_point, :cluster)
        push!(pTraces, PlotlyJS.scattermapbox(
            lat=gd.latitude, lon=gd.longitude, text=gd.id,
            marker=PlotlyJS.attr(size=gd.scattersize),
            opacity=opacitymarker,
            name=first(gd.cluster), mode="markers",
        ))
    end
    showlegend = (length(unique(df_point.cluster)) > 1)
    pLayout = PlotlyJS.Layout(
        hovermode="closest",
        mapbox=PlotlyJS.attr(
            style=mapstyle, zoom=zoom, accesstoken="dummytoken",
            center=PlotlyJS.attr(lat=sum(df_point.latitude)/DataFrames.nrow(df_point), lon=sum(df_point.longitude)/DataFrames.nrow(df_point))
        ),
        margin=PlotlyJS.attr(l=1, r=1, t=1, b=1),
        showlegend=showlegend, legend=PlotlyJS.attr(orientation="h"),
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end

# Data interface functions ====================================================

"""
    save_dftoh5(filename::String, objectname::String, df::DataFrame; col_value=:value)

Save `df` as an object (folder) into a HDF5 file

`df` must be in a long format with only one value column. In order to save memory, indexes are treated and stored as CategoricalArrays, i.e., with levels and levelcode.
"""
function save_dftoh5(filename::String, objectname::String, df::DataFrames.DataFrame; col_value=:value)
    # Establish connection and create group
    fid = HDF5.h5open(filename, "cw")
    objectname ∈ HDF5.keys(fid) && HDF5.delete_object(fid, objectname)
    gid = HDF5.create_group(fid, objectname)
    # Save index columns
    indexcols = setdiff(DataFrames.propertynames(df), [col_value])
    for col ∈ indexcols
        tmpCol = CategoricalArrays.categorical(df[!, col], compress=true)
        HDF5.write_dataset(gid, "idset_" * string(col), CategoricalArrays.levels(tmpCol))
        HDF5.write_dataset(gid, "idlvl_" * string(col), CategoricalArrays.levelcode.(tmpCol))
    end
    # Save the value column
    HDF5.write_dataset(gid, "value_" * string(col_value), df[!, col_value])
    # Close connections
    HDF5.close(gid)
    HDF5.close(fid)
    return nothing
end

"""
    load_h5todf(filename::String, objectname::String)

Load a saved `df` (folder) from a HDF5 file

See also : [`save_dftoh5`](@ref)
"""
function load_h5todf(filename::String, objectname::String)
    fid = HDF5.h5open(filename, "r")
    gid = fid[objectname]
    # reconstruct the dataframe
    df = DataFrames.DataFrame()
    for colkey ∈ [x for x ∈ HDF5.keys(gid) if !contains(x, "idlvl")]
        if contains(colkey, "value")
            colsym = Symbol(replace(colkey, "value_" => ""))
            DataFrames.insertcols!(df, colsym => HDF5.read(gid, colkey))
        elseif contains(colkey, "idset")
            idset = HDF5.read(gid, colkey)
            idlvl = HDF5.read(gid, replace(colkey, "idset_" => "idlvl_"))
            colvalue = idset[idlvl]
            colsym = Symbol(replace(colkey, "idset_" => ""))
            DataFrames.insertcols!(df, colsym => colvalue)
        end
    end
    HDF5.close(fid)
    return df
end

"""
    loadall_h5todf(filename::String)

Load all saved `df` (folder) from a HDF5 file

See also : [`load_h5todf`](@ref)
"""
function loadall_h5todf(filename::String)
    fid = HDF5.h5open(filename, "r")
    listallobject = HDF5.keys(fid)
    HDF5.close(fid)
    data = Dict{String, DataFrames.DataFrame}()
    for objname ∈ listallobject
        df = load_h5todf(filename, objname)
        data[objname] = df
    end
    return data
end

"""
    save_dftodb(dbpath::String, tablename::String, df::DataFrame)

Save `df` as a table in a SQLite database
"""
function save_dftodb(dbpath::String, tablename::String, df::DataFrames.DataFrame)
    db = SQLite.DB(dbpath)
    SQLite.drop!(db, tablename)
    SQLite.load!(df, db, tablename)
end

"""
    load_dbtodf(dbpath::String, tablename::String)

Load a table from a SQLite database as DataFrame
"""
function load_dbtodf(dbpath::String, tablename::String)
    db = SQLite.DB(dbpath)
    df = DataFrames.DataFrame(DBInterface.execute(db, "SELECT * FROM ($tablename)"))
    return df
end

"""
    list_dbtable(dbpath)

List tables in a database
"""
list_dbtable(dbpath) = [x.name for x ∈ SQLite.tables(SQLite.DB(dbpath))]

"""
    appendtxt(filename::String, text::String)

Append `text` a text file
"""
function appendtxt(filename::String, text::String)
    io = open(filename, "a")
    if length(text) == 0
        write(io, "")
    else
        write(io, text * '\n')
    end
    close(io)
end

# Profile modification functions ==============================================

"""
    generatedailypattern(nCycle, angleoffset, ΔT; bfullwave::Bool)

Synthesize a daily profile pattern using a minus cosine function
"""
function generatedailypattern(nCycle, angleoffset, ΔT; bfullwave::Bool=false)
    pattern = -cos.(range(0, nCycle*2*pi, length=Int(24/ΔT)) .+ angleoffset)
    if !bfullwave
        pattern[findall(pattern .< 0)] .= 0
    end
    return pattern
end

"""
    generatepoissonseries(n::Int, λ::Int)

Generate a random Poisson serie with a mean of `λ` whose sum equals to `n`
"""
function generatepoissonseries(n::Int, λ::Int)
    dist = Distributions.Poisson(λ)
    series = rand(dist, n ÷ λ)
    while sum(series) != n
        if sum(series) < n
            push!(series, rand(dist))
        else
            diff = n - sum(series[1:(end-1)])
            if diff > 0
                series[end] = diff
            else
                pop!(series)
            end
        end
    end
    return series
end

"""
    synthesizeprofile(avgprofile::Vector{Float64}, λ::Int; base_rel::Float64, base_fix::Float64)

Synthesize a random profile from `avgprofile` by grouping neighboring values

Neighbors are randomly chosen using a poissonseries whose mean is `λ`

# Keyword Arguments
- `base_rel` : base value as share of the `avgprofile`, default 0.1
- `base_fix` : base value as a fixed quantity, default 0.0

See also: [`generatepoissonseries`](@ref)
"""
function synthesizeprofile(avgprofile::Vector{Float64}, λ::Int; base_rel::Float64=0.1, base_fix::Float64=0.0)
    nts = length(avgprofile)
    # Shift avgprofile
    moveindex = rand(1:nts)
    avgprofile = circshift(avgprofile, moveindex)
    # Generate random sequence
    randomseries = generatepoissonseries(nts, λ)
    deleteat!(randomseries, findall(==(0), randomseries))
    # Synthesize a profile
    synprofile = min.(base_rel * avgprofile .+ base_fix, avgprofile)
    comp_rest = avgprofile .- synprofile
    ind_first = 1
    for interval ∈ randomseries
        ind_last = ind_first + interval - 1
        synprofile[ind_first + (interval ÷ 2)] += sum(comp_rest[ind_first:ind_last])
        ind_first += interval
    end
    # Shift back synprofile
    circshift!(synprofile, -moveindex)
    return synprofile
end

"""
    normalize(vt::Vector; mode::Symbol=:unitpeak, ϵ=1e-6)

Normalize the vector `vt`

# Keyword Arguments
- `mode`: either :unitpeak, :unitsum, or :unitvariance
- `ϵ`: criteria to reject the normalization, e.g., when `vt` is a zero vector.
"""
function normalize(vt::Vector; mode::Symbol=:unitpeak, ϵ=1e-6)
    if mode == :unitpeak
        min_pf = minimum(skipmissing(vt)); max_pf = maximum(skipmissing(vt))
        return (abs(max_pf - min_pf) > ϵ) ? ((vt .- min_pf) / (max_pf - min_pf)) : vt
    elseif mode == :unitsum
        sum_pf = sum(skipmissing(vt))
        return (abs(sum_pf) > ϵ) ? (vt / sum_pf) : vt
    elseif mode == :unitvariance
        std_pf = Statistics.std(skipmissing(vt))
        return (std_pf > ϵ) ? (vt / std_pf) : vt
    else
        @warn "Unrecognized mode"
        return repeat([missing], length(vt))
    end
end

# Miscellaneous functions =====================================================

"""
    clippy(obj)

Copy an object `obj` into system's clipboard
"""
clippy(obj) = Main.clipboard(sprint(show, "text/tab-separated-values", obj))

"""
    createdummydata(nYear, nTime, nRegion, nVariable)

Create a dummy DataFrame with index columns `:year`, `:time`, `:region`, and `:variable`, and the value column `:value`
"""
function createdummydata(nYear, nTime, nRegion, nVariable)
    sYear = 2025 .+ collect(1:nYear)
    sTime = collect(1:nTime)
    sRegion = [Random.randstring(8) for _ ∈ 1:nRegion]
    sVariable = [Random.randstring(8) for _ ∈ 1:nVariable]
    df = DataFrames.crossjoin(
        DataFrames.DataFrame(:year => sYear),
        DataFrames.DataFrame(:time => sTime),
        DataFrames.DataFrame(:region => sRegion),
        DataFrames.DataFrame(:variable => sVariable),
    )
    DataFrames.insertcols!(df, :value => round.(1000*rand(DataFrames.nrow(df)), digits=2))
    return df
end

"""
    mergeposneg(dfpos, dfneg; col_value = :value)

Merge two dataframe representing positive values and negative values into a single bidirectional data.
"""
function mergeposneg(dfpos::DataFrames.DataFrame, dfneg::DataFrames.DataFrame; col_value=:value)
    dfpos_fmt = DataFrames.rename(dfpos, col_value => :pos)
    dfneg_fmt = DataFrames.rename(dfneg, col_value => :neg)
    df = DataFrames.outerjoin(dfpos_fmt, dfneg_fmt, on=setdiff(DataFrames.propertynames(dfpos), [col_value]))
    df[ismissing.(df.pos), :pos] .= 0.0
    df[ismissing.(df.neg), :neg] .= 0.0
    DataFrames.transform!(df, [:pos, :neg] => ((x, y) -> x .- y) => :value)
    DataFrames.select!(df, DataFrames.Not([:pos, :neg]))
    return df
end

"""
    averageprofile(pf::Vector; Δt=1.0, bseason=true, bufferlength=3, idx_winter=missing, idx_summer=missing)

Process average daily profile from an entire year's profile `pf`

# Keyword Arguments
- `Δt` : lenght of a timestep in [hour]
- `bseason` : separate individual profiles for summer, winter, and transition period
- `bufferlength` : empty space between seasonal profiles
- `idx_winter` and `idx_summer` : boolean vectors indicating winter and summer. If missing, the default German winter and summer days are used. The default assumes that one year has 365 days.
"""
function averageprofile(pf::Vector; Δt=1.0, bseason=true, bufferlength=3, idx_winter=missing, idx_summer=missing)
    winterday_in_DE = Bool.(vcat(ones(59), zeros(275), ones(31)))
    summerday_in_DE = Bool.(vcat(zeros(151), ones(92), zeros(122)))
    ts_in_day = Int(24/Δt)
    tod = repeat(1:ts_in_day, outer=length(pf)÷ts_in_day)
    if bseason
        # Process season indexes
        ismissing(idx_winter) && (idx_winter = repeat(winterday_in_DE, inner=ts_in_day))
        ismissing(idx_summer) && (idx_summer = repeat(summerday_in_DE, inner=ts_in_day))
        df = DataFrames.DataFrame(:tod => tod, :value => pf, :season => ".")
        df[idx_winter, :season] .= "winter"
        df[idx_summer, :season] .= "summer"
        df[.!(idx_summer .| idx_winter), :season] .= "transition"
        df = DataFrames.combine(DataFrames.groupby(df, [:tod, :season]), :value => Statistics.mean => :value)
        DataFrames.sort!(df, [:season, :tod])
        pf_avg = vcat(reshape(df.value, :, 3), zeros(bufferlength, 3))[:][1:(end-bufferlength)]
    else
        df = DataFrames.DataFrame(:tod => tod, :value => pf)
        df = DataFrames.combine(DataFrames.groupby(df, :tod), :value => Statistics.mean => :value)
        DataFrames.sort!(df, :tod)
        pf_avg = df.value
    end
    return pf_avg
end

"""
    equipmentloading(df; col_region=:region, col_value=:value)

Calculate loading factors from time-resolved operation profiles

Loading factor is defined as the ratio of average and peak values. `df` contains time, region, and value columns. The calculation uses the absolute values.
"""
function equipmentloading(df; col_region=:region, col_value=:value)
    df = deepcopy(df)
    df[:, col_value] = abs.(df[:, col_value])
    loading = DataFrames.combine(DataFrames.groupby(df, col_region), col_value => maximum => :peak, col_value => Statistics.mean => :avg)
    DataFrames.select!(loading, :region, [:avg, :peak] => ((x,y) -> x ./ y) => :value)
    loading[isnan.(loading.value), :value] .= 0.0
    return loading
end

"""
    df_filter_collapse_plot(df, filter_dc, keep_cols; f_collapse=sum, col_time=:time, col_value=:value, kwargs...)

Basic processing of `df`, then plot either the timeseries or bar plot

Filter `df` using criteria in `filter_dc`, a dictionary mapping a column as a symbol to the corresponding filter function. Collapse the filtered dataframe to only the `keep_cols` and `col_value`, the value column. Finally, the processed `df` is plotted either as a timeseries plot, or a bar plot.

# Keyword Arguments
- `f_collapse` : a function to collapse the value column, default is `sum`
- `col_time` and `col_value` : overwrite the default column names

See also : [`plottimeseries`](@ref), [`plotbar`](@ref)
"""
function df_filter_collapse_plot(df, filter_dc, keep_cols; f_collapse=sum, col_time=:time, col_value=:value, kwargs...)
    df = DataFrames.deepcopy(df)
    [DataFrames.filter!(col => filter_dc[col], df) for col ∈ keys(filter_dc)]
    df = DataFrames.combine(DataFrames.groupby(df, keep_cols), col_value => f_collapse => col_value)
    # Call `plottimeseries` or `plotbar`
    if col_time ∈ keep_cols
        col_variable = setdiff(keep_cols, [col_time])[1]
        p = plottimeseries(df, col_variable=col_variable, col_time=col_time; kwargs...)
    else
        col_axis = keep_cols[1]; col_variable = keep_cols[2]
        p = plotbar(df, col_axis=col_axis, col_variable=col_variable; kwargs...)
    end
    return p
end

"""
    getcolor(key; colorcode=colorcode)

Fetch an RGBA color from `colorcode`, a global dictionary variable

If the key does not exist, then random a color and update the `colorcode`
"""
function getcolor(key; colorcode=colorcode)
    if key ∉ keys(colorcode)
        tmp = lpad.(rand(0:255, 3), 3, "0")
        colorcode[key] = "rgba($(tmp[1]), $(tmp[2]), $(tmp[3]), 0.5)"
    end
    return colorcode[key]
end

# Pending retirement ==========================================================

end
