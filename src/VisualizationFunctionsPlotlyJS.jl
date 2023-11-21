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
