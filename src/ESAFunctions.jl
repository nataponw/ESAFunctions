module ESAFunctions

# Import dependencies =========================================================
import DataFrames, CategoricalArrays
import SQLite, DBInterface, HDF5
import PlotlyJS
import JuMP
import Random, Distributions

# Declare export ==============================================================
# Functions related to the CommOpt model (to be relocated)
export gp, showjo, extractDBTable
# Visualization functions
export plottimeseries, plothistogram, plotcontour, plotheatmap, plotvolume, plotslicevolumn
# Data interface functions
export save_dftoh5, load_h5todf, loadall_h5todf, save_dftodb, load_dbtodf, list_dbtable, writetxt
# Profile modification functions
export generatedailypattern, generatepoissonseries, synthesizeprofile
# Miscellaneous functions
export clippy, createdummydata
# Pending retirement

# Functions related to the CommOpt model ======================================
"""
Select a scalar parameter value from a DataFrame
"""
function gp(df::DataFrames.DataFrame, selcol::Symbol, filter::Vector)
    index = collect(1:size(df)[1])
    for i in 1:length(filter)
        if !ismissing(filter[i])
            intersect!(index, findall(df[:, i] .== filter[i]))
        end
    end
    if isempty(index)
        return missing
    else
        return df[index, selcol][1]
    end
end
# Dispatches of `gp`
function gp(df::DataFrames.DataFrame, selcol::Symbol, filter::Symbol)
    return gp(df, selcol, [filter])
end
function gp(df::DataFrames.DataFrame, selcol::Symbol, filter::Int64)
    return gp(df, selcol, [filter])
end

"""
show a JuMP DenseAxisArray as a DataFrame with appropriate axes
"""
function showjo(obj::JuMP.Containers.DenseAxisArray; namedim1::Symbol=:dim1)
    # Check dimension
    dim = length(obj.axes)
    values = JuMP.value.(obj).data
    if dim == 1
        df = DataFrames.DataFrame(namedim1 => obj.axes[1], :value => values)
    elseif dim == 2
        if isa(obj.axes[2][1], Symbol)
            df = DataFrames.DataFrame(values, obj.axes[2])
        else
            df = DataFrames.DataFrame(values, :auto)
        end
        DataFrames.insertcols!(df, 1, namedim1 => obj.axes[1])
    else
        println("Warning: the display of higher dimensions is not yet supported.")
    end
    return df
end

"""
Extract a table from a database `db` into a DataFrame
"""
function extractDBTable(db, tableName::String, symbolColumn::Vector{Int}=Int[])
    df = DataFrames.DataFrame(DBInterface.execute(db, "SELECT * FROM ($tableName)"))
    if !isempty(symbolColumn)
        for i ∈ symbolColumn
            df[!, i] = Symbol.(df[!, i])
        end
    end
    return df
end

# Visualization functions =====================================================
"""
Plot timeseries from a dataframe `df` containing columns [:time, :variable, :value]
"""
function plottimeseries(df::DataFrames.DataFrame;
    xlab::String="Time", ylab::String="Power (kW)", title::Union{String, Missing}=missing,
    col_time=:time, col_variable=:variable, col_value=:value,
    bstack::Bool=false
    )
    # Plot settings
    stackgroup = (bstack ? "one" : missing)
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    for gd ∈ DataFrames.groupby(df, col_variable)
        push!(pTraces, PlotlyJS.scatter(x=gd[:, col_time], y=gd[:, col_value], name=gd[1, col_variable], mode="lines", stackgroup=stackgroup))
    end
    # Do not show legends when there is only one trace
    showlegend = length(pTraces) > 1
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.10)", # Translucent plot BG
        paper_bgcolor="rgba(0,0,0,0.0)", # No paper BG
        title=title,
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=showlegend, legend=PlotlyJS.attr(orientation="h"),
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end
# Dispatches of `plottimeseries`
function plottimeseries(dt::Dict; kwargs...)
    df = DataFrames.DataFrame()
    [df = vcat(df, DataFrames.DataFrame(:time => 1:length(dt[k]), :variable => k, :value => dt[k])) for k ∈ keys(dt)]
    return plottimeseries(df; kwargs...)
end
plottimeseries(vt::Vector; kwargs...) = plottimeseries(df; kwargs...)

"""
Plot histogram from `dt`, a dictionary mapping between symbols and respective arrays
"""
function plothistogram(dt::Dict; xlab::String="Value", ylab::String="Count")
    # Process data
    df = DataFrames.DataFrame()
    for iLine ∈ keys(dt)
        if dt[iLine] isa Matrix{<: Number}
            DataFrames.insertcols!(df, iLine => sum(dt[iLine], dims=1)[:])
        elseif dt[iLine] isa Vector{<: Number}
            DataFrames.insertcols!(df, iLine => dt[iLine])
        else break
        end
    end
    # Plot
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    lvlOpacity = (DataFrames.ncol(df) == 1 ? 1.00 : 0.60)
    for iLine ∈ names(df)
        push!(pTraces, PlotlyJS.histogram(x=df[:, iLine], name=iLine, opacity=lvlOpacity))
    end
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.10)", # Translucent plot BG
        paper_bgcolor="rgba(0,0,0,0)", # No paper BG
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=true, legend=PlotlyJS.attr(orientation="h"),
        barmode="overlay",
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end
plothistogram(dt::Vector; xlab::String="Value", ylab::String="Count") = plothistogram(Dict(:data => dt); xlab=xlab, ylab=ylab)

"""
A custom contour plot
"""
function plotcontour(X, Y, Z; title=nothing, xlab=nothing, ylab=nothing, zmin=nothing, zmax=nothing, bsave=false, plotname="test.png")
    trace = PlotlyJS.contour(
        x=X, y=Y, z=Z, zmin=zmin, zmax=zmax,
        contours=PlotlyJS.attr(
            showlabels=true, 
            labelfont = PlotlyJS.attr(color="darkgray")
        )
    )
    layout = PlotlyJS.Layout(
        title=title,
        xaxis_title=xlab,
        yaxis_title=ylab,
        showscale=false,
    )
    p = PlotlyJS.plot(trace, layout)
    bsave && PlotlyJS.savefig(p, plotname, width=650, height=600, scale=2.0, format="png")
    return p
end

"""
A custom heatmap plot
"""
function plotheatmap(X, Y, Z; title=nothing, xlab=nothing, ylab=nothing, zmin=nothing, zmax=nothing, bsave=false, plotname="test.png")
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
    bsave && PlotlyJS.savefig(p, plotname, width=650, height=600, scale=2.0, format="png")
    return p
end

"""
A custom volume plot
"""
function plotvolume(X, Y, Z, V; title=nothing, xlab=nothing, ylab=nothing, zlab=nothing, isomin=nothing, isomax=nothing, surface_count=10)
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
            ),
            yaxis_title=ylab,
            yaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
            ),
            zaxis_title=zlab,
            zaxis=PlotlyJS.attr(
                showbackground=false,
                gridcolor="rgba(200,200,200,0.9)",
            ),
        ),
    )
    p = PlotlyJS.plot(trace, layout)
    return p
end

"""
Slice a volume at a specific value of one of the dimension
"""
function plotslicevolumn(X, Y, sliceDim::Int, slicePoint)
    indexlist = findall(X[:, sliceDim] .== slicePoint)
    if sliceDim == 1
        return plotcontour(X[indexlist, 2], X[indexlist, 3], Y[indexlist], title="Volume sliced at dim $(sliceDim) : $(slicePoint)", xlab="Dim 2", ylab="Dim 3")
    elseif sliceDim == 2
        return plotcontour(X[indexlist, 1], X[indexlist, 3], Y[indexlist], title="Volume sliced at dim $(sliceDim) : $(slicePoint)", xlab="Dim 1", ylab="Dim 3")
    else
        return plotcontour(X[indexlist, 1], X[indexlist, 2], Y[indexlist], title="Volume sliced at dim $(sliceDim) : $(slicePoint)", xlab="Dim 1", ylab="Dim 2")
    end
end

# Data interface functions ====================================================

"""
Generic save function of a dataframe into a H5 file
"""
function save_dftoh5(filename::String, objectname::String, df::DataFrames.DataFrame; indexcols=[])
    # Establish connection and create group
    fid = HDF5.h5open(filename, "cw")
    objectname ∈ HDF5.keys(fid) && HDF5.delete_object(fid, objectname)
    gid = HDF5.create_group(fid, objectname)
    # If `indexcols` is empty, then columns whose data are not Float are considered index columns
    if isempty(indexcols)
        for col ∈ DataFrames.propertynames(df)
            if !(eltype(df[!, col]) .<: Union{Missing, Real})
                push!(indexcols, col)
            end
        end
    end
    # Process each columns
    for col ∈ DataFrames.propertynames(df)
        if col ∈ indexcols
            # First implementation: sequencially process and write each index columns
            tmpCol = CategoricalArrays.categorical(df[!, col], compress=true)
            HDF5.write_dataset(gid, "idset_" * string(col), CategoricalArrays.levels(tmpCol))
            HDF5.write_dataset(gid, "idlvl_" * string(col), CategoricalArrays.levelcode.(tmpCol))
            # Second implementation: first process all columns and then write them all at once <<TBD>>
        else
            # HDF5.write_dataset(gid, "value_" * string(col), Float32.(df[!, col]))
            HDF5.write_dataset(gid, "value_" * string(col), df[!, col])
        end
    end
    # Close connections
    HDF5.close(gid)
    HDF5.close(fid)
    return nothing
end

"""
Generic load function of a dataframe from a H5 file
"""
function load_h5todf(filename::String, objectname::String; caindices=true)
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
            if caindices
                DataFrames.insertcols!(df, colsym => CategoricalArrays.categorical(colvalue))
            else
                DataFrames.insertcols!(df, colsym => colvalue)
            end
        end
    end
    HDF5.close(fid)
    return df
end

"""
Load all dataframes from a hdf5 file
"""
function loadall_h5todf(filename::String; caindices=true)
    fid = HDF5.h5open(filename, "r")
    listallobject = HDF5.keys(fid)
    HDF5.close(fid)
    data = Dict{String, DataFrames.DataFrame}()
    for objname ∈ listallobject
        df = load_h5todf(filename, objname, caindices=caindices)
        data[objname] = df
    end
    return data
end

"""
Save a table to a SQLite database
"""
function save_dftodb(dbpath::String, tableName::String, df::DataFrames.DataFrame)
    db = SQLite.DB(dbpath)
    SQLite.drop!(db, tableName)
    SQLite.load!(df, db, tableName)
end

"""
Load a table from a SQLite database
"""
function load_dbtodf(dbpath::String, tableName::String)
    db = SQLite.DB(dbpath)
    df = DataFrames.DataFrame(DBInterface.execute(db, "SELECT * FROM ($tableName)"))
    return df
end

"""
List of tables in a SQLite database
"""
list_dbtable(dbpath) = [x.name for x ∈ SQLite.tables(SQLite.DB(dbpath))]

"""
Simple write to a text file. The file is reset when `text` is empty.
"""
function writetxt(filename::String, text::String)
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
Generate a daily sinusoidal pattern
"""
function generatedailypattern(cycle, offset, dT; bfullwave::Bool=false)
    pattern = -cos.(range(0, cycle*2*pi, length=Int(24/dT)) .+ offset)
    if !bfullwave
        pattern[findall(pattern .< 0)] .= 0
    end
    return pattern
end

"""
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
Given an average profile, synthesize a disaggregate profile by grouping neighboring values
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

# Miscellaneous functions =====================================================

clippy(df) = Main.clipboard(sprint(show, "text/tab-separated-values", df))

"""
Create a dummy dataframe with [:year, :ts, :region, :variable, :value]
"""
function createdummydata(nY, nTS, nRegion, nVariable; ascategoricalarray::Bool=false)
    sY = 2025 .+ collect(1:nY)
    sTS = collect(1:nTS)
    sRegion = [Random.randstring(8) for i ∈ 1:nRegion]
    sVariable = [Random.randstring(8) for i ∈ 1:nVariable]
    if ascategoricalarray
        df = DataFrames.crossjoin(
            DataFrames.DataFrame(:year => sY),
            DataFrames.DataFrame(:time => sTS),
            DataFrames.DataFrame(:region => CategoricalArrays.categorical(sRegion)),
            DataFrames.DataFrame(:variable => CategoricalArrays.categorical(sVariable)),
        )
    else
        df = DataFrames.crossjoin(
            DataFrames.DataFrame(:year => sY),
            DataFrames.DataFrame(:time => sTS),
            DataFrames.DataFrame(:region => sRegion),
            DataFrames.DataFrame(:variable => sVariable),
        )
    end
    DataFrames.insertcols!(df, :value => round.(1000*rand(DataFrames.nrow(df)), digits=6))
    return df
end

# Pending retirement ==========================================================

end
