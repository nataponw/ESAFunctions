module ESAFunctions

import DataFrames, CategoricalArrays
import SQLite, DBInterface, HDF5
import PlotlyJS
import JuMP
import Random

# parameter processing functions
export extractDBTable, gp, showjo
# visualization functions
export plotTimeseries, plotHistogram, plotContour, plotHeatmap, plotVolume, sliceVolume
# Custom I/O
export savedat, loaddat, loadalldat, save_dbtable, load_dbtable, list_dbtable, simpletxtwrite
# Others
export sigmoid, clippy
# potentially obsolete functions
export createDummyDF, generateDailyPattern

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
### Alternative dispatches of `gp`
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
Plot timeseries from `dt`, a dictionary mapping between symbols and respective arrays
"""
function plotTimeseries(dt::Dict; xlab::String="Timesteps", ylab::String="Power (kW)")
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
    for iLine ∈ names(df)
        push!(pTraces, PlotlyJS.scatter(y=df[:, iLine], name=iLine, mode="lines"))
    end
    pLayout = PlotlyJS.Layout(
        xaxis_rangeslider_visible=false,
        plot_bgcolor="rgba(255,255,255,0.10)", # Translucent plot BG
        paper_bgcolor="rgba(0,0,0,0.0)", # No paper BG
        xaxis_title=xlab,
        xaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        yaxis_title=ylab,
        yaxis=PlotlyJS.attr(linecolor="rgba(0,0,0,0.10)"),
        showlegend=true, legend=PlotlyJS.attr(orientation="h"),
    )
    p = PlotlyJS.plot(pTraces, pLayout)
    return p
end
plotTimeseries(dt::Vector; xlab::String="Timesteps", ylab::String="Power (kW)") = plotTimeseries(Dict(:data => dt); xlab=xlab, ylab=ylab)

"""
Plot histogram from `dt`, a dictionary mapping between symbols and respective arrays
"""
function plotHistogram(dt::Dict; xlab::String="Value", ylab::String="Count")
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
plotHistogram(dt::Vector; xlab::String="Value", ylab::String="Count") = plotHistogram(Dict(:data => dt); xlab=xlab, ylab=ylab)

"""
A custom contour plot
"""
function plotContour(X, Y, Z; title=nothing, xlab=nothing, ylab=nothing, zmin=nothing, zmax=nothing, bsave=false, plotname="test.png")
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
function plotHeatmap(X, Y, Z; title=nothing, xlab=nothing, ylab=nothing, zmin=nothing, zmax=nothing, bsave=false, plotname="test.png")
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
function plotVolume(X, Y, Z, V; title=nothing, xlab=nothing, ylab=nothing, zlab=nothing, isomin=nothing, isomax=nothing, surface_count=10)
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
function sliceVolume(X, Y, sliceDim::Int, slicePoint)
    indexlist = findall(X[:, sliceDim] .== slicePoint)
    if sliceDim == 1
        return plotContour(X[indexlist, 2], X[indexlist, 3], Y[indexlist], title="Volume sliced at dim $(sliceDim) : $(slicePoint)", xlab="Dim 2", ylab="Dim 3")
    elseif sliceDim == 2
        return plotContour(X[indexlist, 1], X[indexlist, 3], Y[indexlist], title="Volume sliced at dim $(sliceDim) : $(slicePoint)", xlab="Dim 1", ylab="Dim 3")
    else
        return plotContour(X[indexlist, 1], X[indexlist, 2], Y[indexlist], title="Volume sliced at dim $(sliceDim) : $(slicePoint)", xlab="Dim 1", ylab="Dim 2")
    end
end

"""
Generic save function of a dataframe into a H5 file
"""
function savedat(filename::String, objectname::String, df::DataFrames.DataFrame; indexcols=[])
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
function loaddat(filename::String, objectname::String; caindices=true)
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
function loadalldat(filename::String; caindices=true)
    fid = HDF5.h5open(filename, "r")
    listallobject = HDF5.keys(fid)
    HDF5.close(fid)
    data = Dict{String, DataFrames.DataFrame}()
    for objname ∈ listallobject
        df = loaddat(filename, objname, caindices=caindices)
        data[objname] = df
    end
    return data
end

"""
Load a table from a SQLite database
"""
function load_dbtable(dbpath::String, tableName::String)
    db = SQLite.DB(dbpath)
    df = DataFrames.DataFrame(DBInterface.execute(db, "SELECT * FROM ($tableName)"))
    return df
end

"""
Save a table to a SQLite database
"""
function save_dbtable(dbpath::String, tableName::String, df::DataFrames.DataFrame)
    db = SQLite.DB(dbpath)
    SQLite.drop!(db, tableName)
    SQLite.load!(df, db, tableName)
end

"""
List of tables in a SQLite database
"""
list_dbtable(dbpath) = [x.name for x ∈ SQLite.tables(SQLite.DB(dbpath))]

"""
Simple write to a text file. The file is reset when `text` is empty.
"""
function simpletxtwrite(filename::String, text::String)
    if length(text) == 0
        io = open(filename, "w")
        write(io, "")
    else
        io = open(filename, "a")
        write(io, text * '\n')
    end
    close(io)
end

# Others ----
sigmoid(x) = @. 1 / (1 + exp(-x))
clippy(df) = Main.clipboard(sprint(show, "text/tab-separated-values", df))

# Potentially obsolete ----

"""
Create a dummy dataframe with [:year, :ts, :region, :variable, :value]
"""
function createDummyDF(nY, nTS, nRegion, nVariable; ascategoricalarray::Bool=false)
    sY = 2025 .+ collect(1:nY)
    sTS = collect(1:nTS)
    sRegion = [Random.randstring(8) for i ∈ 1:nRegion]
    sVariable = [Random.randstring(8) for i ∈ 1:nVariable]
    if ascategoricalarray
        df = DataFrames.crossjoin(
            DataFrames.DataFrame(:year => sY),
            DataFrames.DataFrame(:ts => sTS),
            DataFrames.DataFrame(:region => CategoricalArrays.categorical(sRegion)),
            DataFrames.DataFrame(:variable => CategoricalArrays.categorical(sVariable)),
        )
    else
        df = DataFrames.crossjoin(
            DataFrames.DataFrame(:year => sY),
            DataFrames.DataFrame(:ts => sTS),
            DataFrames.DataFrame(:region => sRegion),
            DataFrames.DataFrame(:variable => sVariable),
        )
    end
    DataFrames.insertcols!(df, :value => round.(1000*rand(DataFrames.nrow(df)), digits=6))
    return df
end

"""
Generate a daily sinusoidal pattern
"""
function generateDailyPattern(cycle, offset, dT; bFullwave::Bool=false)
    pattern = -cos.(range(0, cycle*2*pi, length=Int(24/dT)) .+ offset)
    if !bFullwave
        pattern[findall(pattern .< 0)] .= 0
    end
    return pattern
end

end
