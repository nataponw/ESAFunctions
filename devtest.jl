include("src/ESAFunctions.jl")
using .ESAFunctions
###############################################################################

using DataFrames, PlotlyJS

df = DataFrame(
    :variable => repeat(["AAA", "BBB", "CCC"], 6),
    :x => randn(18),
    :y => randn(18),
)

function plotscatter(df::DataFrames.DataFrame;
    xlab::String="X value", ylab::String="Y value", title::Union{String, Missing}=missing,
    col_variable=:variable, col_x=:x, col_y=:y,
    selectcolor=missing, legendorientation="h",
    )
    # Handle when col_variable is missing.
    col_variable ∉ propertynames(df) && (df[!, col_variable] .= "")
    # Color palette
    ismissing(selectcolor) && (selectcolor = (x -> missing))
    pTraces = PlotlyJS.PlotlyBase.GenericTrace[]
    for gd ∈ DataFrames.groupby(df, col_variable)
        push!(pTraces, PlotlyJS.scatter(x=gd[:, col_x], y=gd[:, col_y], name=gd[1, col_variable], mode="markers", marker=PlotlyJS.PlotlyBase.attr(color=selectcolor(gd[1, col_variable]))))
    end
    showlegend = length(pTraces) > 1
    pLayout = PlotlyJS.Layout(
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
    (minimum(df[:, col_x]) < 0.0) && PlotlyJS.add_vline!(p, 0.0, line_width=1.0, line_dash="dash", line_color="rgba(0,0,0,0.10)")
    (minimum(df[:, col_y]) < 0.0) && PlotlyJS.add_hline!(p, 0.0, line_width=1.0, line_dash="dash", line_color="rgba(0,0,0,0.10)")
    return p
end

gc(k) = getcolor(k, colorcode=colorcode)

global colorcode = Dict(
    "AAA" => "255,0,0,1.0",
    "BBB" => "0,255,0,1.0",
    "CCC" => "0,0,255,1.0",
)

plotscatter(df)

