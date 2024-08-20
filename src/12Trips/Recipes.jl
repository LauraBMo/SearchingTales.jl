
_coord(nodes) = [first.(nodes), last.(nodes)]
close_loop!(xs) = push!(xs, first(xs))
_unfold_complex(v) = [real.(v), imag.(v)]


_colors(color) = Colors.distinguishable_colors(12, Colors.parse.(Colorant, color))
## Turn on rainbow-mode in emacs vterm console.
_colorsshow(color) = "#" .* Colors.hex.(_colors(color))

@userplot PlotRealPoints

"""
$(SIGNATURES)

```julia
   using Plots
   using PlanarCurves
```
"""
function plotrealpoints end

@recipe function f(pp::PlotRealPoints)
    _points, partition = pp.args

    # No legend
    legend --> false
    colors = plotattributes[:colors]

    xs, ys = _coord(_points)
    for (col, T) in zip(colors, partition)
        @series begin
            label := nothing
            color := col
            seriestype := :scatter
            xs[T], ys[T]
        end
    end
end

@userplot PlotLoop

"""
$(SIGNATURES)

```julia
   using Plots
   using PlanarCurves
```
"""
function plotloop end

@recipe function f(rP::PlotLoop)
    _points, col = rP.args

    xs, ys = _coord(_points)
    close_loop!.([xs, ys])
    @series begin
        label := nothing
        color := col
        seriestype := :path
        xs, ys
    end
end

@userplot PlotPartition

"""
$(SIGNATURES)

```julia
   using Plots
   using PlanarCurves
```
"""
function plotpartition end

@recipe function f(pp::PlotPartition)
    _points, partition = pp.args

    # No legend
    legend --> false
    colors = plotattributes[:colors]

    for (col, T) in zip(colors, partition)
        @series begin
            label := nothing
            PlotLoop((_points[T], col))
        end
    end
end

@userplot PlotNodes

"""
$(SIGNATURES)

```julia
   using Plots
   using PlanarCurves
```
"""
function plotnodes end

@recipe function f(nodes::PlotNodes;
                   color=:grey,
                   colors=_colors(color),
                   op=0.4)
    _nodes, partition = nodes.args
    affine_nodes = _dehomo.(_nodes)

    ## We plot real1 x real2 and imag1 x imag2
    real_pts, imag_pts = _unfold_complex(affine_nodes)

    layout := @layout [o o]
    legend := false
    colors := colors

    opacity := 1
    @series begin
        subplot := 1
        title := "Real"
        PlotRealPoints((real_pts, partition))
    end
    @series begin
        subplot := 2
        title := "Imag"
        PlotRealPoints((imag_pts, partition))
    end

    opacity := op
    @series begin
        subplot := 1
        PlotPartition((real_pts, partition))
    end
    @series begin
        subplot := 2
        PlotPartition((imag_pts, partition))
    end
end

@userplot PlotFit

"""
$(SIGNATURES)

```julia
   using Plots
   using PlanarCurves
```
"""
function plotfit end

@recipe function f(F::PlotFit;
                   color=:grey,
                   colors=_colors(color),
                   op=0.4)
    _fit, = F.args
    curve = _fit.curve
    multipp = _fit.multiplepoints
    nodes = eval_nodes(curve, multipp)

    # color := color
    colors := colors
    PlotNodes((nodes, _fit.partition))
end
