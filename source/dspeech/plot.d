/// speech plot module on GGplotD
module dspeech.plot;

import std.array : array;
import std.algorithm : map;

import ggplotd.aes : aes;
import ggplotd.axes : xaxisLabel, yaxisLabel;
import ggplotd.ggplotd : GGPlotD, putIn;
import ggplotd.geom : geomPoint, geomRectangle;
import ggplotd.colour : colourGradient;
import ggplotd.colourspace : XYZ;


auto geomPointRect(AES)(AES aesRange)
{
    import ggplotd.aes : aes, Pixel, DefaultValues, merge;
    import ggplotd.range : mergeRange;

    return DefaultValues.mergeRange(aesRange)
        .map!((a) => a.merge(aes!("sizeStore", "width", "height", "fill")
        (a.size, a.width, a.height, a.alpha))).geomRectangle;
}

auto plotMatrix(T)(T array2d)
{
    import std.algorithm : cartesianProduct;
    import std.range : iota;
    auto xstep = 1;
    auto ystep = 1;
    auto xlen = array2d[0].length;
    auto ylen = array2d.length;
    auto xys = cartesianProduct(xlen.iota, ylen.iota);
    auto gg = xys.map!(xy => aes!("x", "y", "colour", "size", "width",
            "height")(xy[0], xy[1], array2d[$-1-xy[1]][xy[0]], 1.0, xstep, ystep))
        .array.geomPointRect.putIn(GGPlotD());
    gg = colourGradient!XYZ("mediumblue-limegreen-orangered").putIn(gg);
    gg = "time".xaxisLabel.putIn(gg);
    gg = "freq".yaxisLabel.putIn(gg);
    return gg;
}

auto plotVector(T)(T array)
{
    import std.range : enumerate;
    import ggplotd.geom : geomLine;

    auto gg = GGPlotD();
    gg = enumerate(array).map!(a => aes!("x", "y", "colour", "size")(a[0], a[1], 0, 0.1))
        .array.geomLine.putIn(gg);
    gg = "time".xaxisLabel.putIn(gg);
    gg = "gain".yaxisLabel.putIn(gg);
    return gg;
}
