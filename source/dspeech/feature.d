/**
   Speech feature extraction module

   See also: "Chapter 5 Speech Input/Output," HTK book 3.5 alpha1, pp. 89 - 113.
*/
module dspeech.feature;

import numir.signal : hann;

/**
Computes a magnitude spectrogram from a time-domain 1-D signal

Params:
    windowFun = window function
    xs = 1-D slice signal
    nperseg = (default 256) short-time frame width for each FFT segment
    noverlap = (default nperseg / 2) short-time frame overlapped length for each FFT segment
*/
auto spectrogram(alias windowFun=hann, S)(S xs, size_t nperseg = 256, size_t noverlap = 128)
{
    import numir.signal : stft;
    import mir.math.common : sqrt;
    import mir.ndslice.topology : map;
    import mir.ndslice : transposed;

    auto stfs = stft!windowFun(xs, nperseg, noverlap);
    auto upper = stfs.transposed[nperseg / 2 .. $];
    auto mag = upper.map!(a => sqrt(a.re * a.re + a.im * a.im));
    return mag;
}

/// example spectrogram of freq-modulated sin(a t + b cos(c t))
/// its plot should be simple sinwave in freq domain
unittest
{
    import mir.ndslice : as, map, iota;
    import mir.math : sin, cos;
    import std.math : PI;
    import dspeech.plot : plotMatrix;
    auto time = iota(50000).as!double / 10e3;
    auto mod = 1e3 * map!cos(2.0 * PI * 0.5 * time);
    auto xs = map!sin(2.0 * PI * 3e3 * time  + mod);
    auto sp = spectrogram(xs);
    auto fig = plotMatrix(sp);
    fig.save("freq.png");
}

auto melMatrix(size_t nfreq, size_t nbank)
{
    
}
