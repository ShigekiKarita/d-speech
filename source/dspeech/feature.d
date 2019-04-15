/**
   Speech feature extraction module

   See also: "Chapter 5 Speech Input/Output," HTK book 3.5 alpha1, pp. 89 - 113.
*/
module dspeech.feature;


import numir.signal : hann;

/**
   Computes a magnitude 2-D spectrogram from a time-domain 1-D signal

Params:
    windowFun = window function
    xs = 1-D slice signal
    nperseg = (default 256) short-time frame width for each FFT segment
    noverlap = (default nperseg / 2) short-time frame overlapped length for each FFT segment

Returns: a magnitude 2D spectrogram
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
/// $(IMG dspeech.feature.spectrogram.png )
unittest
{
    import mir.ndslice : as, map, iota;
    import mir.math : sin, cos;
    import std.math : PI;
    import dspeech.plot : plotMatrix, docDir;

    auto time = iota(50000).as!double / 10e3;
    auto mod = 1e3 * map!cos(2.0 * PI * 0.5 * time);
    auto xs = map!sin(2.0 * PI * 3e3 * time  + mod);
    auto sp = spectrogram(xs);
    auto fig = plotMatrix(sp);
    fig.save(docDir ~ "dspeech.feature.spectrogram.png",
             cast(int) sp.length!1 * 2, cast(int) sp.length!0 * 2);
}

@nogc @safe nothrow pure melScale(double freq)
{
    import mir.math.common : log;
    return 1127.0 * log(1.0 + freq / 700.0);
}



/**
Computes Mel-scale filterbank matrix

Params:
    nFreq: the number of FFT bins
    nMel: the number of filterbank bins
    sampleFreq: sampling frequency of FFT signal
    lowFreq: lowest frequency threshold of the filterbank
    highFreq: highest (relative from Nyquist freq) frequency threshold of the filterbank

TODO:
    make this function @nogc
 */
@safe nothrow pure melMatrix(Dtype = double)(size_t nFreq, size_t nMel, Dtype sampleFreq, Dtype lowFreq= 20, Dtype highFreq = 0)
do
{
    import mir.ndslice : iota, as, sliced, map;
    import numir : zeros;
    const nyquist = 0.5 * sampleFreq;
    highFreq += nyquist;
    const lowMel = lowFreq.melScale;
    const highMel = highFreq.melScale;
    const deltaMel = (highMel - lowMel) / nMel;
    const deltaFreq = nyquist / nFreq;
    const binsMel = iota(nMel + 2).as!Dtype * deltaMel + lowMel;

    return iota(nFreq, nMel).map!(
        (i) {
            const iFreq = i / nMel;
            const iMel = i % nMel;
            const leftMel = binsMel[iMel];
            const centerMel = binsMel[iMel + 1];
            const rightMel = binsMel[iMel + 2];
            const mel = melScale(deltaFreq * iFreq);
            if (mel >= rightMel)
            {
                assert(iFreq > 0, "too large nMel you set");
                return cast(Dtype) 0;
            }
            if (mel > leftMel)
            {
                return (mel <= centerMel ? mel - leftMel : rightMel - mel) / deltaMel;
            }
            return cast(Dtype) 0;
        });
}

/// $(IMG dspeech.feature.mel.png )
unittest
{
    import mir.ndslice : transposed;
    import dspeech.plot : plotMatrix, docDir;
    auto m = melMatrix!double(200, 40, 16000.0);
    m.transposed.plotMatrix.save(
        docDir ~ "dspeech.feature.mel.png",
        cast(int) m.length!0 * 3, cast(int) m.length!1 * 9);
}
