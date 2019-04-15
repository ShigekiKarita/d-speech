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
 */
auto melMatrix(Dtype = float)(size_t nFreq, size_t nMel, Dtype sampleFreq, Dtype lowFreq= 20, Dtype highFreq = 0)
in
{
    // import mir.math.common : log2;
    // assert(log2(nFreq) % 1 == 0);
}
do
{
    import mir.ndslice : iota, as, sliced;
    import numir : zeros;
    const nyquist = 0.5 * sampleFreq;
    highFreq += nyquist;
    const lowMel = lowFreq.melScale;
    const highMel = highFreq.melScale;
    const deltaMel = (highMel - lowMel) / nMel;
    const deltaFreq = nyquist / nFreq;
    const binsMel = iota(nMel + 2).as!Dtype * deltaMel + lowMel;

    auto ret = zeros!Dtype(nFreq, nMel);
    foreach (iMel; 0 .. nMel)
    {
        const leftMel = binsMel[iMel];
        const centerMel = binsMel[iMel + 1];
        const rightMel = binsMel[iMel + 2];
        size_t startFreq = 0;
        Dtype[] binsFreq;
        foreach (iFreq; 0 .. nFreq)
        {
            const mel = melScale(deltaFreq * iFreq);
            if (mel >= rightMel) break;
            if (mel > leftMel)
            {
                if (startFreq == 0)
                {
                    startFreq = iFreq;
                }
                Dtype weight;
                if (mel <= centerMel)
                {
                    weight = (mel - leftMel) / (centerMel - leftMel);
                }
                else
                {
                    weight = (rightMel - mel) / (rightMel - centerMel);
                }
                binsFreq ~= weight;
            }
        }
        assert(binsFreq.length > 0, "too large nMel you set");
        ret[startFreq .. startFreq + binsFreq.length, iMel] = binsFreq.sliced;
    }
    return ret;
}

///
unittest
{
    import mir.ndslice;
    import dspeech.plot : plotMatrix;
    auto m = melMatrix!double(256, 40, 16000.0);
    m.transposed.plotMatrix.save(
        "mel.png", cast(int) m.length!0 * 3, cast(int) m.length!1 * 9);
}
