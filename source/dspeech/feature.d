/**
   Speech feature extraction module

   See also:
   - "Chapter 5 Speech Input/Output," HTK book 3.5 alpha1, pp. 89 - 1
   - https://github.com/jameslyons/python_speech_features

   TODO: preemphasis and dither
*/
module dspeech.feature;

/** Compute Wave -> Spectrogram -> FBANK -> MFCC eventually

- preemphasised and dithered wave

$(IMG dspeech.feature.speech_wav.png )

- power spectrogram

$(IMG dspeech.feature.speech_spectrogram.png )

- log-Mel filterbank (FBANK)

$(IMG dspeech.feature.speech_fbank.png )

- Mel-frequency cepstral coefficients (MFCC)

$(IMG dspeech.feature.speech_mfcc.png )

 */
unittest
{
    import dspeech.plot : docDir, docWav, plotVector, plotMatrix;
    import dffmpeg : Audio;
    import mir.ndslice : map, sliced;

    auto raw = Audio!short().load(docWav);
    auto wav = raw.data.sliced.dither(1.0).preemphasis(0.97);
    wav.plotVector.save(docDir ~ "dspeech.feature.speech_wav.png", 400, 200);
    auto sp = wav.spectrogram(512, 256);
    sp.plotMatrix.save(docDir ~ "dspeech.feature.speech_spectrogram.png", 400, 200);
    auto fbanks = sp.toFbank(melMatrix(sp.length!0, 80, cast(double) raw.sample_rate));
    fbanks.plotMatrix.save(docDir ~ "dspeech.feature.speech_fbank.png", 400, 200);
    auto mfccs = fbanks.toMfcc;
    mfccs.plotMatrix.save(docDir ~ "dspeech.feature.speech_mfcc.png", 400, 200);
}



import numir.signal : hann;

/// Dither quantized wave
auto dither(S)(S wave, double stddev = 1.0)
{
    import mir.ndslice : map;
    import numir.random : RNG;
    import mir.random.variable: NormalVariable;

    return wave.map!(a => a + NormalVariable!double(0, stddev)(RNG.get));
}


/// Preemphasis wave signal
auto preemphasis(S)(S wave, double factor = 0.97)
{
    return wave[1 .. $] - factor * wave[0 .. $-1];
}

/**
   Computes a magnitude 2-D spectrogram from a time-domain 1-D signal

Params:
    windowFun = window function
    xs = 1-D slice signal
    nperseg = (default 256) short-time frame width for each FFT segment
    noverlap = (default nperseg / 2) short-time frame overlapped length for each FFT segment

Returns: a magnitude 2D spectrogram
*/
auto spectrogram(alias windowFun=hann, bool power=false, S)(S xs, size_t nperseg = 256, size_t noverlap = 128)
{
    import numir.signal : stft;
    import mir.math.common : sqrt;
    import mir.ndslice.topology : map;
    import mir.ndslice : transposed;

    auto stfs = stft!windowFun(xs, nperseg, noverlap);
    auto upper = stfs.transposed[nperseg / 2 .. $];
    static if (power)
    {
        return upper.map!(a => (a.re * a.re + a.im * a.im));
    }
    else
    {
        return upper.map!(a => sqrt(a.re * a.re + a.im * a.im));
    }
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

/**
Computes type-II DCT on 1-d signal:

$(MATH y[k] = 2 f \sum_{n=0}^{N-1} x[n] \cos \pi k \frac{2n + 1}{2 N} ),

where the scaling factor f = sqrt(1 / 4 N) if k = 0,  f = sqrt(1 / 2 N) otherwise,

Params:
    xs = input array

Returns:
    DCT transformed 1-d array

See_also: https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.dct.html
 */
auto dct(alias sumKind = "precise", S)(S xs)
{
    import std.math : PI;
    import mir.math.common;
    import mir.math.sum;
    import mir.primitives : DimensionCount;
    static assert(DimensionCount!S == 1);
    import mir.ndslice;
    const N = xs.length;
    const f0 = 2 * sqrt(0.25 / N);
    const f = 2 * sqrt(0.5 / N);
    return  iota(N).map!(k => (k == 0 ? f0 : f) * iota(N).map!(n => xs[n] * cos(PI * k * (2 * n + 1) / (2 * N))).sum!sumKind);
}


/// compare with scipy.fftpack.dct
unittest
{
    import mir.ndslice;
    import numir.testing : approxEqual;
    auto y = dct([4.0, 3.0, 5.0, 10.0].sliced);
    // python: scipy.fftpack.dct([4., 3., 5., 10.], type=2, norm="ortho")
    assert(approxEqual([11.        , -4.46088499,  3.        , -0.31702534].sliced, y.slice));
}


/// Computes log-Mel filterbank (FBANK) feature from spectrogram signal and Mel matrix
auto toFbank(S, M)(S spect, M melmat)
{
    import lubeck : mtimes;
    import mir.ndslice : map, transposed;
    import mir.math.common : log, fmax;

    return melmat.transposed.mtimes(spect).map!(x => fmax(x, typeof(x).epsilon).log);
}


/// Computes MFCC feature from log-Mel filterbank (FBANK) feature
auto toMfcc(S)(S fbanks, size_t nceps = 13, double lifter = 22)
{
    import numir : alongDim;
    import mir.ndslice : map, iota, ndiota, transposed, slice;
    import mir.math.common : sin;
    import std.math : PI;

    auto l = 1.0 + 0.5 * iota(nceps).map!(i => lifter * sin(PI * i / lifter));
    auto N = fbanks.length!1;
    return ndiota(nceps, N).map!(i => dct(fbanks[0 .. $, i[1]])[i[0] + 1] * l[i[0]]);
}
