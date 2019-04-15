# D-SPEECH

[![CircleCI](https://circleci.com/gh/ShigekiKarita/d-speech.svg?style=svg)](https://circleci.com/gh/ShigekiKarita/d-speech)
[![codecov](https://codecov.io/gh/ShigekiKarita/d-speech/branch/master/graph/badge.svg)](https://codecov.io/gh/ShigekiKarita/d-speech)

WIP. speech processing (especially, speech recognition) toolkit for D

## feature

- minimal
- efficient
- faster

## roadmap

- dsp.feature module
- dsp.learn.acoustic module
- dsp.wfst module

## structure

- source/dsp/
  - feature: speech feature (FBANK, MFCC) extraction module
  - learn: machine learning module
    - acoustic: acoustic model module
      - hmm: hidden Markov model based sequence transducer
      - gmm: Gaussian mixture model based acoustic model module
      - neural: nueral network based acoustic model module
    - language: language model module
      - ngram: n-gram based language model module
      - neural: neural network based language model module
  - wfst: weighted finite state transducer module
    - container: sequence acceptor and transducer containers
    - algorithm: wfst related algorithms

- bin/
  - dsp-feature: speech feature (FBANK, MFCC) extractor command
