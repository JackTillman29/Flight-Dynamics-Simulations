author: Jeff Hole (Booz Allen Hamilton)

These trivial oscillators are created to produces simple common waveforms with +/- 1 peaks.
These will find use as modulation sources for waveform generation, as well as utilities 


The core oscillator is the sawtooth oscillator (called "sawOsc.m") and is identified as
the core since the following oscillators are built up using this simple oscillator:
 - rectangular (difference of two saws, vary duty cycle by varying phase between saws)
 - triangle    (integrate rectangular oscillator)

More elaborate oscillators can be built, including band-limited anti-aliased oscillators.


TODO:
 - need to incorporate a DC-compensated triangle waveform for different rising and falling edges