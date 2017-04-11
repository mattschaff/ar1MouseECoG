%% TEST CASE 1: Similar Frequency
    % 2 waves with some frequency (temporal & spatial)
    % check efficacy as lim(d_freq) --> 0
    % randomize other parameters (including wave type)

    
%% TEST CASE 2: Noisiness
    % N waves with varying noise
    % check efficacy as lim(sig/noise) --> 0
    % randomize other parameters (including wave type)

%% TEST CASE 3: Amplitude (Canonical Power Spectrum)
    % N waves with varying amplitudes
    % check efficacy as N --> inf
    % randomize other parameters (including wave type)

%% TEST CASE 4: Temporal Resolution (changes in parameters)
    % 1 wave with varying parameters over time (esp amplitude)
    % check efficacy as freq(parameter change) --> inf
    % randomize other parameters (including wave type)
    
%% TEST CASE 5: Phase Distribution
    % N identical waves with varying positions on a 2D grid
    % check efficacy as N --> inf
    % also check random phase distribution

%% EXPERIMENTAL CASE: Mouse ECoG Data