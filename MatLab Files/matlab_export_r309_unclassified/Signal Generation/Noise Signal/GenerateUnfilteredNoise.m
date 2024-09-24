function y = GenerateUnfilteredNoise(N,Fs_Hz,PSD_wPerHz,method)
% Generates a vector uniform noise the specified power spectral density
% The analyst should pass this output through the appropirate IF bandpass
% filter to capture the correct in-band noise power.
%
% Syntax: y = GenerateUnfilteredNoise(N, ... % # of samples
%                                     Fs_Hz, ... % Sample Rate
%                                     PSD_wPerHz,method); % Pwr Spectral Density
% method: 'uniform_phase' = old method of e^j(rand phase)
%         'gaussian_iq' = randn(i) + randn(q)
% K. Sawmiller 2017

defaultPSD = 1.0 / Fs_Hz;
scaleFactor = sqrt(PSD_wPerHz / defaultPSD);

if(~exist('method'))
    method = 'uniform_phase';
end

switch lower(method)
    case 'uniform_phase'
        y = scaleFactor .* exp(1i*2*pi*rand(1,N,'single'));
    case 'gaussian_iq'
        % explanation for sqrt(0.5) term:
        %    0.5: half watt in the real part and imag part being
        %         summed together gives 1 watt
        %    sqrt(): convert watts to voltage
        y = sqrt(0.5) * scaleFactor .* (randn(1,N,'single')+1i*randn(1,N,'single'));
    otherwise
        error('Unknown noise format');
end


end