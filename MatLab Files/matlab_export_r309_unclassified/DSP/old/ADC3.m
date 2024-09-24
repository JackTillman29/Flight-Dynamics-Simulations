function [adc_out,adc_out_counts] = ADC3(analog_in,nbits,varargin)
% [adc_out,adc_out_counts] = ADC3( analog_in, nbits, varargin )
% varargin{1} = min volts, {2} = max volts
% NOTE: if not provided, these are assumed to be min & max of input
    if(nargin == 2)
        adc_min = min(analog_in);
        adc_max = max(analog_in);
    else
        adc_min = varargin{1};
        adc_max = varargin{2};
    end
    adc_span = adc_max - adc_min;
    disp(['ADC Span: ' num2str(adc_span) 'V']);
    adc_vperb = adc_span / (2^nbits);
    disp(['ADC mV/bit: ' num2str(1e3*adc_vperb) 'mV']);
    
    % compute fraction of span
    m = 1.0 / adc_vperb;
    x1 = adc_min;
    y1 = 0;
    
    counts = floor( m * (analog_in - x1) + y1);
    counts = max(counts,0);
    counts = min(counts,2^nbits-1);
    adc_out_counts = counts;
    adc_out = adc_out_counts * adc_vperb+x1+0.5*adc_vperb;
end