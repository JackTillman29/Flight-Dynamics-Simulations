function adc_out = ADC2(analog_in,nbits,varargin)
    % Syntax: ADC2(analog_in,nbits,varargin=voltage span)
    % Note: code does not care if signal is AC coupled. Only about peak to peak voltage
   
    if(nargin == 2)
        adc_span = max(analog_in) - min(analog_in);
    else
        adc_span = varargin{1};
    end
    
    disp(['ADC Span: ' num2str(adc_span) 'V']);
    adc_vperb = adc_span / (2^nbits);
    disp(['ADC mV/bit: ' num2str(1e3*adc_vperb) 'mV']);
    
    miny = min(analog_in);
    
    
    adc_out = (analog_in-miny)-mod(analog_in-miny+0.5*adc_vperb,adc_vperb)+miny+0.5*adc_vperb;
    
end