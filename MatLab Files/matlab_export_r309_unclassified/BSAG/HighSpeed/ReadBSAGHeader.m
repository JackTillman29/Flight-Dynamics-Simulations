function hdata = ReadBSAGHeader(file)
    % read header in as a string
    header_str = fileread(file);
    
    % regexp parse to get data
    % [match_strings,match_tokens]=regexp(header_str,'(\w+) (=) "(\S+)"','match','tokens');
    hdata.SIGNAL_NUM_BYTES = ParseBSAGHeaderString(header_str,'SIGNAL_NUM_BYTES');
    hdata.DATA_BITSIZE     = ParseBSAGHeaderString(header_str,'DATA_BITSIZE');
    hdata.DATA_BYTES_PER_UNIT       = ParseBSAGHeaderString(header_str,'DATA_BYTES_PER_UNIT');
    hdata.DATA_BYTE_ORDER           = ParseBSAGHeaderString(header_str,'DATA_BYTE_ORDER');
    hdata.DATA_COMPLEX              = ParseBSAGHeaderString(header_str,'DATA_COMPLEX');
    hdata.DATA_SAMPLE_RATE_MHZ      = ParseBSAGHeaderString(header_str,'DATA_SAMPLE_RATE_MHZ');
    hdata.LOADING_FRACTION          = ParseBSAGHeaderString(header_str,'LOADING_FRACTION');
    hdata.LOADING_MAX_VALUE         = ParseBSAGHeaderString(header_str,'LOADING_MAX_VALUE');
    hdata.LOADING_MIN_VALUE         = ParseBSAGHeaderString(header_str,'LOADING_MIN_VALUE');
    hdata.LOADING_OVERRANGE_COUNT   = ParseBSAGHeaderString(header_str,'LOADING_OVERRANGE_COUNT');
    hdata.LOADING_RANGE_MAX         = ParseBSAGHeaderString(header_str,'LOADING_RANGE_MAX');
    hdata.LOADING_RANGE_MIN         = ParseBSAGHeaderString(header_str,'LOADING_RANGE_MIN');
    hdata.MEASUREMENT_GAIN_DB       = ParseBSAGHeaderString(header_str,'MEASUREMENT_GAIN_DB');
    hdata.PATH_GAIN_DB              = ParseBSAGHeaderString(header_str,'PATH_GAIN_DB');
    hdata.SIGNAL_BB_CENTER_FREQ_MHZ = ParseBSAGHeaderString(header_str,'SIGNAL_BB_CENTER_FREQ_MHZ');
    hdata.SIGNAL_BB_FREQ_SPAN_MHZ   = ParseBSAGHeaderString(header_str,'SIGNAL_BB_FREQ_SPAN_MHZ');
    hdata.SIGNAL_NUM_SAMPLES        = ParseBSAGHeaderString(header_str,'SIGNAL_NUM_SAMPLES');
    hdata.SIGNAL_RF_CENTER_FREQ_MHZ = ParseBSAGHeaderString(header_str,'SIGNAL_RF_CENTER_FREQ_MHZ');
    hdata.RECORD_TIME               = ParseBSAGHeaderString(header_str,'RECORD_TIME');
        
end

