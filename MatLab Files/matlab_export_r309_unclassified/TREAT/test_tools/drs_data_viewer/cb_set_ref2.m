function cb_set_ref2( varargin )
    evalin('base','ref_channel=2;');
    disp('Setting reference to channel 2');
    evalin('base', ...
        'set(h.a2,''Color'',''w'',''YColor'',''r'',''XColor'',''r'');');
    evalin('base', ...
        'set(h.a1,''Color'',''w'',''YColor'',''k'',''XColor'',''k'');');
end