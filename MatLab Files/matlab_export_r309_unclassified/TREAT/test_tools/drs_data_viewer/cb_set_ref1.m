function cb_set_ref1( varargin )
    evalin('base','ref_channel=1;');
    disp('Setting reference to channel 1');
    evalin('base', ...
        'set(h.a1,''Color'',''w'',''YColor'',''r'',''XColor'',''r'');');
    evalin('base', ...
        'set(h.a2,''Color'',''w'',''YColor'',''k'',''XColor'',''k'');');
end

