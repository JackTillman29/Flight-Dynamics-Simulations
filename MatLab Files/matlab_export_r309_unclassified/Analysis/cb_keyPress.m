function cb_keyPress(hfig,hkeydata)

currentChar = get(hfig,'CurrentCharacter');

% disp(['current character pressed:',currentChar])

support_strings{1} = 'p - dx,dy,m between two points';
support_strings{2} = 'f - compute FFT of waveforms on plot';
support_strings{3} = '    |-> f - compute FFT over entire signal duration';
support_strings{4} = '    |-> w - compute FFT over a section of the signal';

switch currentChar
    case 'p'
        getp2p;
    case 'f'
        set(hfig,'Interruptible','off')
        hmsg = msgbox({'f - FFT of entire data';'w - Select a portion of the signal to compute FFT'})
        k = waitforbuttonpress;
        if(k == 1) % key press
            k = get(gcf,'CurrentCharacter');
        end
        disp(['you entered:',k])
        close(hmsg)
        switch k
            case 'f'
                cb_WaveFftStruct;
            case 'w'
                zoom xon
                waitfor(gca,'XLim');
                zoom off
                timeWindow = get(gca,'XLim');
%                 timeWindow(1) = box(1)
%                 timeWindow(2) = box(2)
                cb_WaveFftStruct(timeWindow);
        end
        set(hfig,'Interruptible','on')
        
    case 'h'
        hsupp = msgbox(support_strings)
        
    otherwise
%         disp('unsupported...add your own!')
%         hsupp = msgbox(support_strings);
%         set(hsupp,'Name','HELP')
end



end