classdef YATGObj
    
    %define properties
    properties
        
        Name; %descriptive name of transmition (string)
        NumWaveforms; %number of waveforms
        Fc; %system center frequency (Hz)
        TxFreq; %transmitted frequency (Hz)
        PRI; %pulse repitition interval (s)
        PW; %pulse width (s)
        Repeats; %number of times pulse gets repeated (integer)
        Attenuation; %attenuation (dB)
        LFMBW; % lfm bandwidth used to populate FMTable (Hz, + is upchirp, - is downchirp)
        FreqTable; % which frequency table to use (1-15 set by routine)
        ChipLength; %chip length (s)
        PhaseAngle; %phase angle, used to populate PHTable (0, 90, 180 or 270 deg)
        PhaseKey(:, :); %phase key, used to populate PHTable (ints of 1 or 0)
        PhaseTable; % which phase table to use (1-15 set by routine)
        Event; %not currently used (1-4, 0 disables)
        CW; %enable cw pulse (boolean 0 or 1)
        Blank; %used to put blank pulses in (boolean0 or 1)
        
    end
    
    methods

        %% create YATG object
        function this = YATGObj(filename)
            
            if(exist(filename) == 2)
                
                this = this.readPDWYCSV(filename);                   

            else
                
                warning('Creating Empty Flight Path Object!');
                
            end
                        
        end
        
        %% read "csv" formatted PDW file
        function this = readPDWYCSV(this, filename)
            
            % open file
            fid = fopen(filename);
            
            %preallocate & set idx
            freqidx = 0;
            phasidx = 0;
            pdwridx = 0;
            MaxChipLenght = 0;
            
            % read file
            rowdata = fgetl(fid);
            while ischar(rowdata)
                
                %if not a comment line
                if ~strcmp(rowdata(1), '#')
                    
                    %find the idx where the first comma is
                    idxnum = strfind(rowdata, ',') + 1;
                    
                    %if frequency table ($FMT)
                    if strcmp(rowdata(1:4), '$FMT')
                        
                        %set frequency table
                        freqidx = freqidx+1;
                        FreqData = str2num(rowdata(idxnum(1):end));
                        temp(freqidx).FreqTable = FreqData(2:3);
                        
                    %elseif phase table
                    elseif strcmp(rowdata(1:4), '$PHT')
                        
                        %set phase table
                        phasidx = phasidx+1;
                        PhaseData = str2num(rowdata(idxnum(1):end));
                        temp(phasidx).PhaseTable = PhaseData(1:end);
                        MaxChipLenght = max(MaxChipLenght, PhaseData(1));
                        
                    %elseif PDW record
                    elseif ~isempty(str2num(rowdata))
                                              
                        %set pwd table
                        pdwridx = pdwridx + 1;
                        PWDData = str2num(rowdata);
    
                        %set waveform parameters
                        this.TxFreq(pdwridx) = PWDData(1) * 1e6;
                        this.PRI(pdwridx) = PWDData(2) * 1e-6;
                        this.PW(pdwridx) = PWDData(3) * 1e-6;
                        this.Repeats(pdwridx) = PWDData(4);
                        this.Attenuation(pdwridx) = PWDData(5);
                        this.FreqTable(pdwridx) = PWDData(6);
                        this.LFMBW(pdwridx) = 0;
                        if PWDData(6) > 0
                            this.LFMBW(pdwridx) = diff(temp(PWDData(6)).FreqTable) * 1e6; 
                        end
                        this.PhaseTable(pdwridx) = PWDData(7);
                        this.ChipLength(pdwridx) = 0;
                        this.PhaseAngle(pdwridx) = 0;
                        this.PhaseKey(pdwridx, :) = 0;
                        if PWDData(7) > 0
                            this.ChipLength(pdwridx) = temp(PWDData(7)).PhaseTable(3) .* 1e-6;
                            this.PhaseAngle(pdwridx) = max(temp(PWDData(7)).PhaseTable(2:2:end));
                            if (this.PW(pdwridx) / this.ChipLength(pdwridx) < MaxChipLenght)
                                this.PhaseKey(pdwridx, 1:MaxChipLenght) = [temp(PWDData(7)).PhaseTable(2:2:end) zeros(1, MaxChipLenght - (this.PW(pdwridx) / this.ChipLength(pdwridx)))] ./ this.PhaseAngle(pdwridx);
                            else
                                this.PhaseKey(pdwridx, 1:MaxChipLenght) = temp(PWDData(7)).PhaseTable(2:2:end) ./ this.PhaseAngle(pdwridx);
                            end                        
                        end
                        this.Event(pdwridx) = PWDData(8);
                        this.CW(pdwridx) = PWDData(9);
                        this.Blank(pdwridx) = PWDData(10);
                        
                    %endif (frequency table ($FMT))
                    end
                    
                %endif (not a comment line)
                end
                
                %grab the next line
                rowdata = fgetl(fid);
                
            end
            
            this.NumWaveforms = pdwridx;
            
            %% close file
            fclose(fid);
            
        end
               
        %% write "csv" formatted PDW file  
        function writePDWYCSV(this, filename)
            
            % writes a PDW to a csv file in the format for the YATG
            % inputs:
            %   fid: file pointer
            %   this: waveform descriptor word structure containing the following fields:
            %       Name (string)
            %       NumWaveforms (integer)
            %       Fc (Hz) - System Center Frequency
            %       TxFreq (Hz) - Transmitted Frequency
            %       PRI (s)
            %       PW (s)
            %       Repeats (integer)
            %       Attenuation (dB)
            %       FMTable (Frequency Modulation, starting offset & Ending offset Hz)
            %       Phase Modulation (Phase Code)
            %       Events (Not Used Currently)
            %       CW Override (Not Used Currently)
            %       Blank Override (boolean)
            %
            %
            %##########################################################################
            %
            % CSV2YBF Inpute File Format:
            %
            % FM Table Record:
            %
            %   $FMTn, FREQ, START, STOP, TIME, HOLD
            %
            %     Field  Name      Description
            %     -----  --------  ----------------------
            %         1  "$FMTn"   Token (n = 1-7)
            %         2  FREQ      Frequency (MHz, 0 disables)
            %         3  START     Start offset (MHz)
            %         4  STOP      Stop offset (MHz)
            %         5  TIME      Sweep time (usec)
            %         6  HOLD      Hold Enable (0,1)
            %
            % Phase Table Record:
            %
            %   $PHTn, N, [PHASE1, TIME1, PHASE2, TIME2, ... , PHASEN, TIMEN]
            %
            %     Field  Name      Description
            %     -----  --------  ----------------------
            %         1  "$FMTn"   Token (n = 1-7)
            %         2  N         Number of chips (1-32, 0 disables)
            %         3  PHASE1    Chip 1 phase (0, 90, 180, 270 deg)
            %         4  TIME1     Chip 1 time (usec)
            %         *     *            *            *
            %    N +  2  PHASEN    Chip N phase (0, 90, 180, 270 deg)
            %    N +  3  TIMEN     Chip N time (usec)
            %
            % PDW Record:
            %
            %   FREQ, PRI, PW, COUNT, ATTEN, FM, PHASE, EVENT, CW, BLANK
            %
            %     Field  Name      Description
            %     -----  --------  ----------------------
            %         1  FREQ      Frequency (MHz)
            %         2  PRI       Pulse Repitition Interval (usec)
            %         3  PW        Pulse Width (usec)
            %         4  COUNT     Pulse Count
            %         5  ATTEN     Attenuation (dB)
            %         6  FM        FM Table Slot (1-15, 0 disables)
            %         7  PHASE     Phase Table Slot (1-15, 0 disables)
            %         8  EVENT     Event (1-4, 0 disables)
            %         9  CW        CW Enable (boolean)
            %         10 BLANK     Blank Enable (boolean)
            %
            % Notes:
            %   Fields must be separated by a comma.
            %   White space a the beginning & end of lines is ignored.
            %   Blank lines are ignored.
            %   Lines that start with '#' are treated as comments.
            %   Lines starting with '$' are special tokens.
            %   Undefined table records default to zero (disabled).
            %   All other lines are treated as PDW records.
            %
            %##########################################################################
            
            %% open file
            fid = fopen(filename, 'w');
            
            %% write header informatin (comments onlY)
            fprintf(fid, '%s', '################################');
            fprintf(fid, '%s\n', '#################################');
            fprintf(fid, '%s', '# This file has been automatically generated');
            fprintf(fid, '%s\n', ' to be used with YATG');
            fprintf(fid, '%s\n', ['# Name of Waveform: ' this.Name]);
            fprintf(fid, '%s', '################################');
            fprintf(fid, '%s\n', '#################################');
            fprintf(fid, '%s\n', '#');
            
            %% wite frequency table
            fprintf(fid, '%s\n', '# Frequency Table Record:');
            fprintf(fid, '%s\n', '#');
            fprintf(fid, '%s\n', '# $FMTn, FREQ (MHz), START (MHz), STOP (MHz), TIME (us), HOLD');
            
            %get unique lfm bw & pw list
            [uniquepwbw idxm] = unique([this.PW' this.LFMBW'], 'stable', 'rows');

            cnt = 1;
            %loop over the unique lfm bw & pw pairs
            for k = 1:size(uniquepwbw , 1)
                
                if (abs(uniquepwbw(k, 2)) > 0)
                    
                    %throw error if more than 15 unique frequency tables
                    if (cnt > 15)

                        error(' Error Maximum of 15 Unique Frequency Tables')

                    end

                    fprintf(fid, '%s', ['$FMT', num2str(cnt) ',']);
                    %center frequency (MHz)
                    fprintf(fid, '%s', [num2str(this.Fc(idxm(k)) .* 1e-6) ',']);
                    %start offset frequency (MHz)
                    fprintf(fid, '%s', [num2str(-1*this.LFMBW(idxm(k)) / 2 * 1e-6) ',']);
                    %end offset frequency (MHz)
                    fprintf(fid, '%s', [num2str(1*this.LFMBW(idxm(k)) / 2 * 1e-6) ',']);
                    %sweep time (usec)
                    fprintf(fid, '%s', [num2str(this.PW(idxm(k)) * 1e6) ',']);
                    %assuming hold always 0 for now
                    fprintf(fid, '%s\n', '0');

                    %increment count since valid lfm pulse
                    cnt = cnt + 1;
                    
                end
                                    
            end
            
            uniquepwbw = uniquepwbw(abs(uniquepwbw(:,2)) > 0, :);
            %loop over the # of waveforms
            for k = 1:this.NumWaveforms
                
                if (abs(this.LFMBW(k)) > 0)
                    
                    %loop over the unique lfm bw & pw pairs
                    for l = 1:size(uniquepwbw , 1)

                        if (this.PW(k) == uniquepwbw(l, 1) && this.LFMBW(k) == uniquepwbw(l, 2))
                            
                            %set what frequency table to use
                            this.FreqTable(k) = l;
                            
                        end
                        
                    end
                    
                else
                    
                    %set what frequency table to 0
                    this.FreqTable(k) = 0;
                    
                end
                
            end
            
            fprintf(fid, '%s\n', '#');
            
            %% wite phase table
            fprintf(fid, '%s\n', '# Phase Table Record:');
            fprintf(fid, '%s\n', '#');
            fprintf(fid, '%s\n', '# $PHTn, N, [PHASE1 (deg), TIME1(us), PHASE2, TIME2, ... , PHASEN, TIMEN]');
            
            cnt = 1;
            %loop over the # of waveforms
            for k = 1:this.NumWaveforms
                
                %throw error if more than 15 unique phase tables
                if (cnt > 15)
                    
                    error(' Error Maximum of 15 Unique Frequency Tables')
                    
                end
                
                if (this.ChipLength(k) > 0)
                    
                    fprintf(fid, '%s', ['$PHT', num2str(cnt) ',']);
                    %# of chips
                    fprintf(fid, '%s', [num2str(this.PW(k) / this.ChipLength(k)) ',']);
                    %loop over # of chips
                    for l = 1:(this.PW(k) / this.ChipLength(k))
                        fprintf(fid, '%s', [num2str(this.PhaseKey(l, k) * this.PhaseAngle(k)) ',']);
                        if (l ~= (this.PW(k) / this.ChipLength(k)))
                            fprintf(fid, '%s', [num2str(this.ChipLength(k) * 1e6) ',']);
                        else
                            fprintf(fid, '%s\n', [num2str(this.ChipLength(k) * 1e6)]);
                        end
                    end
                    
                    %set what phase table
                    this.PhaseTable(k) = cnt;
                    
                    %increment count since valid lfm pulse
                    cnt = cnt + 1;
                    
                else
                    
                    %set what phase table
                    this.PhaseTable(k) = 0;
                    
                end
                
            %end loop over the # of waveforms    
            end
            
            fprintf(fid, '%s\n', '#');
            
            %% write PDW record
            fprintf(fid, '%s\n', '# PDW Record:');
            fprintf(fid, '%s\n', '#');
            fprintf(fid, '%s\n', '# FREQ (MHz), PRI (us), PW (us), COUNT, ATTEN (dB), FM, PHASE, EVENT, CW, BLANK');
            
            %% close file
            fclose(fid);
            
            OutMat = [this.TxFreq .* 1e-6;this.PRI .* 1e6;this.PW .* 1e6;this.Repeats;...
                this.Attenuation;this.FreqTable;this.PhaseTable;this.Event;...
                this.CW;this.Blank];
            
            %write to file
            dlmwrite(filename, OutMat','-append', 'delimiter',',')
            
        end
        
        %% create signals
        function [tdSignal,Time] = createPDWSignals(this, Fs, FreqOffset, varargin)
                        
            %calculate Ts
            Ts = 1 / Fs;
            
            %calculate time vector & samples needed
            Time = 0:Ts:sum(this.PRI(1:this.NumWaveforms) .* this.Repeats(1:this.NumWaveforms)) - Ts;
            tdSignal = zeros(size(Time));
            
            idxlow = 1;
            %loop over the number of waveforms
            for k = 1:this.NumWaveforms
 
                %calculate high index
                idxhigh = idxlow + length(0:Ts:this.PRI(k)*this.Repeats(k)-Ts) - 1;
                
                %if not a blanked period
                if (this.Blank(k) ~= 1)
                    
                    %calculate dwell time
                    DwellTime = 0:Ts:this.Repeats(k)*this.PRI(k)-Ts;
                    N_Dwell = size(DwellTime, 2);
                    N_PW = size(0:Ts:this.PW(k)-Ts, 2);
                    N_PRI = size(0:Ts:this.PRI(k)-Ts, 2);
                
                    %calculate center freq
                    FreqCenter = this.TxFreq(k) - FreqOffset;
                    
                    %calculate pulsed on/off
%                     OnOffPulse = PWM_Signal(this.PW(k), this.PRI(k), Ts, N_Dwell);
                    
                    %create carrier waveform
                    CarrierSig = SimpleLfm(FreqCenter, FreqCenter, DwellTime(end)+Ts, Ts);
                    
                    %set modulated waveform to 0
                    ModulatedSig = zeros(size(CarrierSig));

                    %if a modulated waveform (lfm)
                    if (abs(this.LFMBW(k)) > 0)
                        
                        %create waveform
                        TempSig = SimpleLfm(-this.LFMBW(k)/2.0, this.LFMBW(k)/2.0, this.PW(k), Ts);

                        %loop over the number of repeats
                        idx = 1:N_PW;
                        for l = 1:this.Repeats(k)
                            %set modulated sig
                            ModulatedSig(idx) = TempSig;
                            idx = idx + N_PRI;
                        end
                    
                    %elseif a modulated waveform (bpsk)
                    elseif (this.ChipLength(k) > 0)
                        
                        %set bitkey & pulse time
                        BitKey = this.PhaseKey(k, :);
                        BitKey = BitKey(1:this.PW(k)/this.ChipLength(k));
                        PulseTime = 0:Ts:this.PW(k)-Ts;

                        %create waveform
                        TempSig = exp(1i .* (SimplePSK(BitKey, PulseTime) .* pi));
                            
                        %loop over the number of repeats
                        idx = 1:N_PW;
                        for l = 1:this.Repeats
                            %set modulated sig
                            ModulatedSig(idx) = TempSig;
                            idx = idx + N_PRI;
                        end
                        
                    %else a unmodulated waveform    
                    else
                        
                        %loop over the number of repeats
                        idx = 1:N_PW;
                        for l = 1:this.Repeats(k)
                            %set sig
                            ModulatedSig(idx) = 1;
                            idx = idx + N_PRI;
                        end
                        
                    end

                    %create samples for the dwell
                    tdSignal(idxlow:idxhigh) = CarrierSig .* ModulatedSig;
                
                end
 
                %increase indxlow
                idxlow = idxhigh+1;       

            %end loop over the # of waveforms    
            end
            
            %if only certain amount of signal requested
            if (nargin == 4)
                idx = find(Time > varargin{1}, 1);
                Time = Time(1:idx);
                tdSignal = tdSignal(1:idx);
            end   
            
        end
        
    end
    
end