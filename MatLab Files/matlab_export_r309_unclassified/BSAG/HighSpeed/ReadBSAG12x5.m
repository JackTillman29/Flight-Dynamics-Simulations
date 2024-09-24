function [data,time,samplen]=ReadBSAG12x5(fid,headerData,mode,modeoption,startPacket,nPackets)
% Syntax: [data,time,samplen]=ReadBSAG12x5(fid,headerData,mode,modeoption,startPacket,nPackets)
% fid: MATLAB fopen'd file handle (as rb)
% headerData: Data structure from ReadBSAGHeader
% mode: 'full' or 'preview'
% modeoption: full has no option, use []. Preview option is # of packets to skip between outputs
% startPacket: Which packet is the first to be read (1 = first 8 byte packet)
% nPackets: # of packets to read
%
% -- output --
% data: real samples
% time: time vector in seconds
% samplen: sample #. in 'full' mode, this is actual samples. In 'preview' mode, it is packets (5 samples / packet)
% Keith Sawmiller, 2015

nPacketsSkip = startPacket-1;
nSamples = nPackets * 5;
nBytesPerPacket = headerData.DATA_BYTES_PER_UNIT;

frewind(fid);
if(startPacket ~= 0)
    fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
end

switch(lower(mode))
    case 'preview'
        previewPacketDecimation = modeoption;
        
        frewind(fid);
        fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
        
        % skip first four unused bits
        fread(fid,4,'*ubit1',0,'ieee-be');
        
        % read values in packet
        data.r = fread(fid, nPackets, '*bit12', previewPacketDecimation*64-12, 'ieee-be');
        temp = hilbert(single(data.r));
        data.i = int16(real(temp));
        data.q = int16(imag(temp));
        data.m = int16( sqrt( single( ...
            int32(data.i) .* int32(data.i) + int32(data.q) .* int32(data.q) ...
            ) ) );
        
        samplen = (0):(length(data.r)-1);
        samplen = samplen*previewPacketDecimation+(nPacketsSkip-1);
        time = single(5*samplen* (1/(1e6*headerData.DATA_SAMPLE_RATE_MHZ)));
        
%         figure(100);
%         plot(time*1e6,[data.r data.m]);
%         xlabel('Packet #');ylabel('Signed Counts');
%         title(['BSAG 12x5 Format Preview (Packet Decimation: ' num2str(previewPacketDecimation) ')']);
%         grid on;
%         axis tight;
%         %set(gca,'ButtonDownFcn',@cb_bsag_preview);
%         
%         % Define a context menu; it is not attached to anything
%         hcmenu = uicontextmenu;
%         % Define callbacks for context menu items that change linestyle
%         %hcb1 = '''';
%         assignin('base','temp_headerData',headerData);
%         assignin('base','temp_startPacket',startPacket);
%         assignin('base','temp_nPackets',nPackets);
%         assignin('base','temp_decimation',previewPacketDecimation);
%         assignin('base','temp_fid',fid);
%         
%         
%         item1 = uimenu(hcmenu,'Label','Resample','Callback', [ ...
%             'alim=get(gca,''XLim''); temp_startPacket = ceil(alim(1)); ReadBSAG12x5(temp_fid,temp_headerData,''preview'',ceil(diff(alim)/temp_nPackets),temp_startPacket,temp_nPackets);' ...
%             ]);
%         item2 = uimenu(hcmenu,'Label','Half Res','Callback', [ ...
%             'alim=get(gca,''XLim''); temp_startPacket = floor(alim(1)); temp_nPackets = temp_nPackets/2;ReadBSAG12x5(temp_fid,temp_headerData,''preview'',' num2str(floor(previewPacketDecimation*2)) ',temp_startPacket,temp_nPackets);' ...
%             ]);
%         item5 = uimenu(hcmenu,'Label','Max Preview Resample','Callback', [ ...
%             'alim=get(gca,''XLim''); temp_startPacket = floor(alim(1)); temp_nPackets = ceil(diff(alim));ReadBSAG12x5(temp_fid,temp_headerData,''preview'',1,temp_startPacket,temp_nPackets);' ...
%             ]);
%         item3 = uimenu(hcmenu,'Label','Full Res Window','Callback', [ ...
%             'alim=get(gca,''XLim''); temp_startPacket = floor(alim(1)); temp_nPackets = floor(alim(2)-alim(1)); [data,time,samplen]=ReadBSAG12x5(temp_fid,temp_headerData,''full'',''[]'',temp_startPacket,temp_nPackets);' ...
%             ]);
%         item4 = uimenu(hcmenu,'Label','Double Time','Callback', [ ...
%             '[data,time,samplen]=ReadBSAG12x5(temp_fid,temp_headerData,''preview'',' num2str(floor(previewPacketDecimation)) ',temp_startPacket,2*temp_nPackets);' ...
%             ]);
%         set(gca,'uicontextMenu',hcmenu);
        

        
    case 'full' % read as efficient interleaved format (no shifting or seeking)
        % skip first four unused bits
        fread(fid,4,'*ubit1',0,'ieee-be');
        disp(['Loading ' num2str(nSamples/1e6) 'M samples']);
        
        data1 = fread(fid, nPackets, '*bit12', 52, 'ieee-be');
        fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
        fread(fid,4+  12,'*ubit1',0,'ieee-be');
        
        data2 = fread(fid, nPackets, '*bit12', 52, 'ieee-be');
        fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
        fread(fid,4+2*12,'*ubit1',0,'ieee-be');
        
        data3 = fread(fid, nPackets, '*bit12', 52, 'ieee-be');
        fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
        fread(fid,4+3*12,'*ubit1',0,'ieee-be');
        
        data4 = fread(fid, nPackets, '*bit12', 52, 'ieee-be');
        fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
        fread(fid,4+4*12,'*ubit1',0,'ieee-be');
        
        data5 = fread(fid, nPackets, '*bit12', 52, 'ieee-be');
        fseek(fid,nBytesPerPacket*nPacketsSkip,'bof');
        
        % interleave samples
        data.r = [data1 data2 data3 data4 data5]';
        data.r = data.r(1:end);
        samplen = 0:(length(data.r)-1);
        time = single(samplen * (1/(1e6*headerData.DATA_SAMPLE_RATE_MHZ)));
        temp = hilbert(single(data.r));
        data.i = int16(real(temp));
        data.q = int16(imag(temp));
        data.m = int16( sqrt( single( ...
            int32(data.i) .* int32(data.i) + int32(data.q) .* int32(data.q) ...
            ) ) );
        

        figure;
        plot(1e6*time,[data.r.' data.m' ]);
        xlabel('Time (\mus)');ylabel('Signed Counts');
        title(['BSAG 12x5 Format Full']);
        grid on;
    otherwise
        error('mode must be either preview or full');
end


end