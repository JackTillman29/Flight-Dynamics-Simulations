function [dataArray] = AntennaPatternReader(fileName)
%%%%%%%%%%%%%%
%ESAMS Antenna Pattern Reader
%AMC January 2011
%Note: does not handle nested parentheses

fid = fopen(fileName);

if(fid > 0)

    processingComment = 0;
    j = 1;
    line = fgetl(fid);
    startComment = 0;

    while(line > -1 )
        
        wordArray = textscan(line, '%s');
        for i = 1:length(wordArray{1})
            wordCell = wordArray{1}(i);
            currentWord = wordCell{1};
            openParen =  strfind(currentWord,'(');
            closeParen =  strfind(currentWord,')');
            if( openParen > 0)
                startComment = true;
                currentWord = currentWord(1:(openParen-1)); 
                if(closeParen > openParen)
                    startComment = 0;
                end
            
            elseif(closeParen > 0)
                processingComment = false;
                currentWord = currentWord((closeParen+1):end);
            end
            
            if(processingComment == false )
                if(startComment == true)
                    processingComment = true;
                    startComment = false;
                end

                if(isempty(currentWord) == 0)
                    if(isnan(str2double(currentWord))== false)
                        dataArray(j) = str2double(currentWord);
                        j = j+1;
                    end
                end
            end
        end
        processingComment = 0;
        line = fgetl(fid);
        
                
    end
    
   
        
end
%
fclose(fid);