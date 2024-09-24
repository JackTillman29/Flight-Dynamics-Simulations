function out = parse_scpi_response(raw_response)
% this function reformats an instrument response into
% a cell column vector of strings.

% last character is integer "10"
% comma is 44
raw_response = raw_response';
% remove spaces and place commas
raw_response(raw_response == 32) = 44;
iCommas = find(raw_response == 44);
nCommas = length(iCommas);

nFields = nCommas + 1;

idx1 = 1;
out = cell(nFields,1);
for k = 1 : nFields
    if(k ~= nFields)
        idx2 = iCommas(k);
    else
        idx2 = length(raw_response)-1;
    end
    
    
    
    out{k} = char(raw_response(idx1:idx2-1));
    idx1 = idx2+1;

end