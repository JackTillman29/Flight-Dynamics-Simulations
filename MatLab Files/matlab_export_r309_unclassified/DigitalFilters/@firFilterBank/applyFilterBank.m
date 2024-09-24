function y = applyFilterBank(FLTB,u,precision)

% =========================================
% BRUTE FORCE METHOD (COMMENTED)
% =========================================

%     n     = length(u);
%     nTaps = length(tapGains);
%     
%     % first output is after nTaps input
%     nOutputs = n - nTaps + 1;
%     
%     % define the first set of input data that will fill the stages
%     firstInput = 1;
%     lastInput  = firstInput + nTaps - 1;
%     
%     y = zeros(size(u,1),size(u,2));
%     
%     for k = 1 : nOutputs
%         % generate output
%         output     = sum(tapGains .* u(lastInput:-1:firstInput));
%         
%          % assign output
%         y(nTaps-1+k) = output;
%         
%         % update the input data set
%         firstInput = firstInput + 1;
%         lastInput  = lastInput  + 1;
%         
%     end

% =========================================
% FAST CONVOLUTION METHOD
% =========================================
if(~exist('precision','var'))
    precision = 'double';
end

if(strcmp(precision,'double'))
    y = zeros(length(u),FLTB.nChannels);
else
    y = single(zeros(length(u),FLTB.nChannels));
    u = single(u);
end

for k = 1 : FLTB.nChannels
    y(:,k) = conv(u,get(FLTB.filters{k},'tapGains'),'same');
end

end