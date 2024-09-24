function out = sampleAndHold(in,clk)
% output vector will be the same length as input vector but will have held
% values whenever a clock leading edge occurs.
%
% author: Jeff Hole (Booz Allen Hamilton), 7-24-2018

[rin,cin]   = size(in);
[rclk,cclk] = size(clk);

% detect clock leading edge.
dclk = diff([0 clk]);
clk = zeros(size(clk));
clk(dclk > 0) = 1;

idx = find(clk == 1);

out = zeros(size(in));

% one-dimensional signal
if(rin == 1 || cin == 1)
    for k = 1:length(idx)-1
        out(idx(k):idx(k+1)) = in(idx(k));
    end
    return
end

if(cin == cclk)
    for k = 1:length(idx)-1
        L = idx(k+1)-idx(k)+1;
        out(:,idx(k):idx(k+1)) = repmat(in(:,idx(k)),[1,L]);
    end
else
    for k = 1:length(idx)-1
        L = idx(k+1)-idx(k)+1;
        out(idx(k):idx(k+1),:) = repmat(in(idx(k),:),[L 1]);
    end
end


end