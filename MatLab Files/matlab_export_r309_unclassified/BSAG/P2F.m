function [packet_idx,sample_idx] = P2F(pidx,preview_hop)
    packet_idx = (pidx-1) .* preview_hop + 1;
    sample_idx = (packet_idx-1) .* 5 + 1;
end