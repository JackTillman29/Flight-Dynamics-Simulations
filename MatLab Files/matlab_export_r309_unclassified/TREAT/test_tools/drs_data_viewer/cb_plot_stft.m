v1 = single(data1(1:2:end))+1i*single(data1(2:2:end));
wo = WaveObj(v1,[],Fs);
figure;
wo.fast_stft(2^15,@rectwin,2^15);
ah = gca;
caxis([-20 50])
colorbar;
colormap(colormap_fade2black(256,.4));

v2 = single(data2(1:2:end))+1i*single(data2(2:2:end));
wo2 = WaveObj(v2,[],Fs);
figure;
wo2.fast_stft(2^15,@rectwin,2^15);
caxis([-20 50]);
ah = [ah gca];
linkaxes(ah,'xy');
colormap(colormap_fade2black(256,.4));