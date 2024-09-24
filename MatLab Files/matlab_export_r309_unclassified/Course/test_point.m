function test_point(tdIn,Fs,descr)
wo = WaveObj(tdIn,[],Fs);
plot(wo,1e-6,'MHz',0,'pwr')
wo_t_f_plot(wo,gcf,descr);
clear wo;
end