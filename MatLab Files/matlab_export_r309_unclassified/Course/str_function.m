function s = str_function(x)

fmt_str = '%10.1f';
if(imag(x) >= 0)
    s = [num2str(real(x),fmt_str) ' + j' num2str(abs(imag(x)),fmt_str)];
else
    s = [num2str(real(x),fmt_str) ' - j' num2str(abs(imag(x)),fmt_str)];
end


end