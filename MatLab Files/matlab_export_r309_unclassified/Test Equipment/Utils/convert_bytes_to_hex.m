function str = convert_bytes_to_hex(u)


if(size(u,1) > size(u,2))
    u = u';
end

x = dec2hex(u);

str = zeros(1,2*length(u(:)),'uint8');

str(1:2:end) = x(:,1);
str(2:2:end) = x(:,2);

str = char(str);

end