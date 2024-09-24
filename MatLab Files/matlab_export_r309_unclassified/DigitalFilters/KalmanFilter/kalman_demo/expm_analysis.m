function[eA]=expm_analysis(A,dt,n)

ratios = zeros(n,1);
ratios(1) = 1;
eA = zeros(size(A));
for k = 1:n
    term = (A^(k-1)) * dt^(k-1)/factorial(k-1);
    ratios(k) = norm(term);
    eA = eA + term;
end

figure;
plot(1:n,10*log10(ratios));
grid on;
xlabel('# of Taylor Series Terms');
ylabel('Term Magnitude (dB)');
end