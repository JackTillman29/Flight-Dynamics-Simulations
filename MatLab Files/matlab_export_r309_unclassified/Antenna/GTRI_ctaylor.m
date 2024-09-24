	function Wr = GTRI_ctaylor(Rr,etadb,nb)
    % Rr = circular taylor npoints = sqrt((Xpos./ max(Xpos)).^2 +  (Ypos./max(Ypos)).^2);
    % etadb = transmit sum taylor weight = input.antenna.apertureTransmitSumTaylorWeightingValue = positive dB
    % nb = transmit sum taylor nbar = input.antenna.apertureReceiveSumTaylorWeightingnBar = integer (5?)
	Rmu = [1.2196699 2.2331302 3.2383154 4.2410628 5.2427643...
	 6.2439216 7.2447598 8.2453948 9.2458927 10.2462933];
	Rmu((nb+1):10) = []; fu0 = .202642367;

	pisq = pi * pi; nbm1 = nb - 1; eta = 10^(etadb / 20);

	a = log(eta + sqrt(eta * eta - 1)) / pi; asq = a * a; Rmusq = Rmu.^2;
	sgsq = Rmusq(nb) / (asq + (nb - .5).^2); Rmu(nb) = []; Rmusq(nb)= [];
	Rmu = pi * Rmu;

	Dnm = sgsq * (asq + ((1:nbm1) - .5).^2); Fu = zeros(1,nbm1);

	[RMUSQ DNM] = meshgrid(Rmusq,Dnm); Fnum = prod(1 - RMUSQ ./ DNM);
	[RMUMSQ RMUNSQ] = meshgrid(Rmusq); FDEN = (RMUNSQ - RMUMSQ) ./ RMUNSQ;
	Id = 1:nb:nbm1^2; FDEN(Id) = ones(size(Id)); Fden = prod(FDEN);

	Fm = -2*Fnum ./ (pisq *besselj(0,Rmu) .* Fden);
	[RR RMU] = meshgrid(Rr,Rmu); [RR FM] = meshgrid(Rr,Fm);
	Wr = fu0 * ones(size(Rr)) + sum(FM .* besselj(0,(RR .* RMU)));
	Wr = Wr / (fu0 + sum(Fm));