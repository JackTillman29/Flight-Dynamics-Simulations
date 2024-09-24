function out = GTRI_cbayliss(position,sll_db,nbar)
%!----------------------------------------------------------------------
%!
%! Function created in FORTRAN77 by Robert L. Howard
%! Adapted to matlab by Robert G. Hemphill Dec 2000
%!
%!----------------------------------------------------------------------
%Begin:  cbayliss
 
rmu=[0.5860670  1.6970509  2.7171939  3.7261370  4.7312271 ...
     5.7345205  6.7368281  7.7385356  8.7398505  9.7408945 ... 
    10.7417435 11.7424475 12.7430408 13.7435477 14.7439856 ...
    15.7443679 16.7447044 17.7450030 18.7452697 19.7455093];
 
sl1 = -1.*abs(sll_db);
sl2 = sl1.*sl1;
sl3 = sl2.*sl1;
sl4 = sl3.*sl1;
 
a = 0.3038753-0.05042922.*sl1-0.00027989.*sl2-0.00000343.*sl3-0.00000002.*sl4;
sl = [sl1 sl2 sl3 sl4];

nbarm1 = nbar - 1;
nbarp1 = nbar + 1;
sigma = rmu(nbarp1)./GTRI_ze(nbar,sl,a);
 
for k = 1:nbar
   pmu = pi.*rmu(k);
   coef = rmu(k).^2./besselj(1,pmu);
 
   prodm = 1;
   for m = 1:nbarm1
      fm = rmu(k)./(sigma.*GTRI_ze(m,sl,a));
      prodm = prodm.*(1-fm.^2);
   end
 
   prodn = 1;
   for n = 1:nbar
      if n ~= k
         fn = rmu(k)./rmu(n);
         prodn = prodn.*(1-fn.^2);
      end 
   end 
 
   b(k) = coef.*prodm./prodn;
end 
 
rr = sqrt(position(:,1).^2.+position(:,2).^2.);    
rr = rr/max(rr);
cphi = position(:,1)./rr;
 
%RLH: Why do we multiply by pi?
rr = pi * rr;

cbaywt = zeros(size(rr));
 
for k = 1:nbar
   rrmu = rr.*rmu(k);
   cbaywt = cbaywt + b(k).*besselj(1,rrmu);
end 
 
out = cbaywt.*cphi;
