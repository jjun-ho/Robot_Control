function dydt = One_link(t,y)

global Iz1 L1 g m1 r1 tq1

th1 = y(1);
dth1 = y(2);

dydt = [dth1;(tq1-g.*m1.*r1.*cos(th1))./(Iz1-L1.^2.*m1+L1.*m1.*r1.*2.0)];
end