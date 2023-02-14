function y = yguessfun_steadystate(s, mss, tss, tns)
 y(1) = tss(s);
 y(2:3) = [0; 0];
 y(4) = mss(s);
 y(5) = tns(s);
 y(6:10) = [0; 0; 0; 0; 0];
end

