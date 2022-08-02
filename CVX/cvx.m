for i = 10
    
dfar = 23;
cvx_begin quiet
   variable dnear
   dual variables var1 var2
   minimize(-log(1+dnear)*dfar)%maximize secrecy
   subject to
      var1: dnear/dfar >= 0.1;
      var2: 1 <= dnear <= 10;
      var3: 10 <= dfar <= 20;
cvx_end

figure(1)
distancegap(i) = dfar - dnear
capacity(i) = log(1+dnear)*dfar
end
plot(distancegap,capacity,'k--')
