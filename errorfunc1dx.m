function errorval=errorfunc1dx(params, xval, yval)
    %Y=params(1)*xval+params(2);
    
m=params;
%offset=params(2);
Y=1./(xval.^m);
    %Y=exp(-xval.*multiplier);
    %Y=xval.*m+b;
    errorval= mean(((abs(Y'-yval))));
    
%     figure(666)
%     clf
%     hold on;
%     plot(xval,Y,'r')
%     plot(xval,yval)
%     axis([1 6000 -0.1 0.2]);
%     pause(0.1)
end
