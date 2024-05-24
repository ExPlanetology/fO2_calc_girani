function x2= NewtRaphZero(func,x0,epsilon1)
%This Newton Raphson Scheme does not allow the dependent variable (x2) to
%be less than zero.

epsilon2=1e-12;
xa=x0;
xaplus=xa+epsilon2;
xaminus=xa-epsilon2;
ya=func(xa);
yaplus=func(xaplus);
yaminus=func(xaminus);

counter=1;

while ((ya>epsilon1) && (counter<300))
    dfunc=(yaplus-yaminus)/(xaplus-xaminus);
    xa=xa-ya/dfunc;
    if (xa<0)
        if counter>50
            xa=1e-8;
            %if loop is having trouble converging, it may be because
            %equilibrium x2 is very small and needs to be set so.
        else
            xa=0.000001;
        end
    end
    ya=func(xa);
    xaplus=xa+epsilon2;
    xaminus=xa-epsilon2;
    yaplus=func(xaplus);
    yaminus=func(xaminus);
    counter=counter+1;
end
if (counter<300)
    x2=xa;
else
    %If N-R scheme doesn't converge, set x2 to a standard small value.
    x2=0.0001;
end
end