function resi=fit_powercurve(theta,x,y)

test=theta(1)*x.^(theta(2));
resi=test-y; %model - observed 

end