function [sigma]=lsq_lightcurve(theta,X,Y)

%lcurve= mu_max*(light_data./(K+light_data));
% lcurve=theta(1)*(1-exp(-theta(2)*X));
%lcurve=theta(1)+(theta(2)*log10(X)); %different functional form
lcurve=theta(1)+(1-exp(theta(2)*X));

sigma=Y-lcurve;

end

