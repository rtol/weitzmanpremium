function RF = RadiativeForcing(CO2,CH4,N2O,SF6,CFC11,CFC12,SO2,O3)
%function RF = RadiativeForcing(CO2,CH4,N2O,SF6,CFC11,CFC12,SO2,O3)
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.0-matlab-global
%
%This function is part of FUND 4.0 MG
%It returns radiative forcing
%
%Richard Tol, 8 August 2014
%This code is protected by the MIT License

global CO2RF CO21750 CH41750 N2O1750 CH4RF N2ORF SF6RF CFC11RF CFC12RF SdirRF SindRF 

RF = CO2RF*log(CO2/CO21750) + CH4RF*(CH4^0.5 - CH41750^0.5) + N2ORF*(N2O^0.5 - N2O1750^0.5) - CH4N2Oint(CH4,N2O1750)  - CH4N2Oint(CH41750,N2O) + 2*CH4N2Oint(CH41750,N2O1750) + SF6RF*SF6 +CFC11RF*CFC11 +CFC12RF*CFC12 - SdirRF*SO2 - SindRF*log(1+SO2*2/34.4)/log(1+14.6/34.4) + O3;

end

function RFi = CH4N2Oint(M,N);

global  CH4N2ORF CH4N2Op1 CH4N2Op2 CH4N2Op3 CH4N2Op4 CH4N2Op5

RFi = CH4N2ORF*(1+CH4N2Op1*M^CH4N2Op3*N^CH4N2Op3+CH4N2Op2*M^CH4N2Op4*N^CH4N2Op5);

end