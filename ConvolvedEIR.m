function retVec = ConvolvedEIR(f0, FWHM, dt, len)
%ConvolvedEIR Returns the result of convolving the electro-mechanical 
%impulse response (EIR) function with the derivative of the delta function,
%i.e. convolution of delta' and EIR.
%
% NOTES: 
%   This form of the EIR is taken from equation (6) in the paper
%   "Effects of light scattering on optical-resolution photoacoustic
%   microscopy" by Liu et al. (2012).
%
% INPUTS:
%   f0      - The center frequence of the transducer.
%   FWHM    - Full width at half-maximum spot size.
%   dt      - Time spacing between points.
%   len     - Length in terms of number of points.
%
% OUTPUTS:
%   retVec  - The sampled EIR

a = sqrt(2*log(2))/(pi*FWHM);
n = floor(len/2)*dt;
t = [-n:dt:n];
retVec = 2*pi*f0*cos(2*pi*f0*t).*exp(-t.^2./(2*a^2)) - t/(2*a^2).*sin(2*pi*f0*t).*exp(-t.^2/(2*a^2));

end
