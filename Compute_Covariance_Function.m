function C_func = Compute_Covariance_Function(Nx, Ny, r)
%COMPUTE_COVARIANCE_FUNCTION  Compute covariance function.
%
%   C_func = COMPUTE_COVARIANCE_FUNCTION(Nx, Ny, r) computes the covariance
%   function according to Equation (14) in [1].
%
% REFERENCES:
%   [1] J. C. Dainty. The statistics of speckle patterns. Progress in
%   Optics XIV, 1976.
        
[X, Y] = meshgrid(-Nx/2:Nx/2-1, -Ny/2:Ny/2-1);
sumsq = X.^2 + Y.^2;

C_func = abs(fft2(sumsq<r^2)).^2;

end
