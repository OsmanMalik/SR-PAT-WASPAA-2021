function C = Compute_Covariance_Matrix(Nx, Ny, r)
%COMPUTE_COVARIANCE_MATRIX  Compute covariance matrix.
%
%   C = COMPUTE_COVARIANCE_MATRIX(Nx, Ny, r) computes the covariance matrix
%   corresponding to the covariance function in Equation (14) of [1].
%
% REFERENCES:
%   [1] J. C. Dainty. The statistics of speckle patterns. Progress in
%   Optics XIV, 1976.

C_func = Compute_Covariance_Function(Nx, Ny, r);
C = create_conv_mat_circ(fftshift(C_func));

end
