function retMat = D2D_Forward_Implicit(P0, H0, E, no_trans, padding)
%D2D_Forward_Implicit Implicit discrete-to-discrete forward model.
%
% NOTES:
%   This function computes the D2D forward operator implicitly and
%   therefore avoids having to store the system matrix as a dense matrix.
%
% INPUTS:
%   P0          - The initial pressure distribution expressed as a vector
%                 or a matrix with each initial pressure distribution along
%                 the columns.
%   H0          - The sparse H0 system matrix.
%   E           - The sampled EIR convolved with any non-delta shaped
%                 illumination.
%   no_trans    - The number of transducers used.
%   padding     - The amount of padding inserted after each transducer
%                 block in H0.
%
% OUTPUTS:
%   retMat      - The vector or matrix with observations in each column
%                 corresponding to the initial pressure distribution in the
%                 same columns in P0.

FE = repmat(fft(E, size(H0,1) + no_trans*padding),1,size(P0,2));
Y = ifft(fft(pad(1, H0*P0, no_trans, padding)) .* FE);
retMat = unpad(2, Y, no_trans, padding);

end
