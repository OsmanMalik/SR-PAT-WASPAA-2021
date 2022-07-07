function SM = Generate_2D_Speckle(Nx, Ny, r)
%GENERATE_2D_SPECKLE    Generates a two-dimensional speckle pattern with
%                       circular frequency support.
%
%   SM = GENERATE_2D_SPECKLE(Nx, Ny, r) returns a speckle pattern SM of
%   size Nx times Ny. The generated speckle pattern has frequency support
%   which is a disk of radius r.

[X, Y] = meshgrid(-Nx/2:Nx/2-1, -Ny/2:Ny/2-1);
sumsq = X.^2 + Y.^2;
SMf = zeros(size(X));
SMf(sumsq<r^2) = exp(1i*(rand(numel(find(sumsq<r^2)),1)*2*pi));
SM = abs(fft2(SMf)).^2;

end
