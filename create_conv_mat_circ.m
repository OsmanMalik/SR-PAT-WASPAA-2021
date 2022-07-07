function ConvMat = create_conv_mat_circ(PSF)
%CREATE_CONV_MAT_CIRC Creates matrix corresponding to 2D circular convolution
% 
%   ConvMat = CREATE_CONV_MAT_CIRC(PSF) returns the 2D matrix corresponding
%   to 2D circular convolution conv2(D, PSD), where D is some density. Note
%   that the size of the matrix describing the PSF must be the same as the
%   size of the matrix describing the density D.

[M, N] = size(PSF);

H = cell(M,1);
for m = 1:M
    H{m} = toeplitz(PSF(m,:), [PSF(m,1) fliplr(PSF(m, 2:end))]);
end

ConvMat = zeros(M*N);
for m = 1:M
    for c = 1:M
        ConvMat(1+(m-1)*M:m*M, 1+(c-1)*N:c*N) = H{1+mod(m-c,M)};
    end
end

% Fix up ConvMat to avoid having to do fftshift on matrix-vector product
I = 1:size(ConvMat,1);
Im = reshape(I, M, N);
Imf = fftshift(Im);
If = Imf(:);
ConvMat(If, :) = ConvMat(I, :);

end
