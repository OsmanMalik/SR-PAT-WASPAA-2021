function Z = create_star_density(no_x_px, no_y_px)
%CREATE_STAR_DENSITY Create star shaped density
%   
%   Z = CREATE_STAR_DENSITY(no_x_px, no_y_px) returns a density Z of size
%   no_x_px times no_y_px of the kind used in [1].
%
% REFERENCES:
%   [1] Idier et al. [IEEE Trans. Comput. Imaging 4(1), 2018].

theta = 20; %40 originally
x = linspace(-1, 1, no_x_px);
y = linspace(-1, 1, no_y_px);
[X, Y] = meshgrid(x,y);
Z = 1 + cos(theta*atan(Y./X));
Z(isnan(Z)) = 0;

end
