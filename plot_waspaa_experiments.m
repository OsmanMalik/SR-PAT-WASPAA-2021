% This script creates Figures 4 and 5 in the paper. Setting plot_circ to
% true creates the plot for the circular array geometry (Fig. 5) and
% setting it to false creates the plot for the square array 
% geometry (Fig. 4). The appropriate scripts need to be run before running
% this script in order to ensure that all the relevant mat files are
% present. See the README for further details.

f1 = figure;
plot_circ = true;

if false
    load N_gen_101_N_rec_81_r_3_circ
    load('kWave_recon_circular')
else
    load N_gen_101_N_rec_81_r_3_surf
    load('kWave_recon_square')
end

t = tiledlayout(2,3,'Padding', 'none', 'TileSpacing','Compact');

font_size = 13;

nexttile(1)
imagesc(Rho)
title('(a) Original', 'fontsize', font_size)
xticks([])
yticks([])

nexttile(2)
imagesc(p0_recon)
title('(b) k-Wave recon.', 'fontsize', font_size)
xticks([])
yticks([])

nexttile(4)
imagesc(Out)
title('(d) 2nd order recon.', 'fontsize', font_size)
xticks([])
yticks([])


if plot_circ
    load N_gen_101_N_rec_81_r_6_circ
else
    load N_gen_101_N_rec_81_r_6_surf
end

nexttile(5)
imagesc(Out)
title('(e) 2nd order recon.', 'fontsize', font_size)
xticks([])
yticks([])

if plot_circ
    load N_gen_101_N_rec_81_r_12_circ
else
    load N_gen_101_N_rec_81_r_12_surf
end

nexttile(3)
imagesc(Out_unif)
title('(c) 1st order recon.', 'fontsize', font_size)
xticks([])
yticks([])

nexttile(6)
imagesc(Out)
title('(f) 2nd order recon.', 'fontsize', font_size)
xticks([])
yticks([])

axesHandles = findobj(get(t,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
colormap hot
