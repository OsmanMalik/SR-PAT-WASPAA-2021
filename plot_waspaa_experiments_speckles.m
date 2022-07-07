% This script is used to generate the plot in Figure 3 of the paper. The
% appropriate scripts need to be run first to ensure that the relevant mat
% files are available. See the README file for details.

t = tiledlayout(1,3,'Padding', 'none', 'TileSpacing','Compact');
font_size = 13;

load N_gen_101_N_rec_81_r_3_circ

nexttile(1)
imagesc(spec)
title('(a) Size: 16.0 \mum', 'fontsize', font_size)
xticks([])
yticks([])

load N_gen_101_N_rec_81_r_6_circ

nexttile(2)
imagesc(spec)
title('(b) Size: 7.8 \mum', 'fontsize', font_size)
xticks([])
yticks([])

load N_gen_101_N_rec_81_r_12_circ

nexttile(3)
imagesc(spec)
title('(c) Size: 3.9 \mum', 'fontsize', font_size)
xticks([])
yticks([])

axesHandles = findobj(get(t,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
colormap hot
