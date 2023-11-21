function vertical_errorbars(x, width, mu, sem, color, linewidth)

line([x x], [mu - sem mu + sem], 'Color', color, 'linewidth', linewidth); hold on;
line([x-width x+width], [mu - sem mu - sem], 'Color', color, 'linewidth', linewidth); hold on;
line([x-width x+width], [mu + sem mu + sem], 'Color', color, 'linewidth', linewidth); hold on;






