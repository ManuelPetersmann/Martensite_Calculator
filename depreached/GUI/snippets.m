to define three edit boxes with white text background and left text alignment:
uicontrol('Units', 'normalized', 'Position', [0.1,0.9,0.7,0.05], ...
'HorizontalAlignment', 'left', 'Style', 'edit', 'BackgroundColor', [1,1,1]);
uicontrol('Units', 'normalized', 'Position', [0.1,0.8,0.7,0.05], ...
'HorizontalAlignment', 'left', 'Style', 'edit', 'BackgroundColor', [1,1,1]);
uicontrol('Units', 'normalized', 'Position', [0.1,0.7,0.7,0.05], ...
'HorizontalAlignment', 'left', 'Style', 'edit', 'BackgroundColor', [1,1,1]);

A vectorized call to set reduces this to
h(1) = uicontrol('Units', 'normalized', 'Position', [0.1,0.9,0.7,0.05]);
h(2) = uicontrol('Units', 'normalized', 'Position', [0.1,0.8,0.7,0.05]);
h(3) = uicontrol('Units', 'normalized', 'Position', [0.1,0.7,0.7,0.05]);
set(h, 'HorizontalAlignment', 'left', 'Style', 'edit','BackgroundColor', [1,1,1]);

