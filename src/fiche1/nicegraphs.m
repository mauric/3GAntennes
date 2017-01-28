set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , -3:0.1:2, ...
  'LineWidth'   , 1         );
set( gca, ...
    'FontName'   , 'Helvetica',...
    'FontSize'   , 8);

if(exist('hTitle'))
    set( hTitle,...
		'FontName'   , 'AvantGarde', ...
        'FontSize'   , 12          , ...
        'FontWeight' , 'bold'      );
end

if(exist('hXLabel'))
    set( hXLabel,...
        'FontName'   , 'AvantGarde',...
        'FontSize'   , 10          );
end
if(exist('hYLabel'))
    set( hYLabel,...
        'FontName'   , 'AvantGarde',...
        'FontSize'   , 10          );
end

set(gcf, 'PaperPositionMode', 'auto');