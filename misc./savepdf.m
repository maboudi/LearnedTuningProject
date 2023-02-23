
function savepdf(gcf, filename, format)

% '-dpdf' for pdf and '-dpng' for png ...

set(gcf,'units','centimeters');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

print(gcf, filename, format, '-r0')

end