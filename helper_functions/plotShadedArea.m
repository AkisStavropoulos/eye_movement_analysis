function plotShadedArea(x_start, x_end, ylimz,colr, transparency)


if isempty(ylimz)
    ylimz=get(gca,'ylim');
end

fill([x_start x_end x_end x_start],[ylimz(1) ylimz(1) ylimz(2) ylimz(2)],colr,'FaceAlpha',transparency,'edgecolor','none')
