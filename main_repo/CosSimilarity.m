function cos_similarity = CosSimilarity(X,Y)

X = squeeze(X); Y = squeeze(Y);
xy = (sum(X.*Y));
xx = (sqrt(sum(X.*X)));
yy = (sqrt(sum(Y.*Y)));
cos_similarity = nanmean(xy./(xx.*yy));