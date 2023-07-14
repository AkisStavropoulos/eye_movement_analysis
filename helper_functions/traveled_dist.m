function travel_d = traveled_dist(trials)

%% find the distance traveled by the subject for each trial
travel_d = [];
for i = 1:length(trials)
        dist = sqrt((trials(i).continuous.xmp(end) - trials(i).continuous.xmp(1)).^2 + (trials(i).continuous.ymp(end) - trials(i).continuous.ymp(1)).^2);
    travel_d = [travel_d dist];
end