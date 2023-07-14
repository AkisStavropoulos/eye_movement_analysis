%% check gaussian and butter filtering
%% define filter (done)
sig = prs.filtwidth; %filter width
sz = prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

%% filter position and velocity channels (done)
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},{'tsi','mrk','yle','yre','zle','zre'}))
        ch1.(chnames{i}) = conv(ch.(chnames{i})(1:MAX_LENGTH),h,'same');
%          ch1.(chnames{i}) = ch1.(chnames{i})(sz/2+1:end);
    end
end
%% Jean's filter
[b,a] = butter(2,30/(SR/2),'low');
for i=1:length(chnames)
        if ~any(strcmp(chnames{i},{'tsi','mrk','yle','yre','zle','zre'}))
            ch2.(chnames{i}) = filtfilt(b,a,ch.(chnames{i}));
        end
end
%% Plot all FILTERED channels
chnames = fieldnames(ch);
for i = 1:length(ch_title)
    if ~any(strcmp(chnames{i},{'tsi','mrk','yle','yre','zle','zre'}))
        d = length(ch.(chnames{i})) - length(ts);
        d2 = length(ch2.(chnames{i})) - length(ts);
        figure;plot(ts,ch.(chnames{i})(1:end-d));title(chnames{i});xlabel('time (s)');vline([t.beg]);
        hold on;plot(ts,ch1.(chnames{i}),'g','LineWidth',2);plot(ts,ch2.(chnames{i})(1:end-d2),'r','LineWidth',2);hold off;
    end
end
