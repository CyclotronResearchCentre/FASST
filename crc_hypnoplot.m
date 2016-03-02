function crc_hypnoplot(hand,handles,windowsize,score,scorer)

% FORMAT crc_hypnoplot(hand,handles,windowsize,score,scorer)
% Plotting hypnogram for scored sleep data set.
%__________________________________________________________________
% Copyright (C) 2013 Cyclotron Research Centre

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liège, Belgium
% $Id$

cla(hand)
delete(get(hand,'Children'));
def = crc_get_defaults('score');

if iscell(handles.Dmeg)
    D = handles.Dmeg{1};
else
    D = handles.Dmeg;
end

if nargin < 5 && isfield(handles,'score')
    score = handles.score{1,handles.currentscore};
elseif nargin<5 && ~isfield(handles,'score')
    score = D.CRC.score{1,handles.currentscore};
    handles.score = D.CRC.score;
end

if nargin < 3
    if isfield(handles,'currentscore')
        windowsize = handles.score{2,handles.currentscore};
    else
        windowsize = def.winsize;
    end
end

vect = -.5:1/windowsize:(.5-1/windowsize);

if isfield(handles,'scoring') && handles.scoring
    Zero = find(score==0);
    mxZero = ones(length(vect),1)*Zero;
    mxvect = vect'*ones(1,length(Zero));
    Zero = reshape(mxZero+mxvect, 1, numel(mxvect));
    
    One = find(score==1);
    mxOne = ones(length(vect),1)*One;
    mxvect = vect'*ones(1,length(One));
    One = reshape(mxOne+mxvect, 1, numel(mxvect));
    
    Two = find(score==2);
    mxTwo = ones(length(vect),1)*Two;
    mxvect = vect'*ones(1,length(Two));
    Two = reshape(mxTwo+mxvect, 1, numel(mxvect));
    
    Three = find(score==3);
    mxThree = ones(length(vect),1)*Three;
    mxvect = vect'*ones(1,length(Three));
    Three = reshape(mxThree+mxvect, 1, numel(mxvect));
    
    Four = find(score==4);
    mxFour = ones(length(vect),1)*Four;
    mxvect = vect'*ones(1,length(Four));
    Four = reshape(mxFour+mxvect, 1, numel(mxvect));
    
    Five = find(score==5);
    mxFive = ones(length(vect),1)*Five;
    mxvect = vect'*ones(1,length(Five));
    Five = reshape(mxFive+mxvect, 1, numel(mxvect));
    
    Six = find(score==6);
    mxSix = ones(length(vect),1)*Six;
    mxvect = vect'*ones(1,length(Six));
    Six = reshape(mxSix+mxvect, 1, numel(mxvect));
    
    Seven = find(score==7);
    mxSeven = ones(length(vect),1)*Seven;
    mxvect = vect'*ones(1,length(Seven));
    Seven = reshape(mxSeven+mxvect, 1, numel(mxvect));
end

labl = def.stnames_S(end:-1:1);
totime = ceil(nsamples(D)/fsample(D));

% Display Wake state (st0 => 7)
if nargin < 2
    d = axes('YTick',[0 1 2 3 4 5 6 7],'YTickLabel',labl,...
        'XTick',sort((totime):-(totime/10):0 ),...
        'Fontsize',8);
else
    set(hand,'YTick',[0 1 2 3 4 5 6 7],'YTickLabel',labl,...
        'XTick',sort((totime):-(totime/10):0 ),...
        'Fontsize',8);
    d = hand;
end

hold on

if isfield(handles,'scoring') && handles.scoring
    % Display awake state
    if size(Zero,2)>0
        bar((Zero*windowsize-windowsize/2),ones(1,length(Zero))*7,1, ...
            'LineStyle','none','FaceColor',[0.2 0.75 0.6])
    end
    % Display Sleep stage 1
    if size(One,2)>0
        bar((One*windowsize-windowsize/2 ),ones(1,length(One))*6,1, ...
            'LineStyle','none','FaceColor',[0 0.8 1])
    end
    % Display Sleep stage 2
    if size(Two,2)>0
        bar((Two*windowsize-windowsize/2),ones(1,length(Two))*5,1, ...
            'LineStyle','none','FaceColor',[0.1 0.5 0.9])
    end
    % Display Sleep stage 3
    if size(Three,2)>0
        bar(Three*windowsize-windowsize/2,ones(1,length(Three))*4,1, ...
            'LineStyle','none','FaceColor',[0.1 0.2 0.8])
    end
    % Display Sleep stage 4
    if size(Four,2)>0
        bar(Four*windowsize-windowsize/2,ones(1,length(Four))*3,1, ...
            'LineStyle','none','FaceColor',[0.1 0.15 0.5])
    end
    % Display REM state
    if size(Five,2)>0
        bar(Five*windowsize-windowsize/2,ones(1,length(Five))*2,1, ...
            'LineStyle','none','FaceColor',[0.5 0.5 0.9])
    end
    % Display MT
    if size(Six,2)>0
        bar(Six*windowsize-windowsize/2,ones(1,length(Six))*1,1, ...
            'LineStyle','none','FaceColor',[0.9 0.4 0.4])
    end
    % Display Unscorable
    if size(Seven,2)>0
        bar(Seven*windowsize-windowsize/2,ones(1,length(Seven))*0.35,1, ...
            'LineStyle','none','FaceColor',[0.9 0.6 0.3])
    end
    
    % display the artefacts, arousals and events of interest as small bars on
    % the top of the hypnogram
end

art  = handles.score{5,handles.currentscore};
arou = handles.score{6,handles.currentscore};
eoi  = handles.score{7,handles.currentscore};

if ~isempty(art)
    art = art(:,1);
    plot(art, ones(length(art),1)*8, '+', 'MarkerSize',8, ...
        'MarkerFaceColor',[0.4 0.4 0.4],...
        'MarkerEdgeColor',[0.4 0.4 0.4],'tag','undart')
end
if ~isempty(arou)
    arou = arou(:,1);
    plot(arou, ones(length(arou),1)*8, '+', 'MarkerSize',8, ...
        'MarkerFaceColor',[1 0 0],...
        'MarkerEdgeColor',[1 0 0],'tag','arou')
end
if ~isempty(eoi)
    eoi = eoi(:,1);
    plot(eoi, ones(length(eoi),1)*8, '+', 'MarkerSize',8, ...
        'MarkerFaceColor',[0.75 0.2 0.2],...
        'MarkerEdgeColor',[0.75 0.2 0.2],'tag','eoi')
end

if isfield(handles,'type')
    evman = handles.type;
else
    evman = {};
end
evtot = events(D);
if ~isempty(evtot)
    type_evtot = cellstr(char(evtot(:).type));
    nevm = size(evman,1);
    
    for evm = 1 : nevm
        evt = find(strcmpi(type_evtot,evman(evm,1)));
        evdis = cell2mat({evtot(evt).time});
        col = char(evman(evm,2));
        indcol = col(1);
        if ~isempty(evdis)
            plot(evdis, ones(length(evdis),1)*8, 'V', ...
                'MarkerSize',5,'MarkerFaceColor',indcol,...
                'MarkerEdgeColor',indcol)
            for vert = 1 : length(evdis)
                plot(ones(1,2)*evdis(vert),[0 8],'-.','Color',indcol)
            end
        end
    end
end

if isfield(D,'info')
    if isfield(D.info,'hour')
        xtick = mod(get(d,'XTick') + handles.offset,24*60^2);
    else
        xtick = get(d,'XTick');
    end
elseif isfield(handles,'offset')
    xtick = get(d,'XTick');
    xtick = mod(xtick + handles.offset,24*60^2);
else
    xtick = get(d,'XTick');
end
[dumb, times] = crc_time_converts(xtick);
set(d,'XTickLabel',times)

ylim([0 8])
xlim([0 totime])
grid on

if isfield(handles,'slider1')
    hold on
    slidval = get(handles.slider1,'Value');
    pos = min(slidval+handles.winsize/2,nsamples(D)/fsample(D));
    handles.cursor = plot(pos,8,'v','Color',[0.2 0.2 0.2], ...
        'LineWidth',2.5,'tag','cursor');
    hold off
end

if isfield(handles,'scoring') && handles.scoring
    hold on
    if nargin<5
        scorer = handles.score{2,handles.currentscore};
    end
    text(totime-(9*totime/100),7.5,scorer, ...
        'Color',[0.3 0.3 0.3],'Fontweight','demi')
end

return
