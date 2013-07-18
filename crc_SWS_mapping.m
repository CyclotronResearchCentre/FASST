function [hb,htitl]=crc_SWS_mapping (D,value,choice,zebr,dim,fig,ax)

% FUNCTION crc_SWS_mapping
%
% maps slow waves of deep human sleep thanks to two views :
%
%  - electrical signal averaged on 4 ROIs plotted against time in the
%[-500ms, +1500ms] interval around the concerned (value) SW.
%
%  - 2D flattened map of the scalp showing either a hot colormap of delays
%  (choice 1) or a jet movie showing potentials on the scalp in the same
%  time interval (choice 2).
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Written by J. Schrouff & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if nargin<1
    D = crc_eeg_load;
end
if nargin<2, 
    value = spm_input ('Number of the wave to show:',1,'r',1,[]);
    close gcf 
end

if nargin<3, 
    choice = spm_input('Map of:',+1,'b','delays|potentials',[1;2],1);
    close gcf
end
if nargin<5
    dim=spm_input('Electrodes positions:',+1,'b','file|auto',[0,1],0);
    close gcf
end
if nargin<4 && dim==0
    zebr=spm_select(Inf, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
    close gcf
end

if dim==0
    el_set_sphuni = crc_elecpos_sph(D,zebr);
else
    el_set_sphuni=[];
end
origin_count= pm_origin_count(D);
[hb,htitl]=pm_map(D.CRC.SW, D, origin_count, choice, value, dim, el_set_sphuni,fig,ax);

%--------------------------------------------------------------------------
%-----------------    SUBFUNCTION TO MAP THE DELAYS   ---------------------
%--------------------------------------------------------------------------

function [hb,htitl]=pm_map(SW, data, origin_count, choice, value, dim, el_set_sphuni,fig,ax)


% choice 1 shows delays on a 'hot' colour map of the scalp. the first
% electrode is a pink star. Gray represents electrodes that don't detect the
% considered wave.
% choice 2 draws the data on the scalp in a jet colormap of potentials and
% traces the main trajectory of the wave.

crcdef=crc_get_defaults('swsd');
bad=badchannels(data);
notbad=setdiff(1:nchannels(data),bad);
eeg_chan=find(strcmpi(crcdef.type,chantype(data)));
data_channels_eeg=notbad(eeg_chan);
all_names = data.chanlabels;
channels_eeg_labels = all_names(data_channels_eeg)';
el_label = channels_eeg_labels;

if nargin<8
    fig=0;
end
if nargin<9
    ax=0;
end

dr=3;

switch choice
    case 1

        ed=ones(size(el_label,1),1);
        sv2= zeros(size(el_label,1),1);
        for ielec = 1:size(SW.SW(value).electrodes,2)
            for ichannel = 1:size(el_label,1)
                if strcmpi(char(SW.SW(value).electrodes(ielec)), ...
                        char(el_label{ichannel}))
                    ed(ichannel,1)=SW.SW(value).delays(ielec)-SW.SW(value).delays(1)+dr;
                    sv2(ichannel,1)= 1;
                end
            end
        end
        if dim==0
            %uses interpolation on the unit sphere
            SV = crc_CreateSV(ed,[],0,...
                el_set_sphuni',1,7*pi/12,40,el_label);
            SV2= crc_CreateSV(sv2,[],0,...
                el_set_sphuni',1,7*pi/12,40,el_label);
            mask = abs(SV2.interpolated)>0.85;
            SV.interpolated=SV.interpolated .* mask;
            if SW.SW(value).code<400
                typsw='SW';
            else
                typsw='delta';
            end                
            label_fig=[char(SW.SW(value).electrodes(1)), ' , type ', typsw];
            [SV]=crc_DrawSV(SV,fig,0,{label_fig},0,0,10,ax);
            hb=SV.hColourBar;
            htitl=SV.hTitl;
            colormap(hot)
            caxis([0 max(ed)]);
            axis off
            cc=colormap;
            cc(1:end,:)=cc(end:-1:1,:);
            cc(1,:)=[0.5 0.5 0.5];
            colormap(cc)
            axis off
        else
            % uses interpolation on a unit disk
            load('CRC_electrodes.mat')
            try
                pos_eeg_chan=(coor2D(data,data_channels_eeg))';
            catch
                pos_eeg_chan=zeros(numel(eeg_chan,2));
                dumb2=zeros(numel(eeg_chan,2));
                for ich=1:numel(eeg_chan)
                    iselec=strcmpi(chanlabels(data,eeg_chan(ich)),names);
                    dumb2(ich)=find(iselec);
                    pos_eeg_chan(ich,:)=pos(find(iselec),:);
                end
            end
            [zv,hb]=interp2D(pos_eeg_chan,ed,sv2,1);
            xlim([0,1])
            ylim([0,1])
            if SW.SW(value).code<400
                typsw='SW';
            else
                typsw='delta';
            end                
            label_fig=[char(SW.SW(value).electrodes(1)), ' , type ', typsw];
            htitl=title(label_fig);
        end



    case 2
        if dim==0
            %uses spline interpolation on the unit sphere
            SV = crc_CreateSV(data(data_channels_eeg, ...
                              SW.SW(value).negmax_tp-50:SW.SW(value).negmax_tp+250), ...
                              20:220,0, el_set_sphuni',1,pi/2,10,el_label);
            if SW.SW(value).code<400
                typsw='SW';
            else
                typsw='delta';
            end                
            label_fig=[char(SW.SW(value).electrodes(1)), ' , type ', typsw];
            [SV]=crc_DrawSV(SV,fig,0,{label_fig},1,0,10,ax);
            hb=SV.hColourBar;
            htitl=SV.hTitl;
            axis off
        else
            disp('Option not available yet, use delay maps')
        end

end

%--------------------------------------------------------------------------
%-----------  SUBFUNCTION TO INITIALIZE COUNT ON CHANNELS   ---------------
%--------------------------------------------------------------------------
function [origin_count] = pm_origin_count(data)

% data_channels_eeg=find(strcmpi('eeg',chantype(data)));
all_names = data.chanlabels;
data_channels_eeg=find(strcmpi('eeg',chantype(data)));
origin_count=all_names(data_channels_eeg)';
for i=1:size(data_channels_eeg,2)
    origin_count(i,2)={0};
    origin_count(i,3)={0};
end
origin_count(:,2)={0};
origin_count(:,3)={0};

%--------------------------------------------------------------------------
%-----------  SUBFUNCTION TO INTERPOLATE in 2D   --------------------------
%--------------------------------------------------------------------------

function [zv,hb]=interp2D(el_set_2d,delval,maskval,dra,dst)

% Script for 2D interpolation using matlab built in functions
% Written by Y. Leclercq & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium

if nargin<5
    dst=0.02;
end

% Create space
X=el_set_2d(:,1);
Y=el_set_2d(:,2);
R1 = max(X);
R2=max(Y);
if R1>=R2
    [XI,YI] = meshgrid(-R1:dst:R1);
else
    [XI,YI] = meshgrid(-R2:dst:R2);
end

V=delval;
M=maskval;

% Interpolation & display
[xl,yl,zl] = griddata(X,Y,V,XI,YI,'linear');
[xv,yv,zv] = griddata(X,Y,V,XI,YI,'v4');
[xma,yma,zma]=griddata(X,Y,M,XI,YI,'cubic');
zm = zl; zm(~isnan(zl)) = 1;
zv=zv.*zm;
zmask=zma;
zmask(zma<0.8)=nan;
zv=zv.*zmask;
if dra==1
    pcolor(xv,yv,zv)
    axis image
    hold on
    plot(X,Y,'x')
    colormap(jet)
    spm_figure('Colormap','Invert')
    cc=colormap;
    cc(1,:)=[0 0 0];
    colormap(cc)
    hb=colorbar;
    shading interp
else
end
return

%_________________________________________________________________________________
function [ xf, yf ] = CartToFlat( x, y, z )
% Convert Cartesian coordinates on surface of unit sphere to 2D Cartesian
% Flat Map coordinates


rh = sqrt( x.*x + y.*y );

rh( rh == 0 ) = 1e-20;

costheta = x ./ rh;
sintheta = y ./ rh;

% Modification by chrisp@fil, 02/02/24,
% to take into account electrodes with phi<-90 or phi>90
%---------------------------------------------------------------------
% lz = find(z<0) ;
phi   = atan2( sqrt( x.*x + y.*y ) , z );
% phi(lz) = pi + atan(  sqrt( x(lz).*x(lz) + y(lz).*y(lz) )  ./  z(lz)  );

xf = costheta .* phi;
yf = sintheta .* phi;

%--------------------------------------------------------------------------
%---------  SUBFUNCTION TO DISPLAY DATA AND POINT WAVES  ------------------
%--------------------------------------------------------------------------
% 
% function fig = pm_display(dispdata, position, displstep, t, sss, ...
%                             negmax, yLim, uiscroll, disptext, ...
%                             delays, roi_selec, namroi, crcdef)
% 
% if ~nargin
%     position = [10 10 800 800];
%     displstep = 75;
%     t = 1:size(dispdata,2);
% end
% 
% %figure
% fig=spm_figure('Create','1','Display of detected Slow Waves','on');
% subplot(2,1,1)
% axis off
% hold;
% WS = spm('WinScale');
% position=position.*WS;
% 
% for jj = 1:size(dispdata,1)
%     if ~isempty(disptext)
%         d=(position(4))/(position(4)*10);
%         axes('Position',[0.1,0.85-(jj-1)*d,0.8,d]);
%         plot(t,dispdata(jj,:))
%         plot(delays(jj),min(dispdata(jj,:)),'o')
%         axis ([t(1) t(end) yLim(1) yLim(end)]);
%         str = disptext(jj);
%         text(t(end)+.01,-(jj-1)*displstep,str)
%         grid on
%     else
%         d=(position(4))/(position(4)*10);
%         axes('Position',[0.1,0.85-(jj-1)*d,0.8,d]);
%         plot(t,dispdata(jj,:),'LineWidth',2)
%         axis ([t(1) t(end) yLim(1) yLim(end)]);
%         if roi_selec==1
%             if jj==1
%                 c='frontal ROI';
%             elseif jj==2
%                 c='left ROI';
%             elseif jj==3
%                 c='right ROI';
%             else
%                 c='rear ROI';
%             end
%         elseif roi_selec==2
%             c=namroi{jj};
%         end
%         ylabel(c)
%         grid on
%     end
%     hold on;
%     if ~isempty(sss)
% 
%         hold on;
%         pmin=max(negmax-150,1);
%         pmax=min(negmax+200,length(t));
%         px=[t(pmin)  t(pmin)  t(pmax)  t(pmax)];
%         py=[crcdef.dispscale -crcdef.dispscale -crcdef.dispscale crcdef.dispscale];
%         col=[0.9 0.9 0.1];
%         p=patch(px,py,col,'FaceAlpha', 0.5);
%         hold on;
%         plot(t,0,'k','LineWidth',2)
%     end
% end
% if uiscroll == 1
%     uiscroll([1 size(dispdata,2)],10,[],[],yLim);
% end
