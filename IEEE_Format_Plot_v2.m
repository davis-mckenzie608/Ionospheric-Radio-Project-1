function IEEE_Format_Plot_v2(varargin)
% Desired changes:
% make 'hwr' input that is height to width ratio
% this can just go into _v0.
% write out a good way to incorporate this code.
% make calls for geo plot stuff
% some defaults
scale = 4;
% get current figure whr (width height ratio = height/width)
fig = gcf;
width = 8.8 * scale;
width_height_ratio = fig.Position(3)/fig.Position(4);
height = width / width_height_ratio;
% find where types is a certain class
% scale factor option
for ii = 1:nargin
    if strcmp('scale', varargin{ii})
        scale_index = ii+1;
    end
    % figure handle option
    if strcmp('figure', varargin{ii})
        figure_index = ii+1;
    end
    % width height ratio option
    if strcmp('whr', varargin{ii})
        whr_index = ii+1;
    end
    % axis handle option
    if strcmp('axis', varargin{ii})
        axis_index = ii+1;
    end
    % curve handle option
    if strcmp('curve', varargin{ii})
        curve_index = ii+1;
    end
    % legend handle option
    if strcmp('legend', varargin{ii})
        legend_index = ii+1;
    end
end
%edit scale per the input if it exists
if exist('scale_index', 'var')
    scale = varargin{scale_index};
end
if exist('whr_index', 'var')
    width_height_ratio = varargin{whr_index};
    width = 8.8 * scale;
    height = width / width_height_ratio;
end
%edit figure per the input if it exists
if exist('figure_index', 'var')
    fig = varargin{figure_index};
    set(fig, 'Units','centimeters');
    set(fig, 'Position', [1 2 width height]);
    %edit legend per the input if it exists
end
%edit axis per the input if it exists
if exist('axis_index', 'var')
    ax = varargin{axis_index};
    if iscell(ax)
        for jj = 1:length(ax)
            set(ax{jj},'GridLineStyle','--');
            set(ax{jj},'YGrid','on');
            set(ax{jj},'XGrid','on');
            set(ax{jj},'FontName','Times New Roman');
            set(ax{jj},'FontSize',8*scale);
            set(ax{jj},'FontName','Times New Roman');
            set(ax{jj},'FontSize',8*scale);
            set(ax{jj},'Units','normalized');
            % set(ax{jj},'color','none');
        end
    else
        set(ax,'GridLineStyle','--');
        set(ax,'YGrid','on');
        set(ax,'XGrid','on');
        set(ax,'FontName','Times New Roman');
        set(ax,'FontSize',8*scale);
        set(ax,'FontName','Times New Roman');
        set(ax,'FontSize',8*scale);
        set(ax,'Units','normalized');
        % set(ax,'color','none');
    end
end
%edit curve per the input if it exists
if exist('curve_index', 'var')
    curve = varargin{curve_index};
    if iscell(curve)
        for jj = 1:length(curve)
            set(curve{jj},'LineWidth',0.5*scale);
            % set(curve{jj},'MarkerSize',4*scale);
        end
    else
        set(curve,'LineWidth',0.5*scale);
        set(curve,'MarkerSize',2*scale);
    end
end
%edit legend per the input if it exists
if exist('legend_index', 'var')
    lgd = varargin{legend_index};
    set(lgd,'Visible','on');
    set(lgd,'FontName','Times New Roman');
    set(lgd,'FontSize',6*scale);
    set(lgd,'Location','best');
end
end