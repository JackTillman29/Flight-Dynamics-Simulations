function varargout = PrepDataCallout(varargin)
% There are 2 uses:

%    Method 1.)
%       Get a structure containing default values for callout visual
%       properties to be used in method 2):
%              cp = PrepDataCallout;
%                  ** note it is highly recommended you use the name "cp"
%                     for this output because this is what is assumed by
%                     callbacks to update the visuals using the context
%                     menu approach.
%    Method 2.)
%       Use structure (likely modified) obtained in method 1 to update
%       callout visual properties in figure defined by passed in figure
%       handle:
%              PrepDataCallout(gcf,cp)
%              PrepDataCallout(1,cp)
%              PrepDataCallout(hfigure,cp)

if(nargin == 0)
    callout_props.axes_border_linewidth = 1;
    callout_props.axes_border_color     = 'k';
    callout_props.shadow_color          = 'k';
    callout_props.shadow_alpha          = 0.5;
    callout_props.shadow_edge_color     = 'none';
    callout_props.shadow_edge_linewidth = 1;
    
    varargout{1} = callout_props;
    return
elseif(nargin == 2)
    
    % update visual properties of the existing figure's callouts with
    % "callout_props"
    hfig = varargin{1};
    callout_props = varargin{2};
    
    % update visuals
    cb_make_data_callout__update_properties(hfig);
else
    error('wrong number of input arguments (0 or 2)')
end






end