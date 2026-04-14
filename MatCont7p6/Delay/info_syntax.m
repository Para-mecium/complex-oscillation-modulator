function varargout = info_syntax(varargin)
% INFO_SYNTAX MATLAB code for info_syntax.fig
%      INFO_SYNTAX, by itself, creates a new INFO_SYNTAX or raises the existing
%      singleton*.
%
%      H = INFO_SYNTAX returns the handle to a new INFO_SYNTAX or the handle to
%      the existing singleton*.
%
%      INFO_SYNTAX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INFO_SYNTAX.M with the given input arguments.
%
%      INFO_SYNTAX('Property','Value',...) creates a new INFO_SYNTAX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before info_syntax_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to info_syntax_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help info_syntax

% Last Modified by GUIDE v2.5 15-Oct-2023 15:45:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @info_syntax_OpeningFcn, ...
                   'gui_OutputFcn',  @info_syntax_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before info_syntax is made visible.
function info_syntax_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to info_syntax (see VARARGIN)

% Choose default command line output for info_syntax
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes info_syntax wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.text7,'String',textwrap(handles.text7,{
    'Information on the input syntax and the discretization of delay equations.'}));
set(handles.text6,'String',textwrap(handles.text6,{
    ['The importer can handle delay differential equations and renewal equations ', ...
    'with discrete and distributed delays.'], ...
    '', ...
    ['- Delay differential equations, like ordinary equations, are indicated ', ...
    'by using the prime ('') on the left-hand side, ', ...
    'which is not relevant for renewal equations.'], ...
    ['- The time dependency of variables is specified with square brackets: ', ...
    'e.g. x(t-2) is input as x[t-2]. ', ...
    'The time dependency can optionally be specified also for current-time terms, ', ...
    'e.g. x[t], or x''[t] in the left-hand side.'], ...
    ['- Distributed delays are specified with the DE_int function: ', ...
    ' e.g. the integral from a to b of x[t+theta] in dtheta is input as ', ...
    'DE_int(@(theta)x[t+theta],a,b). '], ...
    ['- Delayed time dependencies are specified in the form [t-expr], ', ...
    'where expr may depend on the parameters (and on the integration variable ', ...
    'in the distributed case) but not on t or on the coordinates. ', ...
    'The evaluation of expr should result in a positive numeric value. ', ...
    'You can also input [t+expr2], with expr2 negative.'], ...
    '', ...
    ['The discretization of the delay equation results in a system of ordinary differential equations, ', ...
    'obtained via pseudospectral collocation on Chebyshev extremal nodes on the delay interval. ', ...
    'Distributed delays are approximated with the Clenshaw-Curtis quadrature formula.']
    }));


% --- Outputs from this function are returned to the command line.
function varargout = info_syntax_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
