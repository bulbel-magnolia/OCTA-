function ini=iniset(varargin)


[ppath,~,~] = fileparts(mfilename('fullpath'));

if nargin==0
    inilist0=dir([pwd,'/config*.ini']);
    inilist1=dir([pwd,'/../config*.ini']);
    inilist2=dir([pwd,'/../../config*.ini']);
    iniFallback=dir([ppath,'/config*.ini']);
    if ~isempty(inilist0)
        ininame=inilist0(1).name;
        if nargout>0,ini=ini2struct([pwd,'/',ininame]);else, open([pwd,'/',ininame]);end
        fprintf('ini found in current directory\n')
    elseif ~isempty(inilist1)
        ininame=inilist1(1).name;
        if nargout>0,ini=ini2struct([pwd,'/../',ininame]);else, open([pwd,'/../',ininame]);end
        fprintf('ini found in //.. directory\n')
    elseif ~isempty(inilist2)
        ininame=inilist2(1).name;
        if nargout>0,ini=ini2struct([pwd,'/../../',ininame]);else, open([pwd,'/../../',ininame]);end
        fprintf('ini found in ..//.. directory\n')
    elseif ~isempty(iniFallback)
        ininame=iniFallback(1).name;
        if nargout>0,ini=ini2struct([ppath,'/',ininame]);else, open([ppath,'/',ininame]);end
        fprintf('ini found in program directory\n')
    else
        error('ini not found\n')
    end
    else
    if 2==exist(varargin{1},'file')
        if nargout>0,ini=ini2struct(varargin{1});else, open(varargin{1});end
        fprintf('ini loaded\n')
    else
        error('ini not found\n')
    end
end
