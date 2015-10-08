\NeedsTeXFormat{LaTeX2e}[1996/12/01]%
\ProvidesClass{dccpaper}[1996/09/19]
\typeout{*** 1996/09/19 Class for DCC conference papers ***}

%% Options:
%    notesoff: turn off \Notes (but not \Support) appearance at bottom
%              of first page
%    smallcaptions: typeset captions as small
%
\newif\if@notesoff\@notesofffalse
\DeclareOption{notesoff}{\@notesofftrue}
\newif\if@smallcaptions\@smallcaptionsfalse
\DeclareOption{smallcaptions}{\@smallcaptionstrue}
\newif\if@smallabstract\@smallabstractfalse
\DeclareOption{smallabstract}{\@smallabstracttrue}

%  Pass unknown options to parent classd
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

%  Base everything on the parent class 'article'
%  Must be 12pt type
\LoadClass[12pt]{article}

\RequirePackage{epsfig}
\RequirePackage{amsmath}
\RequirePackage{calc}
\RequirePackage{ifthen}

% DCC page margins (new as of DCC1998):
%    top: 1 in
%    bottom: 1 in
%    sides: 1.25 in
%    text: 9in high and 6in wide
\setlength{\topmargin}{0in}
\setlength{\headheight}{0in}
\setlength{\headsep}{0in}
\setlength{\footskip}{0in}
\setlength{\textheight}{9in}
\setlength{\evensidemargin}{0.25in}
\setlength{\oddsidemargin}{0.25in}
\setlength{\textwidth}{6in}

\renewcommand{\floatpagefraction}{0.95}
\pagestyle{empty}

% For small typeset
\def\smallstretch{0.925}

% To adjust pages longer or shorter by one line
\newcommand{\longpage}{\enlargethispage{\baselineskip}}
\newcommand{\shortpage}{\enlargethispage{-\baselineskip}}

% Use \Section and \SubSection
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\centering\normalfont\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\itshape}}
\def\Section#1{\vspace*{-0.125in}\section*{#1}}
\def\SubSection#1{\vspace*{-0.125in}\subsection*{#1}}

% New title
\newif\if@SupportDefined\@SupportDefinedfalse
\def\Support#1{\def\@Support{#1}\@SupportDefinedtrue}
\newif\if@NotesDefined\@NotesDefinedfalse
\def\Notes#1{\def\@Notes{#1}\@NotesDefinedtrue}
\renewcommand{\maketitle}{\par
  \begingroup
%     \spacingset{1}
    \renewcommand{\thefootnote}{}%
    \def\@makefnmark{}%
    \setcounter{footnote}{1}
    \if@twocolumn
        \twocolumn[\@maketitle]%
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{empty}\@thanks
    \if@SupportDefined
      \if@notesoff\footnote{\@Support}\else\if@NotesDefined\footnote{\@Support\\\@Notes}\else\footnote{\@Support}\fi\fi
    \else
      \if@notesoff\else\if@NotesDefined\footnote{\@Notes}\else\fi\fi
    \fi
  \endgroup
  \setcounter{footnote}{0}%
  \let\thanks\relax
  \let\maketitle\relax\let\@maketitle\relax
  \gdef\@thanks{}\gdef\@author{}\gdef\@title{}
}
\def\@maketitle{%
  \newpage
  \null
  \vskip 1em%
  \begin{center}%
    {\large \@title \par}%
    \vskip 1.5em%
    {\normalsize
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
%    \vskip 1.5em%
%     \@date%
  \end{center}%
  \par
  \vskip 1.5em}

% Abstract
\renewcommand{\abstract}
{\if@twocolumn
   \section*{\abstractname}
 \else 
   \vspace*{-0.5in}
   \begin{center}
     {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}} 
   \end{center}
  \fi
  \if@smallabstract\renewcommand{\baselinestretch}{\smallstretch}\small\else\fi}
{\if@twocolumn\else\relax\fi}

% Redefine captions
\if@smallcaptions
  \newsavebox{\tempbox}
  \renewcommand{\@makecaption}[2]
  {\vspace{10pt}\renewcommand{\baselinestretch}{\smallstretch}
   \small\sbox{\tempbox}{#1: #2}
   \ifthenelse{\lengthtest{\wd\tempbox > \linewidth}}
   { #1: #2\par}{\begin{center}#1: #2\end{center}}}
\fi


% Bibliography
%   make text small and closely separated
\renewcommand{\@biblabel}[1]{[#1]}
\renewenvironment{thebibliography}[1]
     {%\section{\refname
      %  \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      \renewcommand{\baselinestretch}{\smallstretch}\small
      \setlength{\itemsep}{-0.125cm}
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@\LaTeX@warning{Empty `thebibliography' environment}}%
      \endlist}

