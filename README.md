# Mycobiome Project Repository

*Author: Lucy Goodyear*  
*Created: 06/12/19*

This repository contains all relevant files for my MRes project.

**Analysis** contains all my scripts.

**Proposal** contains all documents relating to my proposal.

At present, I am unable to upload my LaTeX file with my thesis as it contains data that is yet to be formally published.

For reference here is my Latex code to match my referencing style to the journal Envnvironmental Microbiology:

```% bibliography set up to match journal
\usepackage[style=authoryear,citestyle=authoryear,dashed=false,doi=false,isbn=false,url=false,eprint=false,maxbibnames=8,minbibnames=6,maxcitenames=2,backend=bibtex]{biblatex}
% seperate name and year by comma in citations
\DeclareDelimFormat{nameyeardelim}{\addcomma\space}
% remove default italics for journal title
\DeclareFieldFormat{journaltitle}{\textbf{#1}}
\DeclareFieldFormat{journaltitle}{#1}
% remove quotations marks from paper title
\DeclareFieldFormat*{title}{\textbf{#1}}
\DeclareFieldFormat{title}{#1}
% set up formatting for "et al."
\DefineBibliographyStrings{english}{%
	andothers = {\em et\addabbrvspace al\adddot}}
% remove "In:" before journal name
\renewbibmacro{in:}{}
\addbibresource{MRes.bib} 
```

## Requirements

All code has been created for Mac so there may be a few differences in commands with respect to Linux. MacTeX or equivalent is required for LaTeX related files.