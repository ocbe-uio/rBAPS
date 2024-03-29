% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/laskeMuutokset12345.R
\name{spatialMixture_muutokset}
\alias{spatialMixture_muutokset}
\title{Calculate changes (spatial mixture class)}
\description{
Calculate changes (spatial mixture class)

Calculate changes (spatial mixture class)
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-spatialMixture_muutokset-laskeMuutokset}{\code{spatialMixture_muutokset$laskeMuutokset()}}
\item \href{#method-spatialMixture_muutokset-laskeMuutokset2}{\code{spatialMixture_muutokset$laskeMuutokset2()}}
\item \href{#method-spatialMixture_muutokset-laskeMuutokset3}{\code{spatialMixture_muutokset$laskeMuutokset3()}}
\item \href{#method-spatialMixture_muutokset-laskeMuutokset5}{\code{spatialMixture_muutokset$laskeMuutokset5()}}
\item \href{#method-spatialMixture_muutokset-clone}{\code{spatialMixture_muutokset$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialMixture_muutokset-laskeMuutokset"></a>}}
\if{latex}{\out{\hypertarget{method-spatialMixture_muutokset-laskeMuutokset}{}}}
\subsection{Method \code{laskeMuutokset()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialMixture_muutokset$laskeMuutokset(
  ind,
  rowsFromInd,
  data,
  adjprior,
  priorTerm,
  logml,
  cliques,
  separators
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{ind}}{ind}

\item{\code{rowsFromInd}}{rowsFromInd}

\item{\code{data}}{data}

\item{\code{adjprior}}{adjprior}

\item{\code{priorTerm}}{priorTerm}

\item{\code{logml}}{logml}

\item{\code{cliques}}{cliques}

\item{\code{separators}}{separators}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialMixture_muutokset-laskeMuutokset2"></a>}}
\if{latex}{\out{\hypertarget{method-spatialMixture_muutokset-laskeMuutokset2}{}}}
\subsection{Method \code{laskeMuutokset2()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialMixture_muutokset$laskeMuutokset2(
  i1,
  rowsFromInd,
  data,
  adjprior,
  priorTerm,
  logml,
  cliques,
  separators
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{i1}}{i1}

\item{\code{rowsFromInd}}{rowsFromInd}

\item{\code{data}}{data}

\item{\code{adjprior}}{adjprior}

\item{\code{priorTerm}}{priorTerm}

\item{\code{logml}}{logml}

\item{\code{cliques}}{cliques}

\item{\code{separators}}{separators}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialMixture_muutokset-laskeMuutokset3"></a>}}
\if{latex}{\out{\hypertarget{method-spatialMixture_muutokset-laskeMuutokset3}{}}}
\subsection{Method \code{laskeMuutokset3()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialMixture_muutokset$laskeMuutokset3(
  T2,
  inds2,
  rowsFromInd,
  data,
  adjprior,
  priorTerm,
  i1,
  logml,
  cliques,
  separators
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{T2}}{T2}

\item{\code{inds2}}{inds2}

\item{\code{rowsFromInd}}{rowsFromInd}

\item{\code{data}}{data}

\item{\code{adjprior}}{adjprior}

\item{\code{priorTerm}}{priorTerm}

\item{\code{i1}}{i1}

\item{\code{logml}}{logml}

\item{\code{cliques}}{cliques}

\item{\code{separators}}{separators}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialMixture_muutokset-laskeMuutokset5"></a>}}
\if{latex}{\out{\hypertarget{method-spatialMixture_muutokset-laskeMuutokset5}{}}}
\subsection{Method \code{laskeMuutokset5()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialMixture_muutokset$laskeMuutokset5(
  inds,
  rowsFromInd,
  data,
  adjprior,
  priorTerm,
  logml,
  cliques,
  separators,
  i1,
  i2
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{inds}}{inds}

\item{\code{rowsFromInd}}{rowsFromInd}

\item{\code{data}}{data}

\item{\code{adjprior}}{adjprior}

\item{\code{priorTerm}}{priorTerm}

\item{\code{logml}}{logml}

\item{\code{cliques}}{cliques}

\item{\code{separators}}{separators}

\item{\code{i1}}{i1}

\item{\code{i2}}{i2}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialMixture_muutokset-clone"></a>}}
\if{latex}{\out{\hypertarget{method-spatialMixture_muutokset-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialMixture_muutokset$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
